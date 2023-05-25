#include "precond.hpp"
#include <TPZSimpleTimer.h>
#include <pzcmesh.h>
#include <TPZEquationFilter.h>
#include <tpznodesetcompute.h>
#include <TPZParallelUtils.h>
#include <pzvec_extras.h>
#include <numeric> //for std::accumulate

using namespace wgma::precond;

void wgma::precond::CreateAFWBlocks(TPZAutoPointer<TPZCompMesh> cmesh,
                                    const std::set<int> dirichlet_mats,
                                    const TPZEquationFilter &eqfilt,
                                    TPZVec<int64_t> &eqgraph,
                                    TPZVec<int64_t> &eqgraphindex
                                    )
{
  TPZSimpleTimer timer("CreateAFWBlocks");
  if(!eqfilt.IsActive() && dirichlet_mats.size() > 0){
    DebugStop();
  }
  
  auto gmesh = cmesh->Reference();
  const auto dim = gmesh->Dimension();

  const auto nnodes = gmesh->NNodes();
  /*
    the AFW precond is based on vertex patches, each block has all the equations
    associated with the edges connected to each vertex

    therefore number of blocks = number of vertices

    we need to remove vertices that are associated with dirichlet bc and are filtered out
  */


  /*
    this array will map node indices and block indices
  */
  TPZVec<int64_t> nodemap(nnodes,0);

  const int nel = gmesh->NElements();
  pzutils::ParallelFor(0,nel, [&](int iel){
    const auto el = gmesh->Element(iel);
    if(el && el->Reference() && dirichlet_mats.count(el->MaterialId())){
      TPZManVector<int64_t,20> nodeindices;
      el->GetNodeIndices(nodeindices);
      const auto sz = nodeindices.size();
      for(int i = 0; i < sz; i++){
        nodemap[nodeindices[i]] = -1;
      }
    }
  });

  TPZVec<int64_t> facemap;
  if(dim==3){
    //upper bound for sequence number of face connects
    const auto nindep = cmesh->NIndependentConnects();
    //we will set to 1 all face connects
    TPZVec<int64_t> allcons(nindep,0);
    pzutils::ParallelFor(0, nel, [&](int iel){
      const auto el = gmesh->Element(iel);
      if(el && el->Reference() && el->Dimension() == dim){
        const int ne = el->NSides(1);
        const int nf = el->NSides(2);
        const int ff = el->FirstSide(2);
        for(auto itf = 0; itf < nf; itf++){
          const auto face = ff+itf;
          TPZGeoElSide gelside(el,face);
          auto neigh = gelside.Neighbour();
          bool bcface{false};
          const auto facecon = el->Reference()->Connect(itf+ne);
          if(facecon.IsCondensed() || facecon.LagrangeMultiplier()){continue;}
          const auto seqnum = facecon.SequenceNumber();
          while(neigh!=gelside){
            if(neigh.Element() && neigh.Element()->Dimension() < dim &&
               dirichlet_mats.count(neigh.Element()->MaterialId())){
              bcface=true;
              break;
            }
            neigh++;
          }
          if(!bcface){
#ifdef PZDEBUG
            if(seqnum < 0){
              DebugStop();
            }
#endif
            allcons[seqnum] = 1;
          }
        }
      }
    });
    
    const int n_vol_faces =
    std::accumulate(allcons.begin(),allcons.end(),0,
                    [](const int64_t a, const int64_t b){return a + b;});
    //we allocate the correct size
    facemap.Resize(n_vol_faces);
    int64_t count = 0;
    //now we set the correct sequence numbers
    for(int i = 0; i < nindep; i++){
      if(allcons[i]){facemap[count++] = i;}
    }
  }
  const int nvolfaces = facemap.size();
  
  //the sum of -1 positions is the number of boundary nodes
  const int n_bnd_nodes =
    -1*std::accumulate(nodemap.begin(),nodemap.end(),0,
                       [](const int64_t a, const int64_t b){return a + b;});
  //we store the indices of boundary nodes
  TPZVec<int64_t> bd_nodes(n_bnd_nodes);

  {
    int bdcount = 0;
    int volcount = 0;
    for(int i = 0; i < nnodes; i++){
      if(nodemap[i]<0){
        bd_nodes[bdcount++] = i;
      }
      else{
        nodemap[i] = volcount++;
      }
    }
    if(volcount + bdcount != nnodes){DebugStop();}
  }
  //now, nodemap maps our block indices

  //number of blocks
  const auto nvolnodes = nnodes - bd_nodes.size();
  const auto nbl = dim == 3 ? nvolnodes + nvolfaces : nvolnodes;

  //number of edges associated with each block
  TPZVec<int64_t> bsize(nbl,0);
  //total number of equations in all blocks
  int64_t blockgraphsize = 0;
  //this will keep track of visited edges
  std::set<std::pair<int64_t,int64_t>> edgeset;
  
  for(auto el : gmesh->ElementVec()){
    if(el && el->Dimension() == dim){
      TPZCompEl* cel = el->Reference();
      if(!cel){continue;}
      //first edge
      const auto fe = el->FirstSide(1);
      //number of edges
      const auto ne = el->NSides(1);
      //iterate on element edges
      for(auto ie = 0; ie < ne; ie++){
        //we first check if this edge has been visited already
        const TPZManVector<int64_t,2> vvec =
          {el->SideNodeIndex(fe+ie,0),el->SideNodeIndex(fe+ie,1)};
        const auto pair = vvec[0] > vvec[1] ?
          std::make_pair(vvec[1],vvec[0]) : std::make_pair(vvec[0],vvec[1]);
        if(edgeset.count(pair) == 0){
          edgeset.insert(pair);
          for(int iv = 0; iv < 2; iv++){
            const auto v = vvec[iv];
            const auto blindex = nodemap[v];
            if(blindex > -1){
              bsize[blindex]++;
              blockgraphsize++;
            }
          }
        }
      }
    }
  }
  //we clear the visited edges
  edgeset.clear();

  //face blocks will always have only one connect
  pzutils::ParallelFor(0,nvolfaces,[&bsize,nvolnodes](int itf){
    bsize[nvolnodes+itf] = 1;
  });
  blockgraphsize+=nvolfaces;
  
  
  TPZVec<int64_t> blockgraphindex(nbl+1,0);
  TPZVec<int64_t> blockgraph(blockgraphsize,0);
  for(int ibl = 0; ibl < nbl; ibl++){
    blockgraphindex[ibl+1] = blockgraphindex[ibl] + bsize[ibl];
  }
  //we will use bsize as a counter for current eq of each block
  bsize = 0;
  for(auto el : gmesh->ElementVec()){
    if(el && el->Dimension() == dim){
      TPZCompEl *cel = el->Reference();
      if(!cel){continue;}
      //first edge
      const auto fe = el->FirstSide(1);
      //number of edges
      const auto ne = el->NSides(1);
      //iterate on element edges
      for(auto ie = 0; ie < ne; ie++){
        const auto edge = fe+ie;
        const TPZCompElSide celside(cel, edge);
        const auto cindex = celside.ConnectIndex();
        //gets seqnum of connect
        const auto seqnum = cmesh->ConnectVec()[cindex].SequenceNumber();
        //we first check if this edge has been visited already
        const TPZManVector<int64_t,2> vvec =
          {el->SideNodeIndex(fe+ie,0),el->SideNodeIndex(fe+ie,1)};
        const auto pair = vvec[0] > vvec[1] ?
          std::make_pair(vvec[1],vvec[0]) : std::make_pair(vvec[0],vvec[1]);
        if(edgeset.count(pair)==0){
          edgeset.insert(pair);
          for(int iv = 0; iv < 2; iv++){
            const auto v = vvec[iv];
            const auto blindex = nodemap[v];
            if(blindex > -1){
              //this really shouldnt happen
              if(seqnum < 0){
                DebugStop();
              }
              const auto first = blockgraphindex[blindex];
              const auto pos = first + bsize[blindex]++;
              if(pos >= blockgraphindex[blindex+1]){
                DebugStop();
              }
              blockgraph[pos] = seqnum;
            }
          }
        }
      }
    }
  }

  //now we add the face blocks
  pzutils::ParallelFor(0,nvolfaces,[&](int itf){
    const auto pos = blockgraphindex[nvolnodes+itf];
    blockgraph[pos] = facemap[itf];
  });
  
  

  //we convert the graph from seqnum to eqs
  TPZManVector<int64_t> removed_blocks;
  TPZNodesetCompute::ExpandGraph(blockgraph,blockgraphindex,cmesh->Block(),
                                 eqgraph,eqgraphindex,removed_blocks);
  TPZNodesetCompute::UpdateGraph(blockgraph,blockgraphindex,removed_blocks);
  //now we take filtered equations into account
  if(eqfilt.IsActive()){
    //removed blocks will be resized inside filtergraph, but let us prepare
    //against future changes
    removed_blocks.resize(0);
    TPZNodesetCompute::FilterGraph(eqfilt,eqgraph,eqgraphindex, removed_blocks);
    TPZNodesetCompute::UpdateGraph(blockgraph, blockgraphindex, removed_blocks);
  }
}


int wgma::precond::ColorEqGraph(const TPZVec<int64_t> &graph,
                                const TPZVec<int64_t> &graphindex,
                                const TPZMatrix<CSTATE> &mat,
                                const int64_t neqs, TPZVec<int> &colors)
{
  TPZSimpleTimer timer("precond::ColorEqGraph");
  TPZVec<int> eqcolor(neqs);
  int color = 0;
  bool hasuncolored = true;
  const int64_t nblocks = graphindex.NElements()-1;
  colors.Resize(nblocks);
  colors.Fill(-1);
  while(hasuncolored)
  {
    hasuncolored = false;
    eqcolor.Fill(-1);
    //we iterate over the blocks of the input graph
    for(auto ibl=0; ibl<nblocks; ibl++)
    {
      if(colors[ibl] != -1) continue;
      const int64_t first = graphindex[ibl];
      const int64_t last = graphindex[ibl+1];
      bool is_free = true;
      for(auto ieq=first; ieq<last; ieq++)    
      {
        const auto roweq = graph[ieq];
        if(eqcolor[roweq] == color){is_free = false;break;}
        TPZManVector<int64_t, 300> indices;
        mat.GetRowIndices(roweq,indices);
        for(auto ieq : indices){
          if(eqcolor[ieq] == color){is_free = false;break;}
        }
        if(!is_free){break;}
      }
      if(!is_free)
      {
        hasuncolored = true;
      }
      else
      {
        colors[ibl] = color;
        for(auto ieq=first; ieq<last; ieq++)    
        {
          const auto roweq = graph[ieq];
          TPZManVector<int64_t, 300> indices;
          mat.GetRowIndices(roweq,indices);
          for(auto ieq : indices){
            eqcolor[ieq] = color;
          }
        }
      }
    }
    color++;
  }
  return color;
}


BlockPrecond::BlockPrecond(TPZAutoPointer<TPZMatrix<CSTATE>> refmat,
                           TPZVec<int64_t> &&block,
                           TPZVec<int64_t> &&block_index,
                           const TPZVec<int> &colors,
                           const int nc) :
  m_blockgraph(std::move(block)), m_blockindex(std::move(block_index))
{
  const int nbl = m_blockindex.size()-1;
  for(int ibl = 0; ibl < nbl; ibl++){
    const auto bs = this->BlockSize(ibl);
    if(m_max_bs < bs){m_max_bs = bs;}
  }
  ComputeColorVec(colors,nc);
  SetReferenceMatrix(refmat);
  ComputeInfluence();
}

BlockPrecond::BlockPrecond(TPZAutoPointer<TPZMatrix<CSTATE>> refmat,
                           const TPZVec<int64_t> &block,
                           const TPZVec<int64_t> &block_index,
                           const TPZVec<int> &colors,
                           const int nc) :
  m_blockgraph(block), m_blockindex(block_index)
{
  const int nbl = m_blockindex.size()-1;
  for(int ibl = 0; ibl < nbl; ibl++){
    const auto bs = this->BlockSize(ibl);
    if(m_max_bs < bs){m_max_bs = bs;}
  }
  ComputeColorVec(colors,nc);
  SetReferenceMatrix(refmat);
  ComputeInfluence();
}

BlockPrecond::BlockPrecond(const BlockPrecond &cp) :
  TPZMatrixSolver<CSTATE>(cp),
  m_sym_gs(cp.m_sym_gs),
  m_max_bs(cp.m_max_bs),
  m_blockindex(cp.m_blockindex),
  m_blockgraph(cp.m_blockgraph),
  m_colors(cp.m_colors),
  m_infl(cp.m_infl)
{
  if(cp.m_blockinv.size()){
    const int nbl = m_blockindex.size()-1;
    m_blockinv.Resize(nbl);
    for(int ibl = 0; ibl < nbl; ibl++){
      auto mat = dynamic_cast<TPZFMatrix<CSTATE>*>(cp.m_blockinv[ibl]->Clone());
      m_blockinv[ibl] = mat;
    }
  }
}

BlockPrecond&
BlockPrecond::operator=(const BlockPrecond &cp)
{
  TPZMatrixSolver<CSTATE>::operator=(cp);
  m_sym_gs = cp.m_sym_gs;
  m_max_bs = cp.m_max_bs;
  m_blockindex = cp.m_blockindex;
  m_blockgraph = cp.m_blockgraph;
  m_colors = cp.m_colors;
  m_infl = cp.m_infl;
  if(cp.m_blockinv.size()){
    const int nbl = m_blockindex.size()-1;
    m_blockinv.Resize(nbl);
    for(int ibl = 0; ibl < nbl; ibl++){
      auto mat = dynamic_cast<TPZFMatrix<CSTATE>*>(cp.m_blockinv[ibl]->Clone());
      m_blockinv[ibl] = mat;
    }
  }
  return *this;
}

void
BlockPrecond::UpdateFrom(TPZAutoPointer<TPZBaseMatrix> ref_base)
{
  TPZSimpleTimer timer("precond::BlockPrecond::UpdateFrom");
  auto refmat = TPZAutoPointerDynamicCast<TPZMatrix<CSTATE>>(ref_base);
  if(refmat == this->fReferenceMatrix){
    const int nbl = this->NBlocks();
    if(this->m_blockinv.size()!=nbl){
      this->m_blockinv.Resize(nbl,nullptr);
      for(int ibl = 0; ibl < nbl; ibl++){
        const int bs = BlockSize(ibl);
        m_blockinv[ibl] = new TPZFMatrix<CSTATE>(bs,bs);
      }
    }

    std::atomic<int> blcount{0};
    std::cout<<__PRETTY_FUNCTION__;
    std::cout<<"\nDecomposing blocks...";
    //decompose blocks (coloring need not be taken into account)
    std::mutex mymut;
    pzutils::ParallelFor(0,nbl, [&](int ibl){
      TPZManVector<int64_t,400> indices;
      const int bs = BlockSize(ibl);
      indices.Resize(bs);
      for(int ieq = 0; ieq < bs; ieq++){
        const auto pos = m_blockindex[ibl]+ieq;
        indices[ieq] = m_blockgraph[pos];
      }
      auto &block = *(this->m_blockinv[ibl]);
      refmat->GetSub(indices,block);
      block.Decompose_LU();
      blcount++;
      if(blcount%100==0){
        std::lock_guard lock(mymut);
        std::cout<<"\rcomputed "<<blcount<<" out of "<<nbl<< "blocks"<<std::flush;
      }
    });
    std::cout<<"\rcomputed "<<nbl<<" out of "<<nbl<<" blocks" << std::endl;;
  }
}

void
BlockPrecond::Solve(const TPZFMatrix<CSTATE> &rhs_orig,
                    TPZFMatrix<CSTATE> &du_total,
                    TPZFMatrix<CSTATE> *res)
{
  if(m_blockinv.size() == 0){UpdateFrom(this->fReferenceMatrix);}

  //du and rhs might be same, so we copy rhs
  this->fScratch = rhs_orig;
  const auto &rhs = this->fScratch;
  du_total.Redim(rhs.Rows(),1);
  TPZFMatrix<CSTATE> du = du_total;
  TPZFMatrix<CSTATE> rhscp;
  if(res){
    rhscp = rhs;
  }
  const int nc = this->m_colors.size();

  const auto &refmat = this->fReferenceMatrix;

  
  
  constexpr bool mt{true};
  for(int ic = 0; ic < nc; ic++){
    const int nbl_color = this->m_colors[ic].size();
    if constexpr (mt){
      pzutils::ParallelFor(0,nbl_color, [&](int ibl){
        const auto bl = this->m_colors[ic][ibl];
        SmoothBlock(bl,du,rhs);
      });
    }
    else{
      for(int ibl = 0; ibl < nbl_color; ibl++){
        const auto bl = this->m_colors[ic][ibl];
        SmoothBlock(bl,du,rhs);
      }
    }
    
    //now we compute alpha
    if constexpr (mt){
        pzutils::ParallelFor(0,nbl_color, [&](int ibl){
        const auto bl = this->m_colors[ic][ibl];
        ComputeCorrectionFactor(bl,du,rhs);
      });
    }
    else{
      for(int ibl = 0; ibl < nbl_color; ibl++){
        const auto bl = this->m_colors[ic][ibl];
        ComputeCorrectionFactor(bl,du,rhs);
      }
    }
    refmat->Residual(du,rhs,fScratch);
    du_total += du;
    du.Zero();
  }
  if(m_sym_gs){
    for(int ic = nc-2; ic >= 0; ic--){
      const int nbl_color = this->m_colors[ic].size();
      if constexpr (mt){
        pzutils::ParallelFor(0,nbl_color, [&](int ibl){
          const auto bl = this->m_colors[ic][ibl];
          SmoothBlock(bl,du,rhs);
        });
      }
      else{
        for(int ibl = 0; ibl < nbl_color; ibl++){
          const auto bl = this->m_colors[ic][ibl];
          SmoothBlock(bl,du,rhs);
        }
      }
      
      //now we compute alpha
      if constexpr (mt){
        pzutils::ParallelFor(0,nbl_color, [&](int ibl){
          const auto bl = this->m_colors[ic][ibl];
          ComputeCorrectionFactor(bl,du,rhs);
        });
      }
      else{
        for(int ibl = 0; ibl < nbl_color; ibl++){
          const auto bl = this->m_colors[ic][ibl];
          ComputeCorrectionFactor(bl,du,rhs);
        }
      }
      refmat->Residual(du,rhs,fScratch);
      du_total += du;
      du.Zero();
    }
  }

  if(res){
    refmat->Residual(du_total, rhscp, *res);
  }
}

void
BlockPrecond::SmoothBlock(const int bl, TPZFMatrix<CSTATE> &du,
                          const TPZFMatrix<CSTATE>&rhs)
{
  TPZManVector<int64_t,100> indices;
  TPZFNMatrix<10000,CSTATE> resloc;

  const auto &refmat = this->fReferenceMatrix;
  
  const auto bs = this->BlockSize(bl);
  indices.Resize(bs);
  resloc.Resize(bs, 1);
  for(int ieq = 0; ieq < bs; ieq++){
    const auto pos = m_blockindex[bl]+ieq;
    indices[ieq] = m_blockgraph[pos];
  }
  for (size_t ieq = 0; ieq < bs; ieq++)
  {
    const auto eq = indices[ieq];
    resloc(ieq,0) = rhs.GetVal(eq,0);// - refmat->RowTimesVector(eq, du);
  }
        
  const auto &block_inv = *(this->m_blockinv[bl]);
  block_inv.Substitution(&resloc);
  for (size_t ieq = 0; ieq < bs; ieq++)
  {
    const auto eq = indices[ieq];
    du(eq,0) += resloc.GetVal(ieq,0);
  }
}

void
BlockPrecond::ComputeCorrectionFactor(const int bl, TPZFMatrix<CSTATE> &du,
                                      const TPZFMatrix<CSTATE>&rhs)
{
  TPZManVector<int64_t,100> indices;
  TPZFNMatrix<1000,CSTATE> resloc, kduloc;

  const auto &refmat = this->fReferenceMatrix;
  //block size
  const auto bsz = this->BlockSize(bl);
  //size of influenced eqs
  const auto &infl_eqs = this->m_infl[bl];
  const auto inflsz = infl_eqs.size();
  //total size
  const auto sz = bsz+inflsz;
  indices.Resize(bsz);
  resloc.Resize(sz, 1);
  kduloc.Resize(sz, 1);
  for(int ieq = 0; ieq < bsz; ieq++){
    const auto pos = m_blockindex[bl]+ieq;
    const auto eq = m_blockgraph[pos];
    indices[ieq] = eq;
  }

  for(int ieq = 0; ieq < bsz; ieq++){
    const auto eq = indices[ieq];
    resloc.PutVal(ieq,0,rhs.GetVal(eq,0));
    kduloc.PutVal(ieq,0,refmat->RowTimesVector (eq, du));
  }
  
  for (size_t ieq = 0; ieq < inflsz; ieq++)
  {
    const auto eq = infl_eqs[ieq];
    kduloc.PutVal(bsz+ieq,0,refmat->RowTimesVector (eq, du));
    resloc.PutVal(bsz+ieq,0,rhs.GetVal(eq,0));
  }

  auto normkdu = Norm(kduloc);
  normkdu *= normkdu;
  
  if(normkdu == 0){
    return;
  }
  const auto alpha = Dot(kduloc,resloc).real()/normkdu;
  for (size_t ieq = 0; ieq < bsz; ieq++)
  {
    const auto eq = indices[ieq];
    const auto val = alpha*du.GetVal(eq,0);
    du.PutVal(eq,0,val);
  }
}

void
BlockPrecond::ComputeColorVec(const TPZVec<int> &colors,
                              const int nc)
{
  TPZVec<int64_t> colorcount(nc,0);
  const int nbl = NBlocks();
  
  for(int ibl = 0; ibl < nbl; ibl++){
    colorcount[colors[ibl]]++;
  }
  m_colors.Resize(nc);
  for(int ic = 0; ic < nc; ic++){
    const auto size  = colorcount[ic];
    if(size == 0){
      DebugStop();
    }
    m_colors[ic].Resize(size);
  }
  //now this will keep track of blockcount for each color
  colorcount = 0;
  for(int ibl = 0; ibl < nbl; ibl++){
    const int c = colors[ibl];
    m_colors[c][colorcount[c]++] = ibl;
  }
  //just to make sure we sort all blocks
  for(int ic = 0; ic < nc; ic++){
    std::sort(m_colors[ic].begin(), m_colors[ic].end());
  }
  //#ifdef PZDEBUG
  //sanity check: is the coloring correct?
  for(int ic = 0; ic < nc; ic++){
    if(colorcount[ic]!=m_colors[ic].size()){
      DebugStop();
    }
    std::set<int64_t> eqset;
    int64_t count{0};
    for(auto bl : m_colors[ic]){
      for(auto ieq = m_blockindex[bl]; ieq < m_blockindex[bl+1]; ieq++){
        eqset.insert(m_blockgraph[ieq]);
        count++;
      }
    }
    if(count != eqset.size()){
      DebugStop();
    }
  }
  //#endif
}

void BlockPrecond::ComputeInfluence()
{
  const auto &refmat = this->fReferenceMatrix;
#ifdef PZDEBUG
  if(!refmat){
    DebugStop();
  }
#endif
  const auto nbl = this->NBlocks();
  m_infl.Resize(nbl);
  for(int ibl = 0; ibl < nbl; ibl++){
    TPZManVector<int64_t,400> indices;
    const int bs = BlockSize(ibl);
    indices.Resize(bs);
    TPZManVector<int64_t,800> infl;
    for(int ieq = 0; ieq < bs; ieq++){
      TPZManVector<int64_t, 400> loc_infl;
      const auto pos = m_blockindex[ibl]+ieq;
      const auto eq = m_blockgraph[pos];
      indices[ieq] = eq;
      refmat->GetRowIndices(eq, loc_infl);
      const auto oldsz = infl.size();
      const auto locsz = loc_infl.size();
      const auto newsz = oldsz+locsz;
      infl.Resize(newsz);
      for(auto it = 0; it < locsz; it++){
        infl[oldsz+it] = loc_infl[it];
      }
    }
    RemoveDuplicates(infl);
    const auto sz_orig = infl.size();
    //we must have the memory already allocated
    m_infl[ibl].Resize(sz_orig);
    const auto pos = std::set_difference(infl.begin(),infl.end(),indices.begin(),indices.end(),m_infl[ibl].begin());
    const auto newsz = pos - m_infl[ibl].begin();
    m_infl[ibl].Resize(newsz);
  }

  // //sanity check: is the coloring correct?
  // //this is SLOW
  // const auto nc = m_colors.size();
  // TPZManVector<int64_t,1000> intervec;
  // for(int ic = 0; ic < nc; ic++){
  //   for(auto ibl : m_colors[ic]){
  //     const auto & infl_eqs = m_infl[ibl];
  //     for(auto jbl : m_colors[ic]){
  //       if(ibl == jbl){continue;}
  //       const auto first = m_blockindex[jbl];
  //       const auto sz = this->BlockSize(jbl);
  //       const auto pos =
  //         std::set_intersection(&m_blockgraph[first], &m_blockgraph[first]+sz,
  //                               infl_eqs.begin(),infl_eqs.end(),
  //                               intervec.begin());
  //       const auto inter_sz = pos - intervec.begin();
  //       if(inter_sz!=0){DebugStop();}
  //     }
  //   }
  // }
}