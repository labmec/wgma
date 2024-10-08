#include "scattering.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "cmeshtools_impl.hpp"

#include <TPZSpStructMatrix.h>
#include <TPZSSpStructMatrix.h>
#include <TPZStructMatrixOMPorTBB.h>
#include <TPZCutHillMcKee.h>

#include <pzstepsolver.h>
#include <Electromagnetics/TPZPlanarWgScatt.h>
#include <Electromagnetics/TPZPlanarWgScattSrc.h>
#include <Electromagnetics/TPZScattering.h>
#include <Electromagnetics/TPZScatteringSrc.h>
#include <Electromagnetics/TPZWgma.h>
#include <TPZNullMaterial.h>
#include <TPZSimpleTimer.h>
#include <pzcompelwithmem.h>
#include <pzaxestools.h>
#include <TPZPardisoSolver.h>
#include <tpzsparseblockdiagonal.h>
#include <TPZCompMeshTools.h>
#include <cassert>



namespace wgma::scattering{

  bool using_tbb_mat{false};
  
  Analysis::Analysis(TPZAutoPointer<TPZCompMesh> mesh,
                     const int n_threads,
                     const bool reorder_eqs,
                     const bool filter_bound,
                     const bool is_sym) :
    TPZLinearAnalysis(),
    m_filter_bound(filter_bound), m_sym(is_sym){

#ifdef PZ_USING_METIS
    const auto renumtype = RenumType::EMetis;
#else
    const auto renumtype = RenumType::ECutHillMcKee;
#endif
    this->CreateRenumberObject(renumtype);
    this->SetCompMesh(mesh.operator->(), reorder_eqs);
    m_cmesh = mesh;

    TPZAutoPointer<TPZStructMatrix> strmtrx = nullptr;

    if(m_sym){
      if(using_tbb_mat){
        auto mtrx = new TPZSSpStructMatrix<CSTATE,TPZStructMatrixOMPorTBB<CSTATE>>(m_cmesh);
        mtrx->SetTBBorOMP(true);
        mtrx->SetShouldColor(false);
        strmtrx = mtrx;
      }else{
        auto mtrx = new TPZSSpStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>(m_cmesh);
        strmtrx = mtrx;
      }
    }else{
      if(using_tbb_mat){
        auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOMPorTBB<CSTATE>>(m_cmesh);
        mtrx->SetTBBorOMP(true);
        mtrx->SetShouldColor(false);
        strmtrx = mtrx;
      }else{
        auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>(m_cmesh);
        strmtrx = mtrx;
      }
    }

    strmtrx->SetNumThreads(n_threads);
    
    TPZVec<int64_t> activeEquations;
  
    
    
    m_n_dofs = m_cmesh->NEquations();
    if(filter_bound){
      int n_dofs_before = m_n_dofs;
      std::set<int64_t> boundConnects;
      wgma::cmeshtools::FilterBoundaryEquations(m_cmesh, activeEquations,
                                                boundConnects);
      m_n_dofs = activeEquations.size();
      std::cout<<"neq(before): "<<n_dofs_before
               <<"\tneq(after): "<<m_n_dofs<<std::endl;
      strmtrx->EquationFilter().SetActiveEquations(activeEquations);
    }else{
      std::cout<<"neq: "<<m_n_dofs<<std::endl;;
    }
    SetStructuralMatrix(strmtrx);

    ///Setting a direct solver
    TPZStepSolver<CSTATE> step;
    if(m_sym){
      step.SetDirect(ELDLt);
    }else{
      step.SetDirect(ELU);
    }
    SetSolver(step);
  }


  void Analysis::Assemble(){
    TPZSimpleTimer assemble("Assemble");
    //assembles the system
    TPZLinearAnalysis::Assemble();
    ///by default, matrices are set as symmetric
    //user can change it before calling solve if needed
    //(i.e., lossless systems can result in hermitian mats)
    if(m_sym){
      auto mat = this->GetSolver().Matrix();
      mat->SetSymmetry(SymProp::Sym);
      mat->SetDefPositive(false);
    }
  }
  
  void Analysis::AssembleRhs(std::set<int> matids){
    auto strmat = this->StructMatrix();
    auto matids_cp = strmat->MaterialIds();
    strmat->SetMaterialIds(matids);
    TPZLinearAnalysis::AssembleResidual();
    strmat->SetMaterialIds(matids_cp);
  }
  
  void Analysis::Solve(){
    TPZSimpleTimer solve("Solve");
    TPZLinearAnalysis::Solve();
  }
  void Analysis::Run(){
    Assemble();
    Solve();
  }


  int ColorEqGraph(const TPZVec<int64_t> &graph, const TPZVec<int64_t> &graphindex,
                   const TPZMatrix<CSTATE> &mat,
                   const int64_t neqs, TPZVec<int> &colors)
  {
  
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
      int64_t ibl;
      //we iterate over the blocks of the input graph
      for(ibl=0; ibl<nblocks; ibl++)
      {
        if(colors[ibl] != -1) continue;
        const int64_t first = graphindex[ibl];
        const int64_t last = graphindex[ibl+1];
        bool is_free = true;
        for(auto ieq=first; ieq<last; ieq++)    
        {
          const auto roweq = graph[ieq];
          TPZManVector<int64_t, 300> indices;
          mat.GetRowIndices(roweq,indices);
          for(auto ieq : indices){
            if(eqcolor[ieq] == color){is_free = false;}
          }
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
  
  TPZAutoPointer<TPZMatrixSolver<CSTATE>>
  Analysis::BuildBlockPrecond(const TPZVec<int64_t> &eqgraph,
                              const TPZVec<int64_t> &graphindex,
                              const bool overlap)
  {
    TPZSimpleTimer timer("BuildBlockPrecond");
    const int64_t neq = this->StructMatrix()->EquationFilter().NActiveEquations();
    auto mySolver = dynamic_cast<TPZMatrixSolver<CSTATE>*>(this->Solver());
    std::cout<<"Building "<<graphindex.size()-1<<" blocks"<<std::endl;
    if(overlap){
      auto sp = new TPZSparseBlockDiagonal<CSTATE>(eqgraph, graphindex, neq);
      TPZStepSolver<CSTATE> *step = new TPZStepSolver<CSTATE>(sp);
      step->SetDirect(ELU);
      //this will allow the TPZAnalysis::Solve to call UpdateFrom and insert values
      step->SetReferenceMatrix(mySolver->Matrix());
      return step;
    }else{
      TPZVec<int> colors(neq,0);
      const int numcolors =
        ColorEqGraph(eqgraph, graphindex, mySolver->Matrix(), neq, colors);
      std::cout<<"Blocks divided into "<<numcolors<<" colors"<<std::endl;
      return this->BuildSequenceSolver<CSTATE>(eqgraph, graphindex, neq, numcolors, colors);
    }
  }


  void LoadPrescribedSource1D(TPZPlanarWgScattSrc *mat,
                            const TPZVec<int64_t> &mem_indices,
                            const TPZIntPoints &intrule,
                            TPZCompEl *cel,
                            wgma::scattering::Source1D src
                            )
  {

    auto *gel = cel->Reference();
    //number of integration points
    const auto npts = mem_indices.size();
    //integration point (reference element)
    TPZManVector<REAL,3> pt_qsi(intrule.Dimension());
    TPZManVector<REAL,3> pt_x(3);
    REAL w{0};//wont be used
    for(auto ipt = 0; ipt < npts; ipt++){
        TPZScatteredSol2D ptsol;
        intrule.Point(ipt, pt_qsi, w);
        gel->X(pt_qsi, pt_x);
        src.func(pt_x, ptsol.sol, ptsol.dsol);
        ptsol.x = pt_x;
        (*(mat->GetMemory()))[mem_indices[ipt]] = ptsol;
      }
  }

  void LoadPrescribedSource2D(TPZScatteringSrc *mat,
                            const TPZVec<int64_t> &mem_indices,
                            const TPZIntPoints &intrule,
                            TPZCompEl *cel,
                            wgma::scattering::Source2D src
                            )
  {

    PZError<<__PRETTY_FUNCTION__
           <<"\nNot yet implemented! Must also take grad_ez!"<<std::endl;
    DebugStop();
    auto *gel = cel->Reference();
    //number of integration points
    const auto npts = mem_indices.size();
    //integration point (reference element)
    TPZManVector<REAL,3> pt_qsi(intrule.Dimension());
    TPZManVector<REAL,3> pt_x(3);
    REAL w{0};//wont be used
    for(auto ipt = 0; ipt < npts; ipt++){
        TPZScatteredSol3D ptsol;
        intrule.Point(ipt, pt_qsi, w);
        gel->X(pt_qsi, pt_x);
        src.func(pt_x, ptsol.et);
        ptsol.x = pt_x;
        (*(mat->GetMemory()))[mem_indices[ipt]] = ptsol;
      }
  }

  template<int CELDIM>
  void LoadModalAnalysisSource(TPZPlanarWgScattSrc *mat,
                               const TPZVec<int64_t> &mem_indices,
                               const TPZIntPoints &intrule,
                               TPZInterpolationSpace *cel,
                               TPZTransform<> t,
                               CSTATE coeff)
  {

    //number of integration points
    const auto npts = mem_indices.size();
    /*
      the integration rule will be from a 1d element
      the computational element might be the 1d
      el itself (1d modal analysis) or a 2d neighbour(2d periodic ma)
     */

    
    //integration point (reference element)
    TPZManVector<REAL,3> pt_qsi1(intrule.Dimension());
    TPZManVector<REAL,CELDIM> pt_qsi2(CELDIM);
    TPZManVector<REAL,3> pt_x(3);
    REAL w{0};//wont be used

    TPZMaterialDataT<CSTATE> data;
    cel->InitMaterialData(data);
    for(auto ipt = 0; ipt < npts; ipt++){
        TPZScatteredSol2D ptsol;
        intrule.Point(ipt, pt_qsi1, w);
        t.Apply(pt_qsi1, pt_qsi2);

        data.fNeedsSol = true;
        cel->ComputeRequiredData(data, pt_qsi2);
        //x (for debugging)
        ptsol.x = data.x;
        //sol
        ptsol.sol = coeff*data.sol[0][0];
        //let us take the derivatives to the global coordinates
        if constexpr(CELDIM==2){
          TPZFNMatrix<3,CSTATE> dsol(CELDIM, 1, 0.);
          for(auto ix = 0; ix < CELDIM; ix++){dsol(ix,0) = data.dsol[0].GetVal(ix,0);}
          TPZFNMatrix<3,CSTATE> dsolx(3, 1, 0.);
      
          TPZAxesTools<CSTATE>::Axes2XYZ(dsol, dsolx, data.axes);

          CSTATE im {0,1};
          for(auto ix = 0; ix < 3; ix++){ptsol.dsol[ix] = coeff*dsolx.GetVal(ix,0);}
        }//otherwise the normal derivative is zero
      
        (*(mat->GetMemory()))[mem_indices[ipt]] = ptsol;
      }
  }

  void LoadModalAnalysisSource2D(TPZScatteringSrc *mat,
                                 const TPZVec<int64_t> &mem_indices,
                                 const TPZIntPoints &intrule,
                                 TPZMultiphysicsElement *mfcel,
                                 CSTATE coeff)
  {

    //number of integration points
    const auto npts = mem_indices.size();

    
    constexpr int celdim {2};
    //integration point (reference element)
    TPZManVector<REAL,3> pt_qsi(intrule.Dimension());
    TPZManVector<REAL,3> pt_x(3);
    REAL w{0};//wont be used

    TPZManVector<TPZMaterialDataT<CSTATE>,2> datavec(2,{});
    mfcel->InitMaterialData(datavec);

    for(auto ipt = 0; ipt < npts; ipt++){
        TPZScatteredSol3D ptsol;
        intrule.Point(ipt, pt_qsi, w);

        for (auto &data : datavec){
          data.fNeedsSol = true;
        }
        TPZManVector<TPZTransform<> > trvec;
        mfcel->AffineTransform(trvec);
        mfcel->ComputeRequiredData(pt_qsi,trvec, datavec);
        //x (for debugging)
        ptsol.x = datavec[0].x;
        //sol
        //TODO: compute sol
        // ptsol.sol = data.sol[0][0];
        /*
          we are only interested in the tangential component of the solution
         */
        const auto &hcurl_sol = datavec[TPZWgma::HCurlIndex()].sol[0];
        const auto &h1_dsol = datavec[TPZWgma::H1Index()].dsol[0];

        TPZFNMatrix<3,CSTATE> dsol(celdim, 1, 0.);
        for(auto ix = 0; ix < celdim; ix++){dsol(ix,0) = h1_dsol.GetVal(ix,0);}
        TPZFNMatrix<3,CSTATE> h1_grad_sol(3, 1, 0.);
      
        TPZAxesTools<CSTATE>::Axes2XYZ(h1_dsol, h1_grad_sol, datavec[0].axes);
        
        ptsol.et[0] = coeff*hcurl_sol[0];
        ptsol.et[1] = coeff*hcurl_sol[1];
        ptsol.grad_ez[0] = -coeff*h1_grad_sol.GetVal(0,0);
        ptsol.grad_ez[1] = -coeff*h1_grad_sol.GetVal(1,0);
        (*(mat->GetMemory()))[mem_indices[ipt]] = ptsol;
      }
  }
  
  void LoadSource1D(
    TPZAutoPointer<TPZCompMesh> scatt_cmesh,
    std::variant<
    wgma::scattering::Source1D,
    wgma::scattering::SourceWgma> source,
    CSTATE coeff
    )
  {
    const bool prescribed_source =
      std::holds_alternative<wgma::scattering::Source1D>(source);
    const std::set<int> src_id_set = [prescribed_source, &source](){
      if(prescribed_source){
        return std::get<wgma::scattering::Source1D>(source).id;
      }else{
        return std::get<wgma::scattering::SourceWgma>(source).id;
      }
    }();
    
    auto * gmesh = scatt_cmesh->Reference();
    //this will make it easier for us to find the correct elements from modal analysis
    if(!prescribed_source){
      gmesh->ResetReference();
      //the geometric mesh will point to the modal analysis mesh
      std::get<wgma::scattering::SourceWgma>(source).modal_cmesh->LoadReferences();
    }
  
    const int modal_dim = [prescribed_source, &source](){
      if(prescribed_source){ return 0;}
      else{
        auto modal_cmesh =
          std::get<wgma::scattering::SourceWgma>(source).modal_cmesh;
        return modal_cmesh->Dimension();
      }
    }();
    
    TPZVec<int64_t> mem_indices;
    for(auto cel : scatt_cmesh->ElementVec()){
      //check if material is source mat
      const auto has_src = src_id_set.count(cel->Material()->Id());
      if(!has_src){continue;}
      
      auto mat =
        dynamic_cast<TPZPlanarWgScattSrc *>(cel->Material());
      
      assert(mat);

      //get integration rule
      auto &intrule = cel->GetIntegrationRule();
      //get memory indices (indices for all integration points)
      cel->GetMemoryIndices(mem_indices);
      //just to make sure
      assert(intrule.NPoints() == mem_indices.size());

      if(prescribed_source){
        auto my_src = std::get<wgma::scattering::Source1D>(source);
        LoadPrescribedSource1D(mat, mem_indices, intrule,
                               cel,my_src);
      }else{
        //1D geometric element
          auto gel = cel->Reference();
        if(modal_dim == 1){
          //1d modal analysis
          TPZTransform<> dummy_t(1);
          auto modal_cel = gel->Reference();
          auto dcel =
            dynamic_cast<TPZInterpolationSpace*>(modal_cel);
          assert(dcel);
          LoadModalAnalysisSource<1>(mat, mem_indices, intrule, dcel, dummy_t,coeff);
        }else{
          //periodic modal analysis
          
          //1D element of modal analysis mesh
          TPZGeoElSide gelside(gel,gel->NSides()-1);

          auto neigh_gelside = gelside.Neighbour();
          while(neigh_gelside!=gelside){
            auto gel_2d = neigh_gelside.Element();
            if(gel_2d->Reference()){//we found the correct element
              break;
            }
            neigh_gelside = neigh_gelside.Neighbour();
          }
          auto neigh_gel = neigh_gelside.Element();
          if(neigh_gel->Dimension() != 2){
            DebugStop();//why we couldnt find the 2D element?
          }
          //transforms from 1D element to 2D element's side
          TPZTransform<> t1(1);
          gelside.SideTransform3(neigh_gelside,t1);
    
          //transforms from 2D element's side to 2D element
          auto t2 =
            neigh_gel->SideToSideTransform(neigh_gelside.Side(), neigh_gel->NSides()-1);
    

          t1 = t2.Multiply(t1);

          //2d computational element
          auto neigh_cel =
            dynamic_cast<TPZInterpolationSpace*>(neigh_gel->Reference());
          assert(neigh_cel);

          LoadModalAnalysisSource<2>(mat,mem_indices,intrule,neigh_cel,t1,coeff);
        }
      }
    }

    if(!prescribed_source){
      gmesh->ResetReference();
      //the geometric mesh now points back to the scattering mesh
      scatt_cmesh->LoadReferences();
    }
  }

  void LoadSource2D(
    TPZAutoPointer<TPZCompMesh> scatt_cmesh,
    std::variant<
    wgma::scattering::Source2D,
    wgma::scattering::SourceWgma> source,
    CSTATE coeff
    )
  {
    const bool prescribed_source =
      std::holds_alternative<wgma::scattering::Source2D>(source);
    const std::set<int> src_id_set = [prescribed_source, &source](){
      if(prescribed_source){
        return std::get<wgma::scattering::Source2D>(source).id;
      }else{
        return std::get<wgma::scattering::SourceWgma>(source).id;
      }
    }();
    
    auto * gmesh = scatt_cmesh->Reference();
    //this will make it easier for us to find the correct elements from modal analysis
    if(!prescribed_source){
      gmesh->ResetReference();
      //the geometric mesh will point to the modal analysis mesh
      std::get<wgma::scattering::SourceWgma>(source).modal_cmesh->LoadReferences();
    }
  
    const int modal_dim = [prescribed_source, &source](){
      if(prescribed_source){ return 0;}
      else{
        auto modal_cmesh =
          std::get<wgma::scattering::SourceWgma>(source).modal_cmesh;
        return modal_cmesh->Dimension();
      }
    }();


    if(modal_dim != 0 && modal_dim != 2){
      PZError<<__PRETTY_FUNCTION__
               <<"\nUnexpected dimension of modal analysis mesh! dim: "
             <<modal_dim<<std::endl;
        DebugStop();
    }
    
    TPZVec<int64_t> mem_indices;
    for(auto cel : scatt_cmesh->ElementVec()){
      //check if material is source mat
      const auto has_src = src_id_set.count(cel->Material()->Id());
      if(!has_src){continue;}
      
      auto mat =
        dynamic_cast<TPZScatteringSrc *>(cel->Material());
      
      if(!mat){
        PZError<<__PRETTY_FUNCTION__
               <<"\nUnexpected source material! id: "<<cel->Material()->Id()
               <<std::endl;
        DebugStop();
      }

      //get integration rule
      auto &intrule = cel->GetIntegrationRule();
      //get memory indices (indices for all integration points)
      cel->GetMemoryIndices(mem_indices);
      //just to make sure
      assert(intrule.NPoints() == mem_indices.size());

      if(prescribed_source){
        auto my_src = std::get<wgma::scattering::Source2D>(source);
        LoadPrescribedSource2D(mat, mem_indices, intrule,
                             cel,my_src);
      }else{
        //2D geometric element
        auto gel = cel->Reference();
        auto modal_cel = gel->Reference();
        auto mfcel =
          dynamic_cast<TPZMultiphysicsElement*>(modal_cel);
        if(mfcel){
          LoadModalAnalysisSource2D(mat, mem_indices, intrule, mfcel,coeff);
        }else{
          std::cout<<__PRETTY_FUNCTION__
                   <<"\nCould not obtain multiphysics element!"
                   <<"\nAborting...";
          DebugStop();
        }
      }
    }

    if(!prescribed_source){
      gmesh->ResetReference();
      //the geometric mesh now points back to the scattering mesh
      scatt_cmesh->LoadReferences();
    }
  }
  
  TPZAutoPointer<TPZCompMesh>
  CMeshScattering2DPeriodic(TPZAutoPointer<TPZGeoMesh> gmesh,
                            const wgma::planarwg::mode mode, int pOrder,
                            wgma::cmeshtools::PhysicalData &data,
                            const std::map<int64_t,int64_t> &periodic_els,
                            const std::set<int> src_id_set,
                            const STATE lambda, const REAL scale, bool verbose)
  {
    static constexpr bool isComplex{true};
    static constexpr int dim{2};
    TPZAutoPointer<TPZCompMesh> scatt_cmesh =
      new TPZCompMesh(gmesh,isComplex);
    scatt_cmesh->SetDimModel(dim);

    const int nvolmats = data.matinfovec.size();
    //all mats (so TPZCompMesh::AutoBuild doesnt break on periodic meshes)
    std::set<int> allmats;
    //volumetric mats
    std::set<int> volmats;
    //volumetric mats - pml
    std::set<int> realvolmats;
    
    TPZPlanarWgScatt::ModeType matmode;
    switch(mode){
    case wgma::planarwg::mode::TE:
      matmode = TPZPlanarWgScatt::ModeType::TE;
      break;
    case wgma::planarwg::mode::TM:
      matmode = TPZPlanarWgScatt::ModeType::TM;
      break;
    case wgma::planarwg::mode::Invalid:
      DebugStop();
      break;
    }
    if(verbose){
      std::cout<<"Creating mesh for propagation of  "
               <<wgma::planarwg::mode_to_string(mode)
               <<" modes"<<std::endl;
    }
    if(verbose && data.matinfovec.size()){std::cout<<"VOLMATS:"<<std::endl;}
    for(auto [id,er,ur] : data.matinfovec){
      auto *mat = new TPZPlanarWgScatt(id,er,ur,lambda,matmode,scale);
      scatt_cmesh->InsertMaterialObject(mat);
      if(verbose){
        std::cout<<"\tid "<<id<<" er "<<er<<" ur "<<ur<<std::endl;
      }
      //for pml
      realvolmats.insert(id);
    }

    volmats = realvolmats;
    if(verbose && data.pmlvec.size()){std::cout<<"PMLs:"<<std::endl;}
    for(auto pml : data.pmlvec){
      //skip PMLs of other dimensions
      if(pml->dim != scatt_cmesh->Dimension()){continue;}
      auto cart_pml = TPZAutoPointerDynamicCast<wgma::pml::cart::data>(pml);
      auto cyl_pml = TPZAutoPointerDynamicCast<wgma::pml::cyl::data>(pml);
      if(cart_pml){
        cart_pml->neigh =
          cmeshtools::AddRectangularPMLRegion<TPZPlanarWgScatt>(*cart_pml, realvolmats, gmesh, scatt_cmesh);
        if(verbose){
          std::cout<<"\tid:";
          for(auto [id,neigh]: pml->neigh){std::cout<<' '<<id<<"("<<neigh<<") ";}
          std::cout<<"type  "<<wgma::pml::cart::to_string(cart_pml->t)
                   <<" ax "<<cart_pml->alphax
                   <<" ay "<<cart_pml->alphay
                   <<" az "<<cart_pml->alphaz
                   <<std::endl;
        }
      }else if (cyl_pml){
        cyl_pml->neigh =
          cmeshtools::AddCylindricalPMLRegion<TPZPlanarWgScatt>(*cyl_pml, realvolmats, gmesh, scatt_cmesh);
        if(verbose){
          std::cout<<"\tid(neigh):";
          for(auto [id,neigh]: pml->neigh){std::cout<<' '<<id<<"("<<neigh<<") ";}
          std::cout<<"type  "<<wgma::pml::cyl::to_string(cyl_pml->t)
                   <<" ar "<<cyl_pml->alphar
                   <<" az "<<cyl_pml->alphaz
                   <<std::endl;
        }
      }else{
        DebugStop();
      }
      
      for(auto [id, _] : pml->neigh){
        volmats.insert(id);
      }
    }

    allmats = volmats;

    if(verbose && data.probevec.size()){std::cout<<"PROBES:"<<std::endl;}
    for(auto [id,matdim] : data.probevec){
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterial<CSTATE>(id,matdim,nstate);
      scatt_cmesh->InsertMaterialObject(mat);
      allmats.insert(id);
      if(verbose){
        std::cout<<"\t id "<<id<<" dim "<<matdim<<std::endl;
      }
    }
  
    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);

    /**let us associate each boundary with a given material.
       this is important for the source boundary*/
    if(verbose){std::cout<<"BC:"<<std::endl;}
    for(auto &bc : data.bcvec){
      auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, bc.id, volmats);
      if(!res.has_value()){
        std::cout<<__PRETTY_FUNCTION__
                 <<"\nwarning: could not find neighbour of bc "<<bc.id<<std::endl;
      }
      bc.volid = res.value();
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volmatid = bc.volid;
      auto *volmat =
        dynamic_cast<TPZMaterialT<CSTATE>*> (scatt_cmesh->FindMaterial(volmatid));
      auto *bcmat = volmat->CreateBC(volmat, id, bctype, val1, val2);
      scatt_cmesh->InsertMaterialObject(bcmat);
      if(verbose){
        std::cout<<"\tid "<<id<<" bctype "<<wgma::bc::to_string(bc.t)
                 <<" neigh "<<volmatid<<std::endl;
      }
      allmats.insert(id);
    }

    scatt_cmesh->SetAllCreateFunctionsContinuous();
    scatt_cmesh->SetDefaultOrder(pOrder);
    scatt_cmesh->AutoBuild(allmats);
    scatt_cmesh->CleanUpUnconnectedNodes();

    if(verbose && src_id_set.size()){std::cout<<"SOURCES:"<<std::endl;}
    //we insert all the materials in the computational mesh
    for(auto src_id : src_id_set){
      auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, src_id, realvolmats);
      if(!res.has_value()){
        //if not found, lets include volumetric PML mats in our search
        res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, src_id, volmats);
      }
      if(!res.has_value()){
        PZError<<__PRETTY_FUNCTION__
               <<"\n could not find a material adjacent to the source "<<src_id
               <<".\nAborting..."
               <<std::endl;
        DebugStop();
      }
      

      const auto src_volid = res.value();

      bool found{false};
      //search in volumetric materials
      for(auto [id,er,ur] : data.matinfovec){
        if(id == src_volid){
          auto srcMat = new TPZPlanarWgScattSrc(src_id,er,ur,lambda,matmode,scale);
          scatt_cmesh->InsertMaterialObject(srcMat);
          if(verbose){
            std::cout<<"\tsrc id : "<<src_id<<" neigh "<<src_volid<<std::endl;
          }
          found = true;
          break;
        }
      }

      if(found){continue;}
      //search in pml
      for(auto neigh_pml : data.pmlvec){
        if(neigh_pml->ids.count(src_volid)){
          //found pml region
          const int volid = neigh_pml->neigh[src_volid];
          if(verbose){
            std::cout<<"\tsrc id(pml) : "<<src_id<<" neigh "<<volid<<std::endl;
          }
          for(auto [neigh_id,er,ur] : data.matinfovec){
            if(neigh_id == volid){
              //we found the volumetric neighbour, now we should create the correct PML
              
              //dummy material just so we can take the parameters
              auto src_mat = new TPZPlanarWgScattSrc(src_id,er,ur,lambda,matmode,scale);
              //the src pml will have the same parameters as its pml neigh
              auto src_pml_mat =
                wgma::cmeshtools::ChangeMaterialToPML(src_id,*neigh_pml,src_mat,gmesh);
              delete src_mat;
              scatt_cmesh->InsertMaterialObject(src_pml_mat);
              found = true;
              break;
            }
          }
          if(found){break;}
        }
      }
      //what happened?
      if(!found){
        std::cout<<"source "<<src_id<<" could not be found"<<std::endl;
        DebugStop();
      }
    }
    
    if(src_id_set.size()){
      //now we insert the proper material
      scatt_cmesh->SetAllCreateFunctionsContinuousWithMem();
      //we want different memory areas for each integration point
      gSinglePointMemory = false;
      //create computational elements with memory for the source
      scatt_cmesh->AutoBuild(src_id_set);
    }
    
    if(periodic_els.size()){
      wgma::cmeshtools::SetPeriodic(scatt_cmesh, periodic_els);
    }
    
    return scatt_cmesh;
  }

  TPZAutoPointer<TPZCompMesh>
  CMeshScattering3DPeriodic(TPZAutoPointer<TPZGeoMesh> gmesh,
                            int pOrder,
                            wgma::cmeshtools::PhysicalData &data,
                            const std::map<int64_t,int64_t> &periodic_els,
                            const std::set<int> src_id_set,
                            const STATE lambda, const REAL scale,
                            bool verbose)
  {
    static constexpr bool isComplex{true};
    static constexpr int dim{3};
    TPZAutoPointer<TPZCompMesh> scatt_cmesh = nullptr;
    {
      TPZSimpleTimer tscatt("Init");
      scatt_cmesh = new TPZCompMesh(gmesh,isComplex);
    }
    scatt_cmesh->SetDimModel(dim);

    const int nvolmats = data.matinfovec.size();
    //all mats (so TPZCompMesh::AutoBuild doesnt break on periodic meshes)
    std::set<int> allmats;
    //volumetric mats
    std::set<int> volmats;
    //volumetric mats - pml
    std::set<int> realvolmats;
    //we keep track of PML neighbours because we might need it for sources
    std::map<int,int> pml_neighs;

    {
      TPZSimpleTimer t("Matdata");
      if(verbose && data.matinfovec.size()){std::cout<<"VOLMATS:"<<std::endl;}
      for(auto [id,er,ur] : data.matinfovec){
        auto *mat =  new TPZScattering(id,er,ur,lambda,scale);
        scatt_cmesh->InsertMaterialObject(mat);
        if(verbose){
          std::cout<<"\tid "<<id<<" er "<<er<<" ur "<<ur<<std::endl;
        }
        //for pml
        realvolmats.insert(id);
      }

    //volumetric mats
    volmats = realvolmats;

    {
      TPZSimpleTimer t("PML");
      if(verbose && data.pmlvec.size()){std::cout<<"PMLs:"<<std::endl;}
      for(auto pml : data.pmlvec){
        //skip PMLs of other dimensions
        if(pml->dim != scatt_cmesh->Dimension()){continue;}
        auto cart_pml = TPZAutoPointerDynamicCast<wgma::pml::cart::data>(pml);
        auto cyl_pml = TPZAutoPointerDynamicCast<wgma::pml::cyl::data>(pml);
        if(cart_pml){
          cart_pml->neigh =
            cmeshtools::AddRectangularPMLRegion<TPZScattering>(*cart_pml, realvolmats, gmesh, scatt_cmesh);
          if(verbose){
            std::cout<<"\tid:";
            for(auto [id,neigh]: pml->neigh){std::cout<<' '<<id<<"("<<neigh<<") ";}
            std::cout<<"type  "<<wgma::pml::cart::to_string(cart_pml->t)
                     <<" ax "<<cart_pml->alphax
                     <<" ay "<<cart_pml->alphay
                     <<" az "<<cart_pml->alphaz
                     <<std::endl;
          }
        }else if (cyl_pml){
          cyl_pml->neigh =
            cmeshtools::AddCylindricalPMLRegion<TPZScattering>(*cyl_pml, realvolmats, gmesh, scatt_cmesh);
          if(verbose){
            std::cout<<"\tid(neigh):";
            for(auto [id,neigh]: pml->neigh){std::cout<<' '<<id<<"("<<neigh<<") ";}
            std::cout<<"type  "<<wgma::pml::cyl::to_string(cyl_pml->t)
                     <<" ar "<<cyl_pml->alphar
                     <<" az "<<cyl_pml->alphaz
                     <<std::endl;
          }
        }else{
          DebugStop();
        }
        for(auto [id, _] : pml->neigh){
          volmats.insert(id);
        }
      }
    }

    allmats = volmats;

    if(verbose && data.probevec.size()){std::cout<<"PROBES:"<<std::endl;}
    for(auto [id,matdim] : data.probevec){
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterial<CSTATE>(id,matdim,nstate);
      scatt_cmesh->InsertMaterialObject(mat);
      allmats.insert(id);
      if(verbose){
        std::cout<<"\t id "<<id<<" dim "<<matdim<<std::endl;
      }
    }

    }
    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);

    /**let us associate each boundary with a given material.
       this is important for the source boundary*/
    {
      TPZSimpleTimer tscatt("find bc neigh");
      if(verbose){std::cout<<"BC:"<<std::endl;}
      for(auto &bc : data.bcvec){
        auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, bc.id, volmats);
        if(!res.has_value()){
          std::cout<<__PRETTY_FUNCTION__
                   <<"\nwarning: could not find neighbour of bc "<<bc.id<<std::endl;
        }
        bc.volid = res.value();
        const int bctype = wgma::bc::to_int(bc.t);
        const int id = bc.id;
        const int volmatid = bc.volid;
        auto *volmat =
          dynamic_cast<TPZMaterialT<CSTATE>*> (scatt_cmesh->FindMaterial(volmatid));
        auto *bcmat = volmat->CreateBC(volmat, id, bctype, val1, val2);
        scatt_cmesh->InsertMaterialObject(bcmat);
        if(verbose){
          std::cout<<"\tid "<<id<<" bctype "<<wgma::bc::to_string(bc.t)
                   <<" neigh "<<volmatid<<std::endl;
        }
        allmats.insert(id);
      }

    }
    scatt_cmesh->SetAllCreateFunctionsHCurl();
    scatt_cmesh->SetDefaultOrder(pOrder);

    {
      TPZSimpleTimer tscatt("AutoBuild1");
      scatt_cmesh->AutoBuild(allmats);
      // scatt_cmesh->CleanUpUnconnectedNodes();
    }


    //whether a source material has been created
    std::set<int> src_mats, src_pml_mats;
    //we insert all the materials in the computational mesh
    {
      if(verbose && src_id_set.size()){std::cout<<"SOURCES:"<<std::endl;}
      TPZSimpleTimer tscatt("SrcMats");
      for(auto src_id : src_id_set){
        bool found_this_src = false;
        auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, src_id, realvolmats);
        if(!res.has_value()){
          //if not found, lets include volumetric PML mats in our search
          res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, src_id, volmats);
        }
        if(!res.has_value()){
          PZError<<__PRETTY_FUNCTION__
                 <<"\n could not find a material adjacent to the source.\nAborting..."
                 <<std::endl;
          DebugStop();
        }

        const auto src_neigh_id = res.value();
        
        for(auto [id,er,ur] : data.matinfovec){
          if(id == src_neigh_id){
            auto srcMat = new TPZScatteringSrc(src_id,er,ur,lambda,scale);
            scatt_cmesh->InsertMaterialObject(srcMat);
            if(verbose){
              std::cout<<"\tsrc id : "<<src_id<<" neigh "<<src_neigh_id<<std::endl;
            }
            found_this_src = true;
            src_mats.insert(src_id);
            break;
          }
        }
        if(found_this_src) {continue;}
        //probably a PML material then
        for(auto &pml : data.pmlvec){
          if (pml->ids.count(src_neigh_id)){
            //found neighbouring PML region
            const int volid = pml->neigh[src_neigh_id];
            if(verbose){
              std::cout<<"\tsrc id(pml) : "<<src_id<<" neigh "<<volid<<std::endl;
            }
            for(auto [neigh_id,er,ur] : data.matinfovec){
              if(neigh_id == volid){
                //dummy material just so we can take the parameters
                auto src_mat = new TPZScatteringSrc(src_id,er,ur,lambda,scale);
                //the src pml will have the same parameters as its pml neigh
                auto src_pml_mat =
                  wgma::cmeshtools::ChangeMaterialToPML(src_id,*pml,src_mat,gmesh);
                delete src_mat;
                scatt_cmesh->InsertMaterialObject(src_pml_mat);
                src_pml_mats.insert(src_id);
                found_this_src=true;
                break;
              }
            }
          }
        }
        if(!found_this_src){
          std::cout<<"source "<<src_id<<" could not be found"<<std::endl;
          DebugStop();
        }
      }
    }
    if(src_mats.size() + src_pml_mats.size() > 0){
      TPZSimpleTimer tscatt("AutoBuild2");
      //now we insert the proper material
      scatt_cmesh->SetAllCreateFunctionsHCurlWithMem();
      //we want different memory areas for each integration point
      gSinglePointMemory = false;
      //create computational elements with memory for the source
      scatt_cmesh->AutoBuild(src_id_set);
    }

    if(periodic_els.size()){
      wgma::cmeshtools::SetPeriodic(scatt_cmesh,periodic_els);
    }

    //this will already call cleanup unconnected nodes
    TPZCompMeshTools::CreatedCondensedElements(scatt_cmesh.operator->(), false, false);
    return scatt_cmesh;
  }

  void
  SetPropagationConstant(TPZAutoPointer<TPZCompMesh> cmesh,
                         const CSTATE beta)
  {
    std::streamsize ss = std::cout.precision();
  
    std::cout.precision(std::numeric_limits<STATE>::max_digits10);

    std::cout<<"Setting beta:\n"
             <<'\t'<<beta<<std::endl;
    std::setprecision(ss);

    auto LoadSrc = [](auto mat,auto beta){
      mat->SetBeta(beta);
    };
    
    for(auto [_,mat] : cmesh->MaterialVec()){
      auto scatt_mat_1d =
        dynamic_cast<TPZPlanarWgScattSrc*>(mat);
      if(scatt_mat_1d){
        LoadSrc(scatt_mat_1d,beta);
      }else{
        auto scatt_mat_2d =
          dynamic_cast<TPZScatteringSrc*>(mat);
        if(scatt_mat_2d){
          LoadSrc(scatt_mat_2d,beta);
        }
      }
    }  
  }
};