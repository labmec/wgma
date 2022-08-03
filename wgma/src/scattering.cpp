#include "scattering.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "cmeshtools_impl.hpp"

#include <TPZSpStructMatrix.h>
#include <pzstepsolver.h>
#include <Electromagnetics/TPZPlanarWgScatt.h>
#include <Electromagnetics/TPZPlanarWgScattSrc.h>
#include <TPZSimpleTimer.h>
#include <pzcompelwithmem.h>
#include <pzaxestools.h>


#include <cassert>


namespace wgma::scattering{
  Analysis::Analysis(TPZAutoPointer<TPZCompMesh> mesh,
                     const int n_threads,
                     const bool reorder_eqs,
                     const bool filter_bound) :
    TPZLinearAnalysis(mesh,reorder_eqs),
    m_filter_bound(filter_bound){
    
    m_cmesh = mesh;

    TPZAutoPointer<TPZStructMatrix> strmtrx =
      new TPZSpStructMatrix<CSTATE>(m_cmesh);

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
    }
    SetStructuralMatrix(strmtrx);

    ///Setting a direct solver
    TPZStepSolver<CSTATE> step;
    step.SetDirect(ELU);
    SetSolver(step);
  }


  void Analysis::Assemble(){
    TPZSimpleTimer assemble("Assemble");
    //assembles the system
    TPZLinearAnalysis::Assemble();    
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
    ///solves the system
    TPZLinearAnalysis::Solve();
  }
  void Analysis::Run(){
    Assemble();
    Solve();
  }


  void LoadPrescribedSource(TPZPlanarWgScattSrc *mat,
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

  template<int CELDIM>
  void LoadModalAnalysisSource(TPZPlanarWgScattSrc *mat,
                               const TPZVec<int64_t> &mem_indices,
                               const TPZIntPoints &intrule,
                               TPZInterpolationSpace *cel,
                               TPZTransform<> t)
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
        ptsol.sol = data.sol[0][0];
        //let us take the derivatives to the global coordinates
        if constexpr(CELDIM==2){
          TPZFNMatrix<3,CSTATE> dsol(CELDIM, 1, 0.);
          for(auto ix = 0; ix < CELDIM; ix++){dsol(ix,0) = data.dsol[0].GetVal(ix,0);}
          TPZFNMatrix<3,CSTATE> dsolx(3, 1, 0.);
      
          TPZAxesTools<CSTATE>::Axes2XYZ(dsol, dsolx, data.axes);

          CSTATE im {0,1};
          for(auto ix = 0; ix < 3; ix++){ptsol.dsol[ix] = dsolx.GetVal(ix,0);}
        }//otherwise the normal derivative is zero
      
        (*(mat->GetMemory()))[mem_indices[ipt]] = ptsol;
      }
  }
  
  void LoadSource(
    TPZAutoPointer<TPZCompMesh> scatt_cmesh,
    std::variant<
    wgma::scattering::Source1D,
    wgma::scattering::SourceWgma> source
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
        LoadPrescribedSource(mat, mem_indices, intrule,
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
          LoadModalAnalysisSource<1>(mat, mem_indices, intrule, dcel, dummy_t);
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

          LoadModalAnalysisSource<2>(mat,mem_indices,intrule,neigh_cel,t1);
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
  CMeshScattering2D(TPZAutoPointer<TPZGeoMesh> gmesh,
                    const wgma::planarwg::mode mode, int pOrder,
                    wgma::cmeshtools::PhysicalData &data,
                    const std::set<int> src_id_set,
                    const STATE lambda, const REAL scale)
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
    //we keep track of PML neighbours because we might need it for sources
    std::map<int,int> pml_neighs;
    
    TPZPlanarWgScatt::ModeType matmode;
    switch(mode){
    case wgma::planarwg::mode::TE:
      matmode = TPZPlanarWgScatt::ModeType::TE;
      break;
    case wgma::planarwg::mode::TM:
      matmode = TPZPlanarWgScatt::ModeType::TM;
      break;
    }
    for(auto [id,er,ur] : data.matinfovec){
      auto *mat = new TPZPlanarWgScatt(id,er,ur,lambda,matmode,scale);
      scatt_cmesh->InsertMaterialObject(mat);
      //for pml
      realvolmats.insert(id);
      volmats.insert(id);
      allmats.insert(id);
    }
  
    for(auto pml : data.pmlvec){
      const auto neighs = wgma::cmeshtools::AddRectangularPMLRegion<TPZPlanarWgScatt>
        (pml, realvolmats, gmesh, scatt_cmesh);
      for(auto id : pml.ids){
        volmats.insert(id);
        allmats.insert(id);
        pml_neighs[id] = neighs.at(id);
      }
    }

  
    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);

    /**let us associate each boundary with a given material.
       this is important for the source boundary*/
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
      allmats.insert(id);
    }

    scatt_cmesh->SetAllCreateFunctionsContinuous();
    scatt_cmesh->SetDefaultOrder(pOrder);
    scatt_cmesh->AutoBuild(allmats);
    scatt_cmesh->CleanUpUnconnectedNodes();


    //we insert all the materials in the computational mesh
    for(auto src_id : src_id_set){
      auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, src_id, volmats);
      if(!res.has_value()){
        PZError<<__PRETTY_FUNCTION__
               <<"\n could not find a material adjacent to the source.\nAborting..."
               <<std::endl;
        DebugStop();
      }
      
      const auto src_volid = [res,pml_neighs](){
        if(pml_neighs.count(res.value())){
          //pmlmat
          return pml_neighs.at(res.value());
        }else{
          //volmat
          return res.value();
        }
      }();
      
      for(auto [id,er,ur] : data.matinfovec){
        if(id == src_volid){
          auto srcMat = new TPZPlanarWgScattSrc(src_id,er,ur,lambda,matmode,scale);
          scatt_cmesh->InsertMaterialObject(srcMat);
        }
      }
    }
    //now we insert the proper material
    scatt_cmesh->SetAllCreateFunctionsContinuousWithMem();
    //we want different memory areas for each integration point
    gSinglePointMemory = false;
    //create computational elements with memory for the source
    scatt_cmesh->AutoBuild(src_id_set);

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
    for(auto [_,mat] : cmesh->MaterialVec()){
      auto scatt_mat =
        dynamic_cast<TPZPlanarWgScattSrc*>(mat);
      if(scatt_mat){
        scatt_mat->SetBeta(beta);
      }
    }
  
  }
};