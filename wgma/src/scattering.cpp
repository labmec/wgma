#include "scattering.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "cmeshtools_impl.hpp"

#include <TPZSpStructMatrix.h>
#include <pzstepsolver.h>
#include <Electromagnetics/TPZPlanarWGScattering.h>
#include <Electromagnetics/TPZPlanarWGScatteringSrc.h>
#include <TPZSimpleTimer.h>
#include <pzcompelwithmem.h>
#include <pzaxestools.h>


#include <cassert>


namespace wgma::scattering{
  Analysis::Analysis(TPZAutoPointer<TPZCompMesh> mesh,
                     const int n_threads,
                     const bool reorder_eqs,
                     const bool filter_bound) :
    m_filter_bound(filter_bound){
    
    m_cmesh = mesh;
  
    m_an = new TPZLinearAnalysis(m_cmesh, reorder_eqs);

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
    m_an->SetStructuralMatrix(strmtrx);

    ///Setting a direct solver
    TPZStepSolver<CSTATE> step;
    step.SetDirect(ELU);
    m_an->SetSolver(step);
  }
    
  void Analysis::SetSolver(const TPZMatrixSolver<CSTATE> & solv){
    m_an->SetSolver(solv);
  }
    
  TPZMatrixSolver<CSTATE> & Analysis::GetSolver() const{
    return m_an->MatrixSolver<CSTATE>();
  }

  void Analysis::Run(){
    TPZSimpleTimer total("Total");
    {
      TPZSimpleTimer assemble("Assemble");
      //assembles the system
      m_an->Assemble();
    }
    {
      TPZSimpleTimer solve("Solve");
      ///solves the system
      m_an->Solve();
    }
  }

  void Analysis::PostProcess(std::string filename,
                               const int vtk_res){


    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Field_re");
    scalnames.Push("Field_abs");
    vecnames.Push("Deriv_re");
    vecnames.Push("Deriv_abs");
    const std::string plotfile = filename+".vtk";
    constexpr int dim{2};
    m_an->DefineGraphMesh(dim, scalnames, vecnames,plotfile);
    m_an->PostProcess(vtk_res);
    std::cout<<"\nFinished post processing"<<std::endl;
    std::cout<<std::endl;
  }


  void LoadPrescribedSource(TPZPlanarWGScatteringSrc *mat,
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
  void LoadModalAnalysisSource(TPZPlanarWGScatteringSrc *mat,
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
          for(auto ix = 0; ix < CELDIM; ix++){dsol(ix,0) = data.dsol[0][ix];}
          TPZFNMatrix<3,CSTATE> dsolx(3, 1, 0.);
      
          TPZAxesTools<CSTATE>::Axes2XYZ(dsol, dsolx, data.axes);


          for(auto ix = 0; ix < 3; ix++){ptsol.dsol[ix] = dsolx(ix,0);}
        }//otherwise the normal derivative is zero
      
        (*(mat->GetMemory()))[mem_indices[ipt]] = ptsol;
      }
  }
  
  void LoadSources(
    TPZAutoPointer<TPZCompMesh> scatt_cmesh,
    std::variant<
    wgma::scattering::Source1D,
    wgma::scattering::SourceWgma> source,
    const bool prescribed_source
    )
  {
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
        dynamic_cast<TPZPlanarWGScatteringSrc *>(cel->Material());
      
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
        if(modal_dim == 1){
          //1d modal analysis
          TPZTransform<> dummy_t(1);
          auto dcel =
            dynamic_cast<TPZInterpolationSpace*>(cel);
          assert(dcel);
          LoadModalAnalysisSource<1>(mat, mem_indices, intrule, dcel, dummy_t);
        }else{
          //periodic modal analysis
          //1D geometric element
          auto gel = cel->Reference();
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
  }
  
  TPZAutoPointer<TPZCompMesh>
  CMeshScattering2D(TPZAutoPointer<TPZGeoMesh> gmesh,
                    const wgma::planarwg::mode mode, int pOrder,
                    wgma::cmeshtools::PhysicalData &data,
                    std::variant<
                    wgma::scattering::Source1D,
                    wgma::scattering::SourceWgma> source,
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
    
    TPZPlanarWGScattering::ModeType matmode;
    switch(mode){
    case wgma::planarwg::mode::TE:
      matmode = TPZPlanarWGScattering::ModeType::TE;
      break;
    case wgma::planarwg::mode::TM:
      matmode = TPZPlanarWGScattering::ModeType::TM;
      break;
    }
    for(auto [id,er,ur] : data.matinfovec){
      auto *mat = new TPZPlanarWGScattering(id,er,ur,lambda,matmode,scale);
      scatt_cmesh->InsertMaterialObject(mat);
      //for pml
      realvolmats.insert(id);
      volmats.insert(id);
      allmats.insert(id);
    }
  
    for(auto pml : data.pmlvec){
      const auto id = pml.id;
      const auto alphax = pml.alphax;
      const auto alphay = pml.alphay;
      const auto type = pml.t;
      wgma::cmeshtools::AddRectangularPMLRegion<TPZPlanarWGScattering>
        (id, alphax, alphay, type, realvolmats, gmesh, scatt_cmesh);
      volmats.insert(id);
      allmats.insert(id);
    }

  
    TPZFNMatrix<1, CSTATE> val1(1, 1, 1);
    TPZManVector<CSTATE,1> val2(1, 0.);

    /**let us associate each boundary with a given material.
       this is important for the source boundary*/
    for(auto &bc : data.bcvec){
      for(auto *gel : gmesh->ElementVec()){
        if(gel->MaterialId() == bc.id){
          const auto maxside = gel->NSides() - 1;
          bc.volid = gel->Neighbour(maxside).Element()->MaterialId();
          break;
        }
      }
    }

    // for(auto bc : data.bcvec){
    //   std::cout<<"bc "<<bc.id<<" mat "<<bc.volid<<std::endl;
    // }
  
    for(auto bc : data.bcvec){
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

    //check whether the solution is analytical or if it comes from modal analysis
    const bool prescribed_source =
      std::holds_alternative<wgma::scattering::Source1D>(source);

    const std::set<int> src_id_set = [prescribed_source, &source](){
      if(prescribed_source){
        return std::get<wgma::scattering::Source1D>(source).id;
      }else{
        return std::get<wgma::scattering::SourceWgma>(source).id;
      }
    }();

    //we insert all the materials in the computational mesh
    for(auto src_id : src_id_set){
      auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, src_id, volmats);
      if(!res.has_value()){
        PZError<<__PRETTY_FUNCTION__
               <<"\n could not find a material adjacent to the source.\nAborting..."
               <<std::endl;
        DebugStop();
      }
  
      const auto src_volid = res.value();
      for(auto [id,er,ur] : data.matinfovec){
        if(id == src_volid){
          auto srcMat = new TPZPlanarWGScatteringSrc(src_id,er,ur,lambda,matmode,scale);
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


    //now we load the source
    LoadSources(scatt_cmesh, source, prescribed_source);
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
        dynamic_cast<TPZPlanarWGScatteringSrc*>(mat);
      if(scatt_mat){
        scatt_mat->SetBeta(beta);
      }
    }
  
  }
};