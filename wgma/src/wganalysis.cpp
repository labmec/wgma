#include "wganalysis.hpp"
#include "cmeshtools.hpp"

#include <TPZSpStructMatrix.h>
#include <TPZKrylovEigenSolver.h>
#include <Electromagnetics/TPZWaveguideModalAnalysis.h>
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>

namespace wgma{
  WGAnalysis::WGAnalysis(const TPZVec<TPZAutoPointer<TPZCompMesh>> &meshvec,
                         const int n_threads, const bool reorder_eqs,
                         const bool filter_bound) : m_filter_bound(filter_bound){
    if(meshvec.size() != 3){
      std::cerr<<__PRETTY_FUNCTION__
               <<"\nThree computational meshes are required."
               <<"Aborting...\n";
      exit(-1);
    }
    //gets the multiphysics mesh (main mesh)
    m_cmesh_mf = meshvec[0];
    m_cmesh_h1 = meshvec[TPZWaveguideModalAnalysis::H1Index()];
    m_cmesh_hcurl = meshvec[TPZWaveguideModalAnalysis::HCurlIndex()];
  
    m_an = new TPZEigenAnalysis(m_cmesh_mf, reorder_eqs);

    TPZAutoPointer<TPZStructMatrix> strmtrx =
      new TPZSpStructMatrix<CSTATE>(m_cmesh_mf);

    strmtrx->SetNumThreads(n_threads);
    
    TPZVec<int64_t> activeEquations;
  
    
    
    if(filter_bound){
      int n_dofs_before{-1};
      wgma::cmeshtools::FilterBoundaryEquations(meshvec, activeEquations,
                                                m_n_dofs_mf, n_dofs_before,
                                                m_n_dofs_h1, m_n_dofs_hcurl);
      std::cout<<"neq(before): "<<n_dofs_before
               <<"\tneq(after): "<<m_n_dofs_mf<<std::endl;
      strmtrx->EquationFilter().SetActiveEquations(activeEquations);
    }else{
      std::set<int64_t> boundConnects;
      wgma::cmeshtools::CountActiveEquations(meshvec,boundConnects,m_n_dofs_mf,
                                             m_n_dofs_h1,m_n_dofs_hcurl);
    }
    m_an->SetStructuralMatrix(strmtrx);
  }
    
  void WGAnalysis::SetSolver(TPZAutoPointer<TPZEigenSolver<CSTATE>> solv){
    m_solver = solv;
    m_an->SetSolver(*m_solver);
  }
    
  TPZAutoPointer<TPZEigenSolver<CSTATE>> WGAnalysis::GetSolver() const{
    return m_solver;
  }

  void WGAnalysis::Run(bool compute_eigenvectors){
    m_an->SetComputeEigenvectors(compute_eigenvectors);

    if(!m_solver){
      std::cerr<<__PRETTY_FUNCTION__
               <<"\nA solver has not been set.\n"
               <<"Check documentation of TPZKrylovEigenSolver"
               <<"\nor wgma::slepc::EPSHandler\nAborting..."<<std::endl;
      exit(-1);
    }

    auto krylov_solver =
      TPZAutoPointerDynamicCast<TPZKrylovEigenSolver<CSTATE>>(m_solver);
    if(krylov_solver && krylov_solver->KrylovInitialVector().Rows() == 0){
      /**this is to ensure that the eigenvector subspace is orthogonal to
         the spurious solutions associated with et = 0 ez != 0*/
      TPZFMatrix<CSTATE> initVec(m_n_dofs_mf, 1, 0.);
      const auto firstHCurl = m_n_dofs_h1 * TPZWaveguideModalAnalysis::HCurlIndex();
      for (int i = 0; i < m_n_dofs_hcurl; i++) {
        initVec(firstHCurl + i, 0) = 1;
      }
      krylov_solver->SetKrylovInitialVector(initVec);
    }
    
    {//scope for timer
      std::cout<<"Assembling..."<<std::flush;
      TPZSimpleTimer assemble("Assemble");
      m_an->Assemble();
      std::cout<<"\rAssembled!"<<std::endl;
    }
    {//scope for timer
      TPZSimpleTimer solv("Solve");
      std::cout<<"Solving..."<<std::flush;
      m_an->Solve();
      std::cout<<"\rSolved!"<<std::endl;
    }

    m_evalues = m_an->GetEigenvalues();

    std::ios cout_state(nullptr);
    cout_state.copyfmt(std::cout);
    
    std::cout << std::setprecision(std::numeric_limits<STATE>::max_digits10);

    for(auto &w : m_evalues){
      std::cout<<w<<std::endl;
    }
    std::cout.copyfmt(cout_state);
    
    if(m_an->ComputeEigenvectors()){
      m_evectors = m_an->GetEigenvectors();
    }
  }

  void WGAnalysis::PostProcess(std::string filename,
                               const int vtk_res,
                               const bool print_real_part){
    if(!m_an->ComputeEigenvectors()){
      std::cout<<__PRETTY_FUNCTION__
               <<"\nOnly eigenvalues were calculated.\n"
               <<"Nothing to do here...\n";
      return;
    }

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Ez");
    vecnames.Push("Et");
    const std::string plotfile = filename+".vtk";
    constexpr int dim{2};
    m_an->DefineGraphMesh(dim, scalnames, vecnames,plotfile);

    const auto neqOriginal = m_evectors.Rows();
    TPZFMatrix<CSTATE> evector(neqOriginal, 1, 0.);
    
    auto &ev = m_evalues;
    auto &eigenvectors = m_evectors;
    
    const auto nev = ev.size();

    TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
    meshVecPost[0] = m_cmesh_h1;
    meshVecPost[1] = m_cmesh_hcurl;
  
    std::cout<<"Post processing..."<<std::endl;
  
    for (int iSol = 0; iSol < ev.size(); iSol++) {
      const CSTATE currentKz = [&ev,iSol](){
        auto tmp = std::sqrt(-1.0*ev[iSol]);
        constexpr auto epsilon = std::numeric_limits<STATE>::epsilon()/
          (10*std::numeric_limits<STATE>::digits10);
        //let us discard extremely small imag parts
        if (tmp.imag() < epsilon)
          {tmp = tmp.real();}
        return tmp;
      }();
      eigenvectors.GetSub(0, iSol, neqOriginal, 1, evector);
      for(auto mat : m_cmesh_mf->MaterialVec()){
        auto id = mat.first;
        auto matPtr =
          dynamic_cast<TPZWaveguideModalAnalysis *>(m_cmesh_mf->FindMaterial(id));
        if(!matPtr) continue;
        matPtr->SetKz(currentKz);
        matPtr->SetPrintFieldRealPart(print_real_part);
      }
      std::cout<<"\rPost processing step "<<iSol+1<<" out of "<<ev.size()
               <<"(kz = "<<currentKz<<")"<<std::flush;
      m_an->LoadSolution(evector);
      TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshVecPost, m_cmesh_mf);
      m_an->PostProcess(vtk_res);
    }
    std::cout<<"\nFinished post processing"<<std::endl;
    std::cout<<std::endl;
  }
};