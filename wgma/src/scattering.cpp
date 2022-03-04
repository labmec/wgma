#include "scattering.hpp"
#include "cmeshtools.hpp"

#include <TPZSpStructMatrix.h>
#include <pzstepsolver.h>
#include <Electromagnetics/TPZPlanarWGScattering.h>
#include <TPZSimpleTimer.h>

namespace wgma{
  ScatteringAnalysis::ScatteringAnalysis(TPZAutoPointer<TPZCompMesh> mesh,
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
    
  void ScatteringAnalysis::SetSolver(const TPZMatrixSolver<CSTATE> & solv){
    m_an->SetSolver(solv);
  }
    
  TPZMatrixSolver<CSTATE> & ScatteringAnalysis::GetSolver() const{
    return m_an->MatrixSolver<CSTATE>();
  }

  void ScatteringAnalysis::Run(){
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

  void ScatteringAnalysis::PostProcess(std::string filename,
                               const int vtk_res){


    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Field_re");
    scalnames.Push("Field_abs");
    const std::string plotfile = filename+".vtk";
    constexpr int dim{2};
    m_an->DefineGraphMesh(dim, scalnames, vecnames,plotfile);
    m_an->PostProcess(vtk_res);
    std::cout<<"\nFinished post processing"<<std::endl;
    std::cout<<std::endl;
  }
};