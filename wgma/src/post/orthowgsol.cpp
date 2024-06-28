#include <post/orthowgsol.hpp>
#include <TPZLinearEigenSolver.h>
#include <wganalysis.hpp>

int wgma::post::OrthoWgSol(wgma::wganalysis::Wgma &an,
                           const STATE tol,
                           const bool conj){
  int n_ortho{0};
  auto BMat = [&an]() -> TPZAutoPointer<TPZMatrix<CSTATE>> {
    auto an_cast = dynamic_cast<wgma::wganalysis::Wgma2D*>(&an);
    if(an_cast){
      return an_cast->GetSolver().MatrixB();
    }else{
      auto an_cast = dynamic_cast<wgma::wganalysis::WgmaPlanar*>(&an);
      if(an_cast){
        return an_cast->GetSolver().MatrixB();
      }else{
        DebugStop();
        return nullptr;
      }
    }
  }();
  const auto neq = BMat->Rows();
  // we get equation filter for gather/scatter
  auto &eqfilt = an.StructMatrix()->EquationFilter();
  // we only orthogonalise eigenvectors if eigenvalue are equal
  const auto &evalues = an.GetEigenvalues();
  auto &evectors = an.GetEigenvectors();
  const auto nev = evalues.size();
  const auto neq_expand = evectors.Rows();
  TPZFMatrix<CSTATE> aux(neq, 1, 0);
  for (int iev = 0; iev < nev; iev++) {
    int count{0};
    for (int iiev = iev + 1; iiev < nev; iiev++) {
      const auto diff = evalues[iev] - evalues[iiev];
      const auto diff_sq = (diff * std::conj(diff)).real();
      if (diff_sq < tol) {
        count++;
      } else {
        break;
      }
    }
    if(count!=0){
      n_ortho+=count+1;
    }

    const int offset_i_scatter = iev * neq_expand;
    TPZFMatrix<CSTATE> evec_scatter(neq_expand, count+1, evectors.Elem() + offset_i_scatter,
                                    neq_expand*(count+1));
    TPZFMatrix<CSTATE> evec(neq, count+1, 0);
    eqfilt.Gather(evec_scatter,evec);
    for(int i = 0; i < count+1;i++){
      const int offset_i = neq*i;
      TPZFMatrix<CSTATE> ei(neq,1,evec.Elem() + offset_i, neq);
      for(int j = 0; j < i; j++){
	const int offset_j = neq*j;
        TPZFMatrix<CSTATE> ej(neq,1,evec.Elem() + offset_j, neq);
        BMat->Multiply(ej, aux);
        const CSTATE dot_ej = Dot(ei,aux,conj);
        //v_j^T B v_j is unity, so we dont need to divide by it
        ei-=dot_ej*ej;
      }
      //now we normalise it
      BMat->Multiply(ei, aux);
      const CSTATE dot_ei = Dot(ei,aux,conj);
      const CSTATE coeff = 1./sqrt(dot_ei);
      ei*=coeff;
      TPZFMatrix<CSTATE> ei_scatter(neq_expand, 1,
                                    evec_scatter.Elem() + neq_expand*i,
                                    neq_expand);
      eqfilt.Scatter(ei,ei_scatter);
    }
    iev+=count;
  }
  an.LoadAllSolutions();
  return n_ortho;
}
