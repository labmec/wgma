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
    // if (count==0) {
    //   continue;
    // }
    n_ortho+=count;

    const int offset = iev * neq_expand;
    TPZFMatrix<CSTATE> evec_scatter(neq_expand, count+1, evectors.Elem() + offset,
                                    neq_expand*(count+1));
    TPZFMatrix<CSTATE> evec(neq, count+1, 0);
    eqfilt.Gather(evec_scatter,evec);
    for(int i = 0; i < count+1;i++){
      const int offset = neq*i;
      TPZFMatrix<CSTATE> ei(neq,1,evec.Elem() + offset, neq);
      BMat->Multiply(ei, aux);
      for(int j = 0; j < i; j++){
        const int offset = neq*j;
        TPZFMatrix<CSTATE> ej(neq,1,evec.Elem() + offset, neq);
        CSTATE dot_ej{0};
        {
          const CSTATE *pej = ej.Elem();
          const CSTATE *pbei = aux.Elem();
	  if(conj){
	    for (int ieq = 0; ieq < neq; ieq++) {
	      dot_ej += std::conj(*pej++) * *pbei++;
	    }
	  }else{
	    for (int ieq = 0; ieq < neq; ieq++) {
	      dot_ej += *pej++ * *pbei++;
	    }
	  }
        }
        //v_j^T B v_j is unity
        const CSTATE &coeff = dot_ej;
        CSTATE *pei = ei.Elem();
        const CSTATE *pej = ej.Elem();
        for (int ieq = 0; ieq < neq; ieq++,pei++) {
          *pei -= coeff * *pej++;
        }
      }
      BMat->Multiply(ei, aux);
      CSTATE dot_ei{0};
      //now we normalise it
      {
        const CSTATE *pei = ei.Elem();
        const CSTATE *pbei = aux.Elem();
	if(conj){
	  for (int ieq = 0; ieq < neq; ieq++) {
	    dot_ei += std::conj(*pei++) * *pbei++;
	  }
	}else{
	  for (int ieq = 0; ieq < neq; ieq++) {
	    dot_ei += *pei++ * *pbei++;
	  }
	}
      }
      dot_ei = sqrt(dot_ei);
      {
        CSTATE *pei = ei.Elem();
        for (int ieq = 0; ieq < neq; ieq++, pei++) {
          *pei /= dot_ei;
        }
      }
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
