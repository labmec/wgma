#include "post/orthosol.hpp"
#include "post/solutionnorm.hpp"

#include <TPZMaterial.h>
#include <TPZMaterialDataT.h>
#include <TPZBndCond.h>
#include <pzvec_extras.h>

namespace wgma::post{

  template<class TSPACE>
  TPZFMatrix<CSTATE> OrthoSol<TSPACE>::Orthogonalise()
  {

    auto mesh = this->Mesh();

    TPZFMatrix<CSTATE> evectors = mesh->Solution();
    const int nev = evectors.Cols();
    const int neq = mesh->NEquations();
    //we iterate through the eigenvectors
    //if there is only one, it will only be normalised
    for(int iev = 0; iev < nev; iev++){
      //load the current solutions at the mesh
      TPZFMatrix<CSTATE> currentsol(neq, iev+1, evectors.Elem(), neq*(iev+1));
      mesh->LoadSolution(currentsol);
      //updates index of solution being normalised
      this->SetSol(iev);
      
      //gets eigenvector being orthogonalised
      const int offset = iev * neq;
      TPZFMatrix<CSTATE> ei(neq,1,evectors.Elem() + offset,neq);
      
      //how many rows the results matrix has
      const int nrows_res = std::max(this->NThreads(),1);
      //if it is the first eigenvector, there is no need to iterate through
      //all the elements, only normalising it
      if(iev>0){
        m_res.Redim(nrows_res, iev + 1);
        this->Integrate(this->m_elvec);
        //this will contain the dot product of the solutions
        TPZVec<CSTATE> res(iev+1,0.);
        //sum the results
        const int nrows = m_res.Rows();
        for(int i = 0; i < nrows; i++){
          for(int jev = 0; jev < iev; jev++){
            res[jev] += m_res(i,jev);
          }
        }
        //gram-schmidt iteration
        for(int jev = 0; jev < iev; jev++){
          const int offset = jev * neq;
          TPZFMatrix<CSTATE> ej(neq,1,evectors.Elem() + offset,neq);
          ei -= res[jev] * ej;
        }
      }

      auto solnorm = SolutionNorm<TSPACE>(mesh, this->m_elvec, this->NThreads());
      auto norm = solnorm.ComputeNorm(iev);
      ei *= 1./norm;
    }

    return std::move(evectors);
  }

  template<class TSPACE>
  void OrthoSol<TSPACE>::Compute(const ElData &eldata,REAL weight, int index)
  {
    const auto &which = WhichSol();
    const TPZMaterialDataT<CSTATE> &data = eldata;
    const auto &cursol = data.sol[which];
    const auto solsize = cursol.size();
    for(auto isol = 0; isol < which; isol++){
      for(auto ix = 0; ix < solsize; ix++){
        const auto val = cursol[ix] * std::conj(data.sol[isol][ix]);
        m_res(index,isol) += weight * fabs(data.detjac) * val;
      }
    }
  }

  template
  class OrthoSol<SingleSpaceIntegrator>;
  template
  class OrthoSol<MultiphysicsIntegrator>;
};
