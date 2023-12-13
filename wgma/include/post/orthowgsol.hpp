#ifndef _ORTHOWGSOL_H_
#define _ORTHOWGSOL_H_
#include <pzreal.h>
namespace wgma{
  namespace wganalysis{
    class Wgma;
  };

  namespace post{
    /*
      @brief Orthogonalises waveguide modes for expansion methods
      Modes with the same corresponding eigenvalues can be orthogonalised
      and are still modes.
      In the context of expansion methods, one often wants modes orthogonalised
      as x_i B x_j = 0 for i != j
      Which is precisely the orthogonalisation performed by this function
      @param [input] an analysis object
      @param [input] tol tolerance for
      @return number of orthogonalised modes
      @note Only valid for modal analysis whose resulting
      algebraic problem is a generalised eigenvalue problem
     */
    int OrthoWgSol(wganalysis::Wgma &an, const STATE tol);
  };
};

#endif /* _ORTHOWGSOL_H_ */
