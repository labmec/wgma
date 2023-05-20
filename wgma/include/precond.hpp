#ifndef _PRECOND_HPP_
#define _PRECOND_HPP_

#include <pzreal.h>
#include <tpzautopointer.h>
#include <TPZMatrixSolver.h>

#include <set>
#include <cstdint>

class TPZEquationFilter;

//! Contains a set of routines for creating and manipulating preconditioners
namespace wgma::precond{
  /**
     @brief This function will create AFW blocks for preconditioning a matrix created
     with hcurl elements
  */
  void CreateAFWBlocks(TPZAutoPointer<TPZCompMesh> cmesh,
                       const std::set<int> dirichlet_mats,
                       const TPZEquationFilter &eqfilt,
                       TPZVec<int64_t> &eqgraph,
                       TPZVec<int64_t> &eqgraphindex
                       );

  int ColorEqGraph(const TPZVec<int64_t> &graph, const TPZVec<int64_t> &graphindex,
                   const TPZMatrix<CSTATE> &mat,
                   const int64_t neqs, TPZVec<int> &colors);

  class BlockPrecond : public TPZMatrixSolver<CSTATE>{
  public:
    BlockPrecond(TPZAutoPointer<TPZMatrix<CSTATE>> refmat,
                 TPZVec<int64_t> &&block,
                 TPZVec<int64_t> &&block_index,
                 const TPZVec<int> &colors,
                 const int numcolors);
    
    BlockPrecond(TPZAutoPointer<TPZMatrix<CSTATE>> refmat,
                 const TPZVec<int64_t> &block,
                 const TPZVec<int64_t> &block_index,
                 const TPZVec<int> &colors,
                 const int numcolors);

    BlockPrecond(const BlockPrecond&);
    BlockPrecond(BlockPrecond&&) = default;
    BlockPrecond& operator=(const BlockPrecond&);
    BlockPrecond& operator=(BlockPrecond&&) = default;
    ~BlockPrecond()=default;

    void UpdateFrom(TPZAutoPointer<TPZBaseMatrix> mat_base) override;

    void Solve(const TPZFMatrix<CSTATE> &F,
               TPZFMatrix<CSTATE> &result,
               TPZFMatrix<CSTATE> *residual = 0) override;

    //! sets symmetric gauss seidel
    inline void SetSymmetricGS(){m_sym_gs = true;}
    //! sets traditional gauss seidel
    inline void SetTraditionalGS(){m_sym_gs = false;}

    inline int NBlocks() const {return m_blockindex.size()-1;}
    inline int BlockSize(const int i) const {return m_blockindex[i+1]-m_blockindex[i];}

    BlockPrecond *Clone() const override{
      return new BlockPrecond(*this);
    }
  protected:
    void SmoothBlock(const int bl, TPZFMatrix<CSTATE> &du,
                     const TPZFMatrix<CSTATE>&rhs);
    //! Computes the list of blocks by color
    void ComputeColorVec(const TPZVec<int> &colors, const int nc);
    //! whether to perform symmetric gauss seidel or not
    bool m_sym_gs{true};
    //! max block size
    size_t m_max_bs{0};
    //! position i will have the first position in blockgraph of block i
    TPZVec<int64_t> m_blockindex;
    //! all equations for the blocks
    TPZVec<int64_t> m_blockgraph;
    //! position i will have all the blocks with color i
    TPZVec<TPZVec<int64_t>> m_colors;
    //! stores the inverses of each block diag matrix
    TPZVec<TPZAutoPointer<TPZFMatrix<CSTATE>>> m_blockinv;
  };
};

#endif /* _PRECOND_HPP_ */
