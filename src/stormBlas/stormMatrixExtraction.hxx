/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
#ifndef _STORM_MATRIX_EXTRACTION_HXX_
#define _STORM_MATRIX_EXTRACTION_HXX_

#include <algorithm>

#include <stormBlas/stormOperator.hxx>
#include <stormBlas/stormSubspace.hxx>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract matrix extractor.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class Matrix, class Vector>
class stormMatrixExtractor : public stormBaseObject {
public:

  /// @brief Extract a matrix from the operator.
  /// 
  /// @param mat Extracted matrix.
  /// @param anyVec Any vector, suitable as an input for the operator.
  /// @param linOp Linear operator.
  virtual void Extract(Matrix& mat,
                       Vector const& anyVec,
                       stormOperator<Vector> const& linOp);

}; // class stormMatrixExtractor

template<class Value, class Index = stormSize_t>
class stormSparseVector final {
public:
  std::vector<Index> RowIndices_;
  std::vector<Value> RowValues_;

public:
  void EmplaceEntry(Index const& index, 
                    Value const& value) {
    RowIndices_.emplace_back(index);
    RowValues_.emplace_back(value);
  }

}; // class stormSparseVector<...>

template<class Value, class Index = stormSize_t>
class stormSparseRowMatrix final {
public:
  std::vector<Index> RowPtrs_;
  std::vector<Index> ColIndices_;
  std::vector<Value> ColValues_;

public:
  stormSparseRowMatrix() {
    RowPtrs_.push_back(0);
  }

  void Clear() {
    RowPtrs_.clear();
    RowPtrs_.push_back(0);
    ColIndices_.clear(), ColValues_.clear();
  }

  void EmplaceRow(stormSparseVector<Value, Index> const& row) {
    ColIndices_.insert(ColIndices_.end(), row.RowIndices_.begin(), row.RowIndices_.end());
    ColValues_.insert(ColValues_.end(), row.RowValues_.begin(), row.RowValues_.end());
    RowPtrs_.push_back(ColIndices_.size());
  }

  template<class Vector>
  void mat_vec(Vector& yVec, Vector const& xVec) const {
    
    stormSize_t const n = xVec.Size();
    for (stormSize_t i = 0; i < n; ++i) {
      yVec(i) = 0.0;
      for (stormSize_t k = RowPtrs_[i]; k < RowPtrs_[i + 1]; ++k) {
        yVec(i) += ColValues_[k]*xVec(ColIndices_[k]);
      }
    }

  }

  void save() const {

  }

}; // class stormSparseRowMatrix<...>

template<class Vector>
void do_the_thing(stormSparseRowMatrix<stormReal_t>& mMat,
                  Vector const& anyVec,
                  stormOperator<Vector> const& linOp) {

  stormSize_t m = 401, n;

  n = anyVec.Size();

  Vector pVec_;
  stormSubspace<Vector> cVecs_;

  pVec_.Assign(anyVec, false);
  cVecs_.Assign(m, anyVec, false);

  // ----------------------
  // Prepare the probe vector and extract the 
  //   next set of the matrix columns:
  // 𝗳𝗼𝗿 𝑘 = 𝟢, 𝑚 - 𝟣 𝗱𝗼:
  //   𝒑 ← {𝟢}ᵀ,
  //   𝗳𝗼𝗿 𝑖 = 𝑘, 𝑛 - 𝟣, 𝑚 𝗱𝗼:
  //     𝒑ᵢ ← 𝟣,
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  //   𝒄ₖ ← 𝓐𝒑.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  for (stormSize_t k = 0; k < m; ++k) {
    stormBlas::Fill(pVec_, 0.0);
    for (stormSize_t i = k; i < n; i += m) {
      pVec_(i) = 1.0;
    }
    linOp.MatVec(cVecs_(k), pVec_);
  }

  // ----------------------
  // Assemble the matrix:
  // 𝓜 ← {}, // Sparse row matrix.
  // 𝗳𝗼𝗿 𝑖 = 𝟢, 𝑛 - 𝟣 𝗱𝗼:
  //   𝒓 ← {}, // Sparse row vector. 
  //   𝗳𝗼𝗿 𝑗 = 𝘮𝘢𝘹{ 𝟢, 𝑖 - ½⋅𝑚 }, 𝘮𝘪𝘯{ 𝑖 + ½⋅𝑚, 𝑛 - 𝟣 } 𝗱𝗼:
  //     𝑘 ← 𝑗 (𝘮𝘰𝘥 𝑚),
  //     𝗶𝗳 |𝒄ₖᵢ| ≥ 𝘛𝘩𝘳𝘦𝘴𝘩𝘰𝘭𝘥:
  //       𝒓ⱼ ← 𝒄ₖᵢ.
  //     𝗲𝗻𝗱 𝗶𝗳
  //   𝗲𝗻𝗱 𝗳𝗼𝗿
  //   𝓜ᵢ,: ← 𝒓.
  // 𝗲𝗻𝗱 𝗳𝗼𝗿
  // ----------------------
  mMat.Clear();
  for (stormSize_t i = 0; i < n; ++i) {
    stormSparseVector<stormReal_t> rVec;
    for (stormSize_t j = std::max<stormPtrDiff_t>(0, i - m/2); 
                     j <= std::min(i + m/2, n - 1); ++j) {
      stormSize_t const k = j % m;
      if (std::abs(cVecs_(k)(i)) > 1.0e-10/*threshold*/) {
        rVec.EmplaceEntry(j, cVecs_(k)(i));
      }
    }
    mMat.EmplaceRow(rVec);
  }
}

#endif // ifndef _STORM_MATRIX_EXTRACTION_HXX_
