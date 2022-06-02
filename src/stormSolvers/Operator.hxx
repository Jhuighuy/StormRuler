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

#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>

#include <stormBase.hxx>
#include <stormSolvers/Vector.hxx>
#include <stormUtils/Object.hxx>

namespace Storm {

class stormArray {
public:

  stormMesh_t Mesh = nullptr;
  stormArray_t Array = nullptr;
  std::shared_ptr<int> RefCounter;

public:

  stormArray() = default;

  stormArray(stormMesh_t mesh, stormArray_t array) : Mesh(mesh), Array(array) {
    RefCounter = std::make_shared<int>(2);
  }

  stormArray(stormArray&& oth) : Mesh(oth.Mesh), Array(oth.Array) {
    RefCounter = std::move(oth.RefCounter);
    oth.Mesh = nullptr, oth.Array = nullptr, oth.RefCounter = nullptr;
  }

  stormArray(stormArray const& oth) : Mesh(oth.Mesh), Array(oth.Array) {
    RefCounter = oth.RefCounter;
    *RefCounter += 1;
  }

  ~stormArray() {
    if (RefCounter) {
      *RefCounter -= 1;
      if (*RefCounter == 0) stormFree(Array);
    }
  }

  stormArray& operator=(stormArray&& oth) {
    this->~stormArray();
    new (this) stormArray(std::forward<stormArray>(oth));
    return *this;
  }
  stormArray& operator=(stormArray const& oth) {
    this->~stormArray();
    new (this) stormArray(oth);
    return *this;
  }

  void Assign(stormArray const& like, bool copy = true) {
    Mesh = like.Mesh;
    Array = stormAllocLike(like.Array);
    if (copy) stormSet(Mesh, Array, like.Array);
  }

  void Swap(stormArray& y) {
    std::swap(*this, y);
  }

  real_t Dot(stormArray const& y) const {
    return stormDot(Mesh, Array, y.Array);
  }
  real_t Norm2() const {
    return stormNorm2(Mesh, Array);
  }

  void Set(stormArray const& y) {
    stormSet(Mesh, Array, y.Array);
  }

  void Fill(real_t a) {
    stormFill(Mesh, Array, a);
  }
  void RandFill() {
    stormRandFill(Mesh, Array);
  }

  void Scale(stormArray const& y, real_t a) {
    stormScale(Mesh, Array, y.Array, a);
  }

  void ScaleAssign(real_t a) {
    Scale(*this, a);
  }

  void Add(stormArray const& y, stormArray const& x) {
    stormAdd(Mesh, Array, y.Array, x.Array);
  }
  void Add(stormArray const& y, stormArray const& x, real_t a) {
    stormAdd(Mesh, Array, y.Array, x.Array, a);
  }
  void Add(stormArray const& y, real_t b, stormArray const& x, real_t a) {
    stormAdd(Mesh, Array, y.Array, x.Array, a, b);
  }

  void AddAssign(stormArray const& y) {
    Add(*this, y);
  }
  void AddAssign(stormArray const& y, real_t a) {
    Add(*this, y, a);
  }

  void Sub(stormArray const& y, stormArray const& x) {
    stormSub(Mesh, Array, y.Array, x.Array);
  }
  void Sub(stormArray const& y, stormArray const& x, real_t a) {
    stormSub(Mesh, Array, y.Array, x.Array, a);
  }
  void Sub(stormArray const& y, real_t b, stormArray const& x, real_t a) {
    stormSub(Mesh, Array, y.Array, x.Array, a, b);
  }

  void SubAssign(stormArray const& y) {
    Sub(*this, y);
  }
  void SubAssign(stormArray const& y, real_t a) {
    Sub(*this, y, a);
  }
};

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator 𝒚 ← 𝓐(𝒙).
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<vector_like InVector, vector_like OutVector = InVector>
class Operator : public Object {
public:

  /// @brief Compute an operator-vector product, 𝒚 ← 𝓐(𝒙).
  ///
  /// @param yVec Output vector, 𝒚.
  /// @param xVec Input vector, 𝒙.
  virtual void MatVec(OutVector& yVec, InVector const& xVec) const = 0;

  /// @brief Compute a chained \
  ///   operator-vector product, 𝒛 ← 𝓐(𝒚 ← 𝓑(𝒙)).
  ///
  /// @param zVec Output vector, 𝒛.
  /// @param yVec Intermediate vector, 𝒚.
  /// @param xVec Input vector, 𝒙.
  template<vector_like InOutVector = InVector>
  void MatVec(OutVector& zVec, InOutVector& yVec,
              Operator<InVector, InOutVector> const& otherOp,
              InVector const& xVec) const {
    otherOp.MatVec(yVec, xVec);
    MatVec(zVec, yVec);
  }

  /// @brief Compute a residual, 𝒓 ← 𝒃 - 𝓐(𝒙).
  ///
  /// @param rVec Residual vector, 𝒓.
  /// @param bVec Input vector, 𝒃.
  /// @param xVec Input vector, 𝒙.
  void Residual(OutVector& rVec, OutVector const& bVec,
                InVector const& xVec) const {
    MatVec(rVec, xVec);
    rVec.Sub(bVec, rVec);
  }

  /// @brief Compute a residual norm, ‖𝒃 - 𝓐𝒙‖.
  ///
  /// @param bVec Input vector, 𝒃.
  /// @param xVec Input vector, 𝒙.
  real_t ResidualNorm(OutVector const& bVec, InVector const& xVec) const {
    OutVector rVec;
    rVec.Assign(bVec, false);
    Residual(rVec, bVec, xVec);
    return rVec.Norm2();
  }

  /// @brief Compute an conjugate operator-vector product, 𝒙 ← 𝓐*(𝒚).
  ///
  /// @param xVec Output vector, 𝒙.
  /// @param yVec Input vector, 𝒚.
  virtual void ConjMatVec(InVector& xVec, OutVector const& yVec) const {
    throw std::runtime_error("`Operator::ConjMatVec` was not overriden");
  }

}; // class Operator

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Operator implementation with external function pointers.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<vector_like InVector, vector_like OutVector = InVector>
class FunctionalOperator final : public Operator<InVector, OutVector> {
private:

  std::function<void(OutVector&, InVector const&)> MatVecFunc_;
  std::function<void(InVector&, OutVector const&)> ConjMatVecFunc_;

public:

  /// @brief Construct the functional operator.
  ///
  /// @param matVecFunc Operator-vector product function, 𝒚 ← 𝓐(𝒙).
  /// @param conjMatVecFunc Conjugate operator-vector product, 𝒙 ← 𝓐*(𝒚).
  /// @{
  template<operator_like<InVector, OutVector> MatVecFunc>
  explicit FunctionalOperator(MatVecFunc&& matVecFunc)
      : MatVecFunc_{std::forward<MatVecFunc>(matVecFunc)} {
    StormAssert(MatVecFunc_);
  }
  template<operator_like<InVector, OutVector> MatVecFunc,
           operator_like<OutVector, InVector> ConjMatVecFunc>
  explicit FunctionalOperator(MatVecFunc&& matVecFunc,
                              ConjMatVecFunc&& conjMatVecFunc)
      : MatVecFunc_{std::forward<MatVecFunc>(matVecFunc)},
        ConjMatVecFunc_{std::forward<ConjMatVecFunc>(conjMatVecFunc)} {
    StormAssert(MatVecFunc_ && ConjMatVecFunc_);
  }
  /// @}

private:

  void MatVec(OutVector& yVec, InVector const& xVec) const override {
    MatVecFunc_(yVec, xVec);
  }

  void ConjMatVec(InVector& xVec, OutVector const& yVec) const override {
    if (!ConjMatVecFunc_) {
      throw std::runtime_error("`FunctionalOperator::ConjMatVec`"
                               " conjugate product function was not set.");
    }
    ConjMatVecFunc_(xVec, yVec);
  }

}; // class FunctionalOperator

/// ----------------------------------------------------------------- ///
/// @brief Make the functional operator.
///
/// @param matVecFunc Operator-vector product function, 𝒚 ← 𝓐(𝒙).
/// @param conjMatVecFunc Conjugate operator-vector product, 𝒙 ← 𝓐*(𝒚).
/// ----------------------------------------------------------------- ///
/// @{
template<vector_like InVector, vector_like OutVector = InVector,
         operator_like<InVector, OutVector> MatVecFunc>
auto MakeOperator(MatVecFunc&& matVecFunc) {
  return std::make_unique<FunctionalOperator<InVector, OutVector>>(
      std::forward<MatVecFunc>(matVecFunc));
}
template<vector_like InVector, vector_like OutVector = InVector,
         operator_like<InVector, OutVector> MatVecFunc,
         operator_like<OutVector, InVector> ConjMatVecFunc>
auto MakeOperator(MatVecFunc&& matVecFunc, ConjMatVecFunc&& conjMatVecFunc) {
  return std::make_unique<FunctionalOperator<InVector, OutVector>>(
      std::forward<MatVecFunc>(matVecFunc),
      std::forward<ConjMatVecFunc>(conjMatVecFunc));
}
/// @}

/// ----------------------------------------------------------------- ///
/// @brief Make the self-adjoint functional operator.
/// ----------------------------------------------------------------- ///
template<vector_like Vector, operator_like<Vector> MatVecFunc>
auto MakeSymmetricOperator(MatVecFunc&& matVecFunc) {
  return std::make_unique<FunctionalOperator<Vector>>(
      matVecFunc, std::forward<MatVecFunc>(matVecFunc));
}

} // namespace Storm
