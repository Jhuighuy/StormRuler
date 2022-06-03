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
};

template<>
class VectorOperations<stormArray> {
public:

  static void Swap(auto& x, auto& y) {
    std::swap(x, y);
  }

  static real_t Dot(stormArray const& z, stormArray const& y) {
    return stormDot(z.Mesh, z.Array, y.Array);
  }
  static real_t Norm2(stormArray const& z) {
    return stormNorm2(z.Mesh, z.Array);
  }

  static void Set(stormArray& z, stormArray const& y) {
    stormSet(z.Mesh, z.Array, y.Array);
  }

  static void Fill(stormArray& z, real_t a) {
    stormFill(z.Mesh, z.Array, a);
  }
  static void RandFill(stormArray& z) {
    stormRandFill(z.Mesh, z.Array);
  }

  static void Scale(stormArray& z, stormArray const& y, real_t a) {
    stormScale(z.Mesh, z.Array, y.Array, a);
  }

  static void ScaleAssign(stormArray& z, real_t a) {
    Scale(z, z, a);
  }

  static void Add(stormArray& z, stormArray const& y, stormArray const& x) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array);
  }
  static void Add(stormArray& z, stormArray const& y, stormArray const& x,
                  real_t a) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a);
  }
  static void Add(stormArray& z, stormArray const& y, real_t b,
                  stormArray const& x, real_t a) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }

  static void AddAssign(stormArray& z, stormArray const& y) {
    Add(z, z, y);
  }
  static void AddAssign(stormArray& z, stormArray const& y, real_t a) {
    Add(z, z, y, a);
  }

  static void Sub(stormArray& z, stormArray const& y, stormArray const& x) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array);
  }
  static void Sub(stormArray& z, stormArray const& y, stormArray const& x,
                  real_t a) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a);
  }
  static void Sub(stormArray& z, stormArray const& y, real_t b,
                  stormArray const& x, real_t a) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }

  static void SubAssign(stormArray& z, stormArray const& y) {
    Sub(z, z, y);
  }
  static void SubAssign(stormArray& z, stormArray const& y, real_t a) {
    Sub(z, z, y, a);
  }
}; // class BlasOp<stormArray>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator 𝒚 ← 𝓐(𝒙).
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<VectorLike InVector, VectorLike OutVector = InVector>
class Operator : public Object {
public:

  /// @brief Compute an operator-vector product, 𝒚 ← 𝓐(𝒙).
  ///
  /// @param yVec Output vector, 𝒚.
  /// @param xVec Input vector, 𝒙.
  virtual void MatVec(OutVector& yVec, InVector const& xVec) const = 0;

  /// @brief Compute a chained
  ///   operator-vector product, 𝒛 ← 𝓐(𝒚 ← 𝓑(𝒙)).
  ///
  /// @param zVec Output vector, 𝒛.
  /// @param yVec Intermediate vector, 𝒚.
  /// @param xVec Input vector, 𝒙.
  template<VectorLike InOutVector = InVector>
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
    Blas::Sub(rVec, bVec, rVec);
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
template<VectorLike InVector, VectorLike OutVector = InVector>
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
template<VectorLike InVector, VectorLike OutVector = InVector,
         operator_like<InVector, OutVector> MatVecFunc>
auto MakeOperator(MatVecFunc&& matVecFunc) {
  return std::make_unique<FunctionalOperator<InVector, OutVector>>(
      std::forward<MatVecFunc>(matVecFunc));
}
template<VectorLike InVector, VectorLike OutVector = InVector,
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
template<VectorLike Vector, operator_like<Vector> MatVecFunc>
auto MakeSymmetricOperator(MatVecFunc&& matVecFunc) {
  return std::make_unique<FunctionalOperator<Vector>>(
      matVecFunc, std::forward<MatVecFunc>(matVecFunc));
}

} // namespace Storm
