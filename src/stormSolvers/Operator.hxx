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
#include <stormSolvers/Matrix.hxx>
#include <stormSolvers/Vector.hxx>
#include <stormUtils/Object.hxx>

namespace Storm {

class stormArray : public BaseMatrix<stormArray> {
public:

  stormMesh_t Mesh = nullptr;
  stormArray_t Array = nullptr;
  std::shared_ptr<int> RefCounter;
  real_t* MyData = nullptr;
  size_t MySize = 0;

public:

  stormArray() = default;

  stormArray(stormMesh_t mesh, stormArray_t array) : Mesh(mesh), Array(array) {
    RefCounter = std::make_shared<int>(2);
    stormArrayUnwrap(Mesh, Array, &MyData, &MySize);
  }

  stormArray(stormArray&& oth) : Mesh(oth.Mesh), Array(oth.Array) {
    RefCounter = std::move(oth.RefCounter);
    oth.Mesh = nullptr, oth.Array = nullptr, oth.RefCounter = nullptr;
    stormArrayUnwrap(Mesh, Array, &MyData, &MySize);
  }

  stormArray(stormArray const& oth) : Mesh(oth.Mesh), Array(oth.Array) {
    RefCounter = oth.RefCounter;
    *RefCounter += 1;
    stormArrayUnwrap(Mesh, Array, &MyData, &MySize);
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

  auto shape() const noexcept {
    return std::pair(MySize, std::integral_constant<size_t, 1>{});
  }

  auto operator()(size_t i, size_t j) const noexcept {
    StormAssert(j == 0);
    return MyData[i];
  }
  auto& operator()(size_t i, size_t j) noexcept {
    StormAssert(j == 0);
    return MyData[i];
  }

  void assign(stormArray const& like, bool copy = true) {
    Mesh = like.Mesh;
    Array = stormAllocLike(like.Array);
    stormArrayUnwrap(Mesh, Array, &MyData, &MySize);
    if (copy) stormSet(Mesh, Array, like.Array);
  }
};

template<>
class VectorOperations<stormArray> {
public:

  static real_t Dot(stormArray const& z, stormArray const& y) {
    return stormDot(z.Mesh, z.Array, y.Array);
  }
  static real_t Norm2(stormArray const& z) {
    return stormNorm2(z.Mesh, z.Array);
  }

  static void Fill(stormArray& z, real_t a) {
    stormFill(z.Mesh, z.Array, a);
  }
  static void RandFill(stormArray& z) {
    stormRandFill(z.Mesh, z.Array);
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
  /// @param y_vec Output vector, 𝒚.
  /// @param x_vec Input vector, 𝒙.
  virtual void MatVec(OutVector& y_vec, InVector const& x_vec) const = 0;

  /// @brief Compute a chained
  ///   operator-vector product, 𝒛 ← 𝓐(𝒚 ← 𝓑(𝒙)).
  ///
  /// @param z_vec Output vector, 𝒛.
  /// @param y_vec Intermediate vector, 𝒚.
  /// @param x_vec Input vector, 𝒙.
  template<VectorLike InOutVector = InVector>
  void MatVec(OutVector& z_vec, InOutVector& y_vec,
              Operator<InVector, InOutVector> const& otherOp,
              InVector const& x_vec) const {
    otherOp.MatVec(y_vec, x_vec);
    MatVec(z_vec, y_vec);
  }

  /// @brief Compute a residual, 𝒓 ← 𝒃 - 𝓐(𝒙).
  ///
  /// @param r_vec Residual vector, 𝒓.
  /// @param b_vec Input vector, 𝒃.
  /// @param x_vec Input vector, 𝒙.
  void Residual(OutVector& r_vec, OutVector const& b_vec,
                InVector const& x_vec) const {
    MatVec(r_vec, x_vec);
    r_vec <<= b_vec - r_vec;
  }

  /// @brief Compute a residual norm, ‖𝒃 - 𝓐𝒙‖.
  ///
  /// @param b_vec Input vector, 𝒃.
  /// @param x_vec Input vector, 𝒙.
  real_t ResidualNorm(OutVector const& b_vec, InVector const& x_vec) const {
    OutVector r_vec;
    r_vec.assign(b_vec, false);
    Residual(r_vec, b_vec, x_vec);
    return r_vec.Norm2();
  }

  /// @brief Compute an conjugate operator-vector product, 𝒙 ← 𝓐*(𝒚).
  ///
  /// @param x_vec Output vector, 𝒙.
  /// @param y_vec Input vector, 𝒚.
  virtual void ConjMatVec(InVector& x_vec, OutVector const& y_vec) const {
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

  void MatVec(OutVector& y_vec, InVector const& x_vec) const override {
    MatVecFunc_(y_vec, x_vec);
  }

  void ConjMatVec(InVector& x_vec, OutVector const& y_vec) const override {
    if (!ConjMatVecFunc_) {
      throw std::runtime_error("`FunctionalOperator::ConjMatVec`"
                               " conjugate product function was not set.");
    }
    ConjMatVecFunc_(x_vec, y_vec);
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
