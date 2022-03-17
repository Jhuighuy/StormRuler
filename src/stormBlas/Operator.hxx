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

#include <memory>
#include <stdexcept>
#include <functional>
#include <cmath>

#include <stormBase.hxx>

namespace Storm {

class BaseObject {
public:
  virtual ~BaseObject() = default;
}; // class BaseObject

class stormArray {
public:
  stormMesh_t Mesh = nullptr;
  stormArray_t Array = nullptr;
  std::shared_ptr<int> RefCounter;
private:
  real_t* Data_ = nullptr;
  size_t Size_ = 0;

public:
  stormArray() = default;
  stormArray(stormMesh_t mesh, stormArray_t array): Mesh(mesh), Array(array) {
    RefCounter = std::make_shared<int>(2);
    stormArrayUnwrap(Array, &Data_, &Size_);
  }
  stormArray(stormArray&& oth): Mesh(oth.Mesh), Array(oth.Array) {
    RefCounter = std::move(oth.RefCounter);
    oth.Mesh = nullptr, oth.Array = nullptr, oth.RefCounter = nullptr;
    Data_ = oth.Data_, Size_ = oth.Size_;
  }
  stormArray(stormArray const& oth): Mesh(oth.Mesh), Array(oth.Array) {
    RefCounter = oth.RefCounter;
    *RefCounter += 1;
    Data_ = oth.Data_, Size_ = oth.Size_;
  }
  ~stormArray() {
    if (RefCounter) {
      *RefCounter -= 1;
      if (*RefCounter == 0) stormFree(Array);
    }
  }

  size_t Size() const noexcept {
    return Size_;
  }
  real_t& operator()(size_t index) {
    //stormAssert(index < Size_);
    return Data_[index];
  }
  real_t const& operator()(size_t index) const {
    //stormAssert(index < Size_);
    return Data_[index];
  }

  stormArray& operator=(stormArray&& oth) {
    this->~stormArray();
    new(this) stormArray(std::forward<stormArray>(oth));
    return *this;
  }
  stormArray& operator=(stormArray const& oth) {
    this->~stormArray();
    new(this) stormArray(oth);
    return *this;
  }

  void Assign(stormArray const& like, bool copy = true) {
    Mesh = like.Mesh;
    Array = stormAllocLike(like.Array);
    stormArrayUnwrap(Array, &Data_, &Size_);
    if (copy) stormSet(Mesh, Array, like.Array);
  }

};

template<class Value>
class Array {
private:
  size_t Size_ = 0;
  std::vector<Value> Data_;

public:

  explicit Array(size_t size) :
      Size_(size), Data_(size) {
  }
  Array(size_t size, size_t fullSize) :
      Size_(size), Data_(fullSize) {
  }

  Array() = default;

  Array(Array const&) = default;

  Array(Array&&) noexcept = default;

  ~Array() = default;

  Array& operator=(Array&&) = default;

  Array& operator=(Array const&) = default;

  void Assign(Array const& array, bool copy = true) {
    Size_ = array.Size_;
    if (copy) {
      Data_ = array.Data_;
    } else {
      Data_.resize(array.Data_.size());
    }
  }

  size_t Size() const noexcept {
    return Size_;
  }

  size_t FullSize() const noexcept {
    return Data_.size();
  }

  Value& operator()(size_t index) noexcept {
    StormAssert(index < Data_.size());
    return Data_[index];
  }
  Value const& operator()(size_t index) const noexcept {
    StormAssert(index < Data_.size());
    return Data_[index];
  }

}; // class Array<...>

namespace Utils {
  
  real_t SafeDivide(real_t x, real_t y) {
    return (y == 0.0) ? 0.0 : (x/y);
  }

  real_t& SafeDivideEquals(real_t& x, real_t y) {
    x = SafeDivide(x, y);
    return x;
  }

} // namespace Utils

#if 0
inline void ToArray(Array<real_t>& to, stormArray const& from) {
  for (size_t i = 0; i < 15906; ++i) {
    to(i) = from(i); 
  }
}
inline void FromArray(stormArray& to, 
                      Array<real_t> const& from) {
  for (size_t i = 0; i < 15906; ++i) {
    to(i) = from(i); 
  }
}

template<class Func>
void ParallelFor(size_t first, size_t last, Func&& func) {
  for (size_t index = first; index <= last; ++index) {
    func(index);
  }
}

template<class Value, class Func>
auto ParallelSum(size_t first, size_t last, Value init, Func&& func) {
  for (size_t index = first; index <= last; ++index) {
    init += func(index);
  }
  return init;
}
#endif

namespace Blas {

  real_t Dot(stormArray const& z, stormArray const& y) {
    return stormDot(z.Mesh, z.Array, y.Array);
  }
  real_t Norm2(stormArray const& z) {
    return stormNorm2(z.Mesh, z.Array);
  }

  void Set(stormArray& z, stormArray const& y) {
    stormSet(z.Mesh, z.Array, y.Array);
  }

  void Fill(stormArray& z, real_t a) {
    stormFill(z.Mesh, z.Array, a);
  }
  void RandFill(stormArray& z) {
    stormRandFill(z.Mesh, z.Array);
  }

  void Scale(stormArray& z, stormArray const& y, real_t a) {
    stormScale(z.Mesh, z.Array, y.Array, a);
  }

  void Add(stormArray& z, 
           stormArray const& y, 
           stormArray const& x) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array);
  }
  void Add(stormArray& z, 
           stormArray const& y, 
           stormArray const& x, real_t a) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a);
  }
  void Add(stormArray& z, 
           stormArray const& y, real_t b,
           stormArray const& x, real_t a) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }

  void Sub(stormArray& z, 
           stormArray const& y, 
           stormArray const& x) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array);
  }
  void Sub(stormArray& z, 
           stormArray const& y, 
           stormArray const& x, real_t a) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a);
  }
  void Sub(stormArray& z, 
           stormArray const& y, real_t b, 
           stormArray const& x, real_t a) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }

#if 0
  template<class Value>
  Value Dot(Array<Value> const& xVec, 
            Array<Value> const& yVec) {
        
    return ParallelSum(0, xVec.Size(), Value(0),
      [&](size_t i) { return xVec(i)*yVec(i); });
  }

  template<class Value>
  auto Norm2(Array<Value> const& xVec) {
        
    return std::sqrt(Dot(xVec, xVec));
  }

  template<class Value>
  void Set(Array<Value>& xVec, 
           Array<Value> const& yVec) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = yVec(i); });

  } // Set<...>

  template<class Value>
  void Fill(Array<Value>& xVec, 
            Value const& alpha) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = alpha; });

  } // Fill<...>
  template<class Value>
  void RandFill(Array<Value>& xVec) {

    _STORM_NOT_IMPLEMENTED_();

  } // RandFill<...>

  template<class Value>
  void Scale(Array<Value>& xVec, 
             Array<Value> const& yVec, real_t alpha) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = alpha*yVec(i); });

  } // Scale<...>

  template<class Value>
  void Add(Array<Value>& xVec, 
           Array<Value> const& yVec,
           Array<Value> const& zVec) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = yVec(i) + zVec(i); });

  } // Add<...>
  template<class Value>
  void Add(Array<Value>& xVec, 
           Array<Value> const& yVec,
           Array<Value> const& zVec, real_t beta) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = yVec(i) + beta*zVec(i); });

  } // Add<...>
  template<class Value>
  void Add(Array<Value>& xVec, 
           Array<Value> const& yVec, real_t alpha,
           Array<Value> const& zVec, real_t beta) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = alpha*yVec(i) + beta*zVec(i); });

  } // Add<...>

  template<class Value>
  void Sub(Array<Value>& xVec, 
           Array<Value> const& yVec,
           Array<Value> const& zVec) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = yVec(i) - zVec(i); });

  } // Sub<...>
  template<class Value>
  void Sub(Array<Value>& xVec, 
           Array<Value> const& yVec,
           Array<Value> const& zVec, real_t beta) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = yVec(i) - beta*zVec(i); });

  } // Sub<...>
  template<class Value>
  void Sub(Array<Value>& xVec, 
           Array<Value> const& yVec, real_t alpha,
           Array<Value> const& zVec, real_t beta) {

    ParallelFor(0, xVec.Size(), 
      [&](size_t i) { xVec(i) = alpha*yVec(i) - beta*zVec(i); });

  } // Sub<...>
#endif

} // namespace Blas

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator 𝒚 ← 𝓐(𝒙).
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class InVector, class OutVector = InVector>
class Operator : public BaseObject {
public:

  /// @brief Compute an operator-vector product, 𝒚 ← 𝓐(𝒙).
  ///
  /// @param yVec Output vector, 𝒚.
  /// @param xVec Input vector, 𝒙.
  virtual void MatVec(OutVector& yVec,
                      InVector const& xVec) const = 0;

  /// @brief Compute a chained \
  ///   operator-vector product, 𝒛 ← 𝓐(𝒚 ← 𝓑(𝒙)).
  ///
  /// @param zVec Output vector, 𝒛.
  /// @param yVec Intermediate vector, 𝒚.
  /// @param xVec Input vector, 𝒙.
  template<class InOutVector = InVector>
  void MatVec(OutVector& zVec,
              InOutVector& yVec,
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
  void Residual(OutVector& rVec,
                OutVector const& bVec,
                InVector const& xVec) const {
    MatVec(rVec, xVec);
    Blas::Sub(rVec, bVec, rVec);
  }

  /// @brief Compute a residual norm, ‖𝒃 - 𝓐𝒙‖.
  ///
  /// @param bVec Input vector, 𝒃.
  /// @param xVec Input vector, 𝒙.
  real_t ResidualNorm(OutVector const& bVec,
                      InVector const& xVec) const {
    OutVector rVec;
    rVec.Assign(bVec, false);
    Residual(rVec, bVec, xVec);
    return Blas::Norm2(rVec);
  }

  /// @brief Compute an conjugate operator-vector product, 𝒙 ← 𝓐*(𝒚).
  ///
  /// @param xVec Output vector, 𝒙.
  /// @param yVec Input vector, 𝒚.
  virtual void ConjMatVec(InVector& xVec,
                          OutVector const& yVec) const {
    throw std::runtime_error(
      "`Operator<...>::ConjMatVec` was not overriden");
  }

}; // class Operator<...>

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Operator implementation with external function pointers.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class InVector, class OutVector = InVector>
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
  template<class MatVecFunc>
  explicit FunctionalOperator(MatVecFunc&& matVecFunc) :
      MatVecFunc_{std::forward<MatVecFunc>(matVecFunc)} {
    stormAssert(MatVecFunc_);
  }
  template<class MatVecFunc, class ConjMatVecFunc>
  explicit FunctionalOperator(MatVecFunc&& matVecFunc,
                              ConjMatVecFunc&& conjMatVecFunc) :
      MatVecFunc_{std::forward<MatVecFunc>(matVecFunc)},
      ConjMatVecFunc_{std::forward<ConjMatVecFunc>(conjMatVecFunc)} {
    stormAssert(MatVecFunc_ && ConjMatVecFunc_);
  }
  /// @}

private:

  void MatVec(OutVector& yVec,
              InVector const& xVec) const override {
    MatVecFunc_(yVec, xVec);
  }

  void ConjMatVec(InVector& xVec,
                  OutVector const& yVec) const override {
    if (!ConjMatVecFunc_) {
      throw std::runtime_error(
        "`FunctionalOperator<...>::ConjMatVec`"
        " conjugate product function was not set.");
    }
    ConjMatVecFunc_(xVec, yVec);
  }

}; // class FunctionalOperator<...>

/// ----------------------------------------------------------------- ///
/// @brief Make the functional operator.
///  
/// @param matVecFunc Operator-vector product function, 𝒚 ← 𝓐(𝒙).
/// @param conjMatVecFunc Conjugate operator-vector product, 𝒙 ← 𝓐*(𝒚).
/// ----------------------------------------------------------------- ///
/// @{
template<class InVector, class OutVector = InVector, 
         class MatVecFunc>
auto MakeOperator(MatVecFunc&& matVecFunc) {

  return std::make_unique<FunctionalOperator<InVector, OutVector>>(
    std::forward<MatVecFunc>(matVecFunc));

} // MakeOperator<...>
template<class InVector, class OutVector = InVector, 
         class MatVecFunc, class ConjMatVecFunc>
auto MakeOperator(MatVecFunc&& matVecFunc,
                  ConjMatVecFunc&& conjMatVecFunc) {

  return std::make_unique<FunctionalOperator<InVector, OutVector>>(
    std::forward<MatVecFunc>(matVecFunc), 
    std::forward<ConjMatVecFunc>(conjMatVecFunc));

} // MakeOperator<...>
/// @}

/// ----------------------------------------------------------------- ///
/// @brief Make the self-adjoint functional operator.
/// ----------------------------------------------------------------- ///
template<class Vector, class MatVecFunc>
auto MakeSymmetricOperator(MatVecFunc&& matVecFunc) {

  return std::make_unique<FunctionalOperator<Vector>>(
    matVecFunc, std::forward<MatVecFunc>(matVecFunc));

} // MakeSymmetricOperator<...>

} // namespace Storm
