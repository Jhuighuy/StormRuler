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
#ifndef _STORM_OPERATOR_HXX_
#define _STORM_OPERATOR_HXX_

#include <memory>
#include <stdexcept>
#include <functional>

#include <StormRuler_API.h>

#define _STORM_NOT_IMPLEMENTED_() do { \
  std::cerr << __FUNCTION__ << " not implemented" << std::endl; exit(1); } while(false)

class stormBaseObject {
public:
  virtual ~stormBaseObject() = default;
};

class stormArray {
public:
  stormMesh_t Mesh = nullptr;
  stormArray_t Array = nullptr;
  std::shared_ptr<int> RefCounter;
private:
  stormReal_t* Data_ = nullptr;
  stormSize_t Size_ = 0;

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

  stormSize_t Size() const noexcept {
    return Size_;
  }
  stormReal_t& operator()(stormSize_t index) {
    stormAssert(index < Size_);
    return Data_[index];
  }
  stormReal_t const& operator()(stormSize_t index) const {
    stormAssert(index < Size_);
    return Data_[index];
  }

  void Assign(stormArray const& like, bool copy = true) {
    this->~stormArray();
    Mesh = like.Mesh;
    Array = stormAllocLike(like.Array);
    stormArrayUnwrap(Array, &Data_, &Size_);
    if (copy) stormSet(Mesh, Array, like.Array);
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
};

namespace stormUtils {
  
  stormReal_t SafeDivide(stormReal_t x, stormReal_t y) {
    return (y == 0.0) ? 0.0 : (x/y);
  }

  stormReal_t& SafeDivideEquals(stormReal_t& x, stormReal_t y) {
    x = SafeDivide(x, y);
    return x;
  }

}

namespace stormBlas {
  stormReal_t Norm2(stormArray const& z) {
    return stormNorm2(z.Mesh, z.Array);
  }
  stormReal_t Dot(stormArray const& z, stormArray const& y) {
    return stormDot(z.Mesh, z.Array, y.Array);
  }

  void Set(stormArray& z, stormArray const& y) {
    stormSet(z.Mesh, z.Array, y.Array);
  }

  void Fill(stormArray& z, stormReal_t a) {
    stormFill(z.Mesh, z.Array, a);
  }
  void RandFill(stormArray& z) {
    stormRandFill(z.Mesh, z.Array);
  }

  void Scale(stormArray& z, stormArray const& y, stormReal_t a) {
    stormScale(z.Mesh, z.Array, y.Array, a);
  }

  void Add(stormArray& z, 
           stormArray const& y, 
           stormArray const& x) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array);
  }
  void Add(stormArray& z, 
           stormArray const& y, 
           stormArray const& x, stormReal_t a) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a);
  }
  void Add(stormArray& z, 
           stormArray const& y, stormReal_t b,
           stormArray const& x, stormReal_t a) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }
  [[deprecated]]
  void Add(stormArray& z, stormArray const& y, stormArray const& x, 
          stormReal_t a, stormReal_t b) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }

  void Sub(stormArray& z, 
           stormArray const& y, 
           stormArray const& x) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array);
  }
  void Sub(stormArray& z, 
           stormArray const& y, 
           stormArray const& x, stormReal_t a) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a);
  }
  void Sub(stormArray& z, 
           stormArray const& y, stormReal_t b, 
           stormArray const& x, stormReal_t a) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }
  [[deprecated]]
  void Sub(stormArray& z, stormArray const& y, stormArray const& x, 
          stormReal_t a, stormReal_t b) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }
}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator 𝒚 ← 𝓐(𝒙).
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class InVector, class OutVector = InVector>
class stormOperator : public stormBaseObject {
public:

  /// @brief Compute an operator-vector product, 𝒚 ← 𝓐(𝒙).
  ///
  /// @param yVec Output vector, 𝒚.
  /// @param xVec Input vector, 𝒙.
  virtual void MatVec(OutVector& yVec,
                      InVector const& xVec) const = 0;

  /// @brief Compute an conjugate operator-vector product, 𝒙 ← 𝓐*(𝒚).
  ///
  /// @param yVec Output vector, 𝒚.
  /// @param xVec Input vector, 𝒙.
  virtual void ConjMatVec(InVector& xVec,
                          OutVector const& yVec) const {
    throw std::runtime_error(
      "`stormOperator<...>::ConjMatVec` was not overriden");
  }

}; // class stormOperator<...>

namespace stormBlas {

  template<class InVector, class tInOutArray, class OutVector>
  void MatVec(OutVector& zVec,
              stormOperator<tInOutArray, OutVector> const& linOp1,
              tInOutArray& yVec,
              stormOperator<InVector, tInOutArray> const& linOp2,
              InVector const& xVec) {

    linOp2.MatVec(yVec, xVec);
    linOp1.MatVec(zVec, yVec);

  } // MatVec<...>

  template<class Vector>
  void ConjMatVec(Vector& yVec,
                  stormOperator<Vector> const& linOp,
                  Vector const& xVec) {

    linOp.ConjMatVec(yVec, xVec);

  } // ConjMatVec

} // namespace stormBlas

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Operator implementation with external function pointers.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class InVector, class OutVector = InVector>
class stormFunctionalOperator final : public stormOperator<InVector, OutVector> {
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
  explicit stormFunctionalOperator(MatVecFunc&& matVecFunc) :
      MatVecFunc_{std::forward<MatVecFunc>(matVecFunc)} {
    stormAssert(MatVecFunc_);
  }
  template<class MatVecFunc, class ConjMatVecFunc>
  explicit stormFunctionalOperator(MatVecFunc&& matVecFunc,
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
        "`stormFunctionalOperator<...>::ConjMatVec` conjugate product function was not set.");
    }
    ConjMatVecFunc_(xVec, yVec);
  }

}; // class stormFunctionalOperator<...>

/// ----------------------------------------------------------------- ///
/// @brief Make the functional operator.
///  
/// @param matVecFunc Operator-vector product function, 𝒚 ← 𝓐(𝒙).
/// @param conjMatVecFunc Conjugate operator-vector product, 𝒙 ← 𝓐*(𝒚).
/// ----------------------------------------------------------------- ///
/// @{
template<class InVector, class OutVector = InVector, 
         class MatVecFunc>
auto stormMakeOperator(MatVecFunc&& matVecFunc) {

  return std::make_unique<stormFunctionalOperator<InVector, OutVector>>(
    std::forward<MatVecFunc>(matVecFunc));

} // stormMakeOperator<...>
template<class InVector, class OutVector = InVector, 
         class MatVecFunc, class ConjMatVecFunc>
auto stormMakeOperator(MatVecFunc&& matVecFunc,
                       ConjMatVecFunc&& conjMatVecFunc) {

  return std::make_unique<stormFunctionalOperator<InVector, OutVector>>(
    std::forward<MatVecFunc>(matVecFunc), 
    std::forward<ConjMatVecFunc>(conjMatVecFunc));

} // stormMakeOperator<...>
/// @}

/// ----------------------------------------------------------------- ///
/// @brief Make the self-adjoint functional operator.
/// ----------------------------------------------------------------- ///
template<class Vector, class MatVecFunc>
auto stormMakeSymmetricOperator(MatVecFunc&& matVecFunc) {

  return std::make_unique<stormFunctionalOperator<Vector>>(
    matVecFunc, std::forward<MatVecFunc>(matVecFunc));

} // stormMakeSymmetricOperator<...>

#endif // ifndef _STORM_OPERATOR_HXX_
