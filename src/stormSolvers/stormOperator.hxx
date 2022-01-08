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

  stormArray() = default;
  stormArray(stormMesh_t mesh, stormArray_t array): Mesh(mesh), Array(array) {
    RefCounter = std::make_shared<int>(2);
  }
  stormArray(stormArray&& oth): Mesh(oth.Mesh), Array(oth.Array) {
    RefCounter = std::move(oth.RefCounter);
    oth.Mesh = nullptr, oth.Array = nullptr, oth.RefCounter = nullptr;
  }
  stormArray(stormArray const& oth): Mesh(oth.Mesh), Array(oth.Array) {
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

  void AllocLike(stormArray const& like, stormArray& z) {
    z.Mesh = like.Mesh;
    z.Array = stormAllocLike(like.Array);
  }
  template<class... tArray>
  void AllocLike(stormArray const& like, stormArray& z, tArray&... zz) {
    AllocLike(like, z);
    AllocLike(like, zz...);
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

  void Add(stormArray& z, stormArray const& y, stormArray const& x, 
          stormReal_t a = 1.0, stormReal_t b = 1.0) {
    stormAdd(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }
  void Sub(stormArray& z, stormArray const& y, stormArray const& x, 
          stormReal_t a = 1.0, stormReal_t b = 1.0) {
    stormSub(z.Mesh, z.Array, y.Array, x.Array, a, b);
  }
}

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Abstract operator 𝒚 ← 𝓐(𝒙).
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tInArray, class tOutArray = tInArray>
class stormOperator : public stormBaseObject {
public:

  /// @brief Apply the operator: 𝒚 ← 𝓐(𝒙).
  ///
  /// @param yArr Output vector, 𝒚.
  /// @param xArr Input vector, 𝒙.
  virtual void MatVec(tOutArray& yArr,
                      tInArray const& xArr) const = 0;

  /// @brief Apply the conjugate operator: 𝒙 ← 𝓐*(𝒚).
  ///
  /// @param yArr Output vector, 𝒚.
  /// @param xArr Input vector, 𝒙.
  virtual void ConjMatVec(tInArray& xArr,
                          tOutArray const& yArr) const {
    throw std::runtime_error(
      "`stormOperator<...>::ConjMatVec` was not overriden");
  }

}; // class stormOperator<...>

namespace stormBlas {

  template<class tInArray, class tInOutArray, class tOutArray>
  void MatVec(tOutArray& zArr,
              stormOperator<tInOutArray, tOutArray> const& linOp1,
              tInOutArray& yArr,
              stormOperator<tInArray, tInOutArray> const& linOp2,
              tInArray const& xArr) {

    linOp2.MatVec(yArr, xArr);
    linOp1.MatVec(zArr, yArr);

  } // MatVec<...>

  template<class tArray>
  void ConjMatVec(tArray& yArr,
                  stormOperator<tArray> const& linOp,
                  tArray const& xArr) {

    linOp.ConjMatVec(yArr, xArr);

  } // ConjMatVec

} // namespace stormBlas

/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// @brief Operator implementation with external function pointers.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
template<class tInArray, class tOutArray = tInArray>
class stormFunctionalOperator final : public stormOperator<tInArray, tOutArray> {
private:
  using tMatVecFunc = std::function<void(tOutArray&, tInArray const&)>;
  using tConjMatVecFunc = std::function<void(tInArray&, tOutArray const&)>;
  tMatVecFunc MatVecPtr;
  tConjMatVecFunc ConjMatVecPtr;

public:

  /// @brief Construct the functional operator.
  ///  
  /// @param matVecPtr 𝒚 ← 𝓐(𝒙) function pointer. 
  /// @param conjMatVecPtr 𝒙 ← 𝓐*(𝒚) function pointer.
  template<class tMatVecFunc>
  explicit stormFunctionalOperator(tMatVecFunc&& matVecPtr) :
      MatVecPtr(std::forward<tMatVecFunc>(matVecPtr)) {
    assert(MatVecPtr);
  }
  template<class tMatVecFunc, class tConjMatVecFunc>
  explicit stormFunctionalOperator(tMatVecFunc&& matVecPtr,
                                   tConjMatVecFunc&& conjMatVecPtr) :
      MatVecPtr(std::forward<tMatVecFunc>(matVecPtr)),
      ConjMatVecPtr(std::forward<tConjMatVecFunc>(conjMatVecPtr)) {
    assert(MatVecPtr && ConjMatVecPtr);
  }

private:
  void MatVec(tOutArray& yArr,
              tInArray const& xArr) const override {
    MatVecPtr(yArr, xArr);
  }

  void ConjMatVec(tInArray& xArr,
                  tOutArray const& yArr) const override {
    if (!ConjMatVecPtr) {
      throw std::runtime_error(
        "`stormFunctionalOperator<...>::ConjMatVec` conjugate product function was not set.");
    }
    ConjMatVecPtr(xArr, yArr);
  }
}; // class stormFunctionalOperator<...>

/// ----------------------------------------------------------------- ///
/// @brief Make the functional operator.
/// ----------------------------------------------------------------- ///
/** @{ */
template<class tInArray, class tOutArray = tInArray, 
         class tMatVecFunc>
std::unique_ptr<stormOperator<tInArray, tOutArray>>
          stormMakeOperator(tMatVecFunc&& matVecPtr) {

  return std::make_unique<stormFunctionalOperator<tInArray, tOutArray>>(
    std::forward<tMatVecFunc>(matVecPtr));

} // stormMakeOperator<...>
template<class tInArray, class tOutArray = tInArray, 
         class tMatVecFunc, class tConjMatVecFunc>
std::unique_ptr<stormOperator<tInArray, tOutArray>>
          stormMakeOperator(tMatVecFunc&& matVecPtr,
                            tConjMatVecFunc&& conjMatVecPtr) {

  return std::make_unique<stormFunctionalOperator<tInArray, tOutArray>>(
    std::forward<tMatVecFunc>(matVecPtr), 
    std::forward<tConjMatVecFunc>(conjMatVecPtr));

} // stormMakeOperator<...>
/** @} */

/// ----------------------------------------------------------------- ///
/// @brief Make the self-adjoint functional operator.
/// ----------------------------------------------------------------- ///
template<class tArray, class tMatVecFunc>
std::unique_ptr<stormOperator<tArray>> 
    stormMakeSymmetricOperator(tMatVecFunc&& matVecPtr) {

  return std::make_unique<stormFunctionalOperator<tArray>>(
    matVecPtr, std::forward<tMatVecFunc>(matVecPtr));

} // stormMakeSymmetricOperator<...>

#endif // ifndef _STORM_OPERATOR_HXX_
