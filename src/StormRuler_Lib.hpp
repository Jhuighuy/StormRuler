// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// Copyright (C) 2021 Oleg Butakov
// 
// Permission is hereby granted, free of charge, to any person 
// obtaining a copy of this software and associated documentation 
// files (the "Software"), to deal in the Software without 
// restriction, including without limitation the rights  to use, 
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the  
// Software is furnished to do so, subject to the following 
// conditions:
// 
// The above copyright notice and this permission notice shall be 
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#include <cassert>
#include <functional>

namespace StormRuler {

class tFieldBase;

template<int rank> 
void DeallocateField(tFieldBase*);

/**
 * Scalar/vector/tensor field class.
 */
template<int rank>
class tField {
private:
  tFieldBase* mData;
  int* mRefCounter;

public:
  explicit tField(tFieldBase* data) :
    mData(data), mRefCounter(new int(1)) {
  }
  explicit tField(tFieldBase* data, std::nullptr_t) :
    mData(data), mRefCounter(nullptr) {
  }
  ~tField() {
    if (mRefCounter != nullptr) {
      *mRefCounter -= 1;
      if (*mRefCounter == 0) {
        DeallocateField<rank>(mData); mData = nullptr; 
        delete mRefCounter; mRefCounter = nullptr;
      }
    }
  }

  tField(const tField& other) :
    mData(other.mData), mRefCounter(other.mRefCounter) {
    if (mRefCounter != nullptr) {
      *mRefCounter += 1;
    }
  }
  tField& operator=(const tField& other) {
    this->~tField();
    new (this) tField(other);
    return *this;
  }
  
  tFieldBase* Data() {
    return mData;
  }
}; // class tField

// Basic math function pointer (for Lib API).
using tMFuncPtr = void(*)(
  int* shape, double* in, double* out, void* env);

// Basic spatial math function pointer (for Lib API).
using tSMFuncPtr = void(*)(
  int dim, double* x, int* shape, double* in, double* out, void* env);

// Basic mesh operator function pointer (for Lib API).
using tMeshOperatorPtr = void(*)(tFieldBase* out, tFieldBase* in, void* env);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

} // namespace StormRuler

// Include the Library API.
#include "StormRuler_Lib.inc"

namespace StormRuler {

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#define _TI_TYPENAME(tExpr) template<int> typename tExpr

struct _tFieldExpr {
  static constexpr bool cHasScaledLikeBody = false;
  static constexpr bool cHasAssignEval = false;
};  // struct _tFieldExpr

template<int rank>
void operator<<(tField<rank> target, double right) {
  BLAS_Fill(target, right);
}
template<int rank>
void operator<<(tField<rank> target, tField<rank> right) {
  BLAS_Set(target, right);
}
// Evaluate expression assignment with temporary storage.
template<int rank, _TI_TYPENAME(tRightExpr)>
void operator<<(tField<rank> target, 
                std::pair<tField<rank>, tRightExpr<rank>> rightExpr) {
  auto [tmp, right] = rightExpr;
  tmp << right;
  target << tmp;
}

#define _SCALED_EXPR_BODY(tScaledFieldExpr, rank) \
  static constexpr bool cHasScaledLikeBody = true; \
  double Scale; \
  tField<rank> Field; \
  \
  explicit tScaledFieldExpr(double scale, tField<rank> field) : \
    Scale(scale), Field(field) { \
  }

template<int rank>
struct tScaledFieldExpr : public _tFieldExpr {
  static constexpr bool cHasAssignEval = false;
  _SCALED_EXPR_BODY(tScaledFieldExpr, rank)
}; // struct tScaledFieldExpr

template<int rank>
void operator<<(tField<rank> left,
                tScaledFieldExpr<rank> right) {
  abort();
}

template<int rank>
auto operator*(double left, tField<rank> right) {
  return tScaledFieldExpr<rank>(left, right);
} 
template<int rank>
auto operator*(tField<rank> left, double right) {
  return tScaledFieldExpr<rank>(right, left);
}

template<int rank, _TI_TYPENAME(tRightExpr)>
auto operator*(double left, 
               tRightExpr<rank> rightExpr) {
  static_assert(tRightExpr<rank>::cHasScaledLikeBody);
  return rightExpr.Scale *= left, rightExpr;
} 
template<int rank, _TI_TYPENAME(tLeftExpr)>
auto operator*(tLeftExpr<rank> leftExpr, double right) {
  static_assert(tLeftExpr<rank>::cHasScaledLikeBody);
  return leftExpr.Scale *= right, leftExpr;
}

// Negate a field.
// TODO:

template<int rank>
struct tFieldSumExpr : public _tFieldExpr {
  static constexpr bool cHasAssignEval = false;
  tScaledFieldExpr<rank> LeftExpr;
  tScaledFieldExpr<rank> RightExpr;

  explicit tFieldSumExpr(tScaledFieldExpr<rank> left,
                         tScaledFieldExpr<rank> right) :
    LeftExpr(left), RightExpr(right) {
  }
}; // struct tFieldSumExpr

template<int rank>
void operator<<(tField<rank> target,
                tFieldSumExpr<rank> rightExpr) {
  BLAS_Add(target, 
           rightExpr.LeftExpr.Field, rightExpr.RightExpr.Field,
           rightExpr.RightExpr.Scale, rightExpr.LeftExpr.Scale);
}

template<int rank>
auto operator+(tScaledFieldExpr<rank> leftExpr,
               tScaledFieldExpr<rank> rightExpr) {
  return tFieldSumExpr<rank>(leftExpr, rightExpr);
}
template<int rank>
auto operator+(tField<rank> left, 
               tScaledFieldExpr<rank> rightExpr) {
  return tFieldSumExpr<rank>(
    tScaledFieldExpr<rank>(1.0, left), rightExpr);
}
template<int rank>
auto operator+(tScaledFieldExpr<rank> leftExpr, 
               tField<rank> right) {
  return tFieldSumExpr<rank>(
    leftExpr, tScaledFieldExpr<rank>(1.0, right));
}
template<int rank>
auto operator+(tField<rank> left, tField<rank> right) {
  return tFieldSumExpr<rank>(
    tScaledFieldExpr<rank>(1.0, left), 
    tScaledFieldExpr<rank>(1.0, right));
}

template<int rank>
struct tFieldDiffExpr : public _tFieldExpr {
  static constexpr bool cHasAssignEval = false;
  tScaledFieldExpr<rank> LeftExpr;
  tScaledFieldExpr<rank> RightExpr;

  explicit tFieldDiffExpr(tScaledFieldExpr<rank> left,
                          tScaledFieldExpr<rank> right) :
    LeftExpr(left), RightExpr(right) {
  }
}; // struct tFieldDiffExpr

template<int rank>
void operator<<(tField<rank> target,
                tFieldDiffExpr<rank> rightExpr) {
  BLAS_Sub(target, 
           rightExpr.LeftExpr.Field, rightExpr.RightExpr.Field,
           rightExpr.RightExpr.Scale, rightExpr.LeftExpr.Scale);
}

template<int rank>
auto operator-(tScaledFieldExpr<rank> leftExpr,
               tScaledFieldExpr<rank> rightExpr) {
  return tFieldDiffExpr<rank>(leftExpr, rightExpr);
}
template<int rank>
auto operator-(tField<rank> left, 
               tScaledFieldExpr<rank> rightExpr) {
  return tFieldDiffExpr<rank>(
    tScaledFieldExpr<rank>(1.0, left), rightExpr);
}
template<int rank>
auto operator-(tScaledFieldExpr<rank> leftExpr, 
               tField<rank> right) {
  return tFieldDiffExpr<rank>(
    leftExpr, tScaledFieldExpr<rank>(1.0, right));
}
template<int rank>
auto operator-(tField<rank> left, tField<rank> right) {
  return tFieldDiffExpr<rank>(
    tScaledFieldExpr<rank>(1.0, left), 
    tScaledFieldExpr<rank>(1.0, right));
}

template<int rank, typename tMFunc>
struct tFieldMapExpr : public _tFieldExpr {
  static constexpr bool cHasAssignEval = false;
  tMFunc MFunc;
  tField<rank> Field;

  explicit tFieldMapExpr(tMFunc&& func, tField<rank> field) :
    MFunc(std::forward<tMFunc>(func)), Field(field) {
  }
}; // struct tFieldDiffExpr

template<int rank, typename tMFunc>
auto MAP(tMFunc&& func, tField<rank> field) {
  return tFieldMapExpr<rank, tMFunc>(std::forward<tMFunc>(func), field);
}

template<int rank, typename tMFunc>
void operator<<(tField<rank> target,
                tFieldMapExpr<rank, tMFunc> rightExpr) {
  BLAS_FuncProd(target, rightExpr.Field, std::move(rightExpr.MFunc));
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

template<int rank>
struct tFieldGradientExpr : public _tFieldExpr {
  static constexpr bool cHasAssignEval = true;
  _SCALED_EXPR_BODY(tFieldGradientExpr, rank-1)
}; // struct tFieldGradientExpr
  
template<int rank>
void operator-=(tField<rank> left, 
                tFieldGradientExpr<rank> rightExpr) {
  FDM_Gradient(left, +rightExpr.Scale, rightExpr.Field, 'c');
}
template<int rank>
void operator+=(tField<rank> left, 
                tFieldGradientExpr<rank> rightExpr) {
  FDM_Gradient(left, -rightExpr.Scale, rightExpr.Field, 'c');
}

template<int rank>
auto GRAD(tScaledFieldExpr<rank> expr) {
  return tFieldGradientExpr<rank+1>(expr);
}
template<int rank>
auto GRAD(tField<rank> field) {
  return tFieldGradientExpr<rank+1>(1.0, field);
}

template<int rank>
struct tFieldDivergenceExpr : public _tFieldExpr {
  static constexpr bool cHasAssignEval = true;
  _SCALED_EXPR_BODY(tFieldDivergenceExpr, rank+1)
}; // struct tFieldDivergenceExpr

template<int rank>
void operator-=(tField<rank> left, 
                tFieldDivergenceExpr<rank> rightExpr) {
  FDM_Divergence(left, +rightExpr.Scale, rightExpr.Field, 'c');
}
template<int rank>
void operator+=(tField<rank> left, 
                tFieldDivergenceExpr<rank> rightExpr) {
  FDM_Divergence(left, -rightExpr.Scale, rightExpr.Field, 'c');
}

template<int rank>
auto DIV(tScaledFieldExpr<rank> expr) {
  return tFieldDivergenceExpr<rank-1>(expr);
}
template<int rank>
auto DIV(tField<rank> field) {
  return tFieldDivergenceExpr<rank-1>(1.0, field);
}

template<int rank>
struct tFieldLaplacianExpr : public _tFieldExpr {
  static constexpr bool cHasAssignEval = true;
  _SCALED_EXPR_BODY(tFieldLaplacianExpr, rank)
}; // struct tFieldLaplacianExpr

template<int rank>
void operator+=(tField<rank> left, 
                tFieldLaplacianExpr<rank> rightExpr) {
  FDM_Laplacian(left, +rightExpr.Scale, rightExpr.Field);
}
template<int rank>
void operator-=(tField<rank> left, 
                tFieldLaplacianExpr<rank> rightExpr) {
  FDM_Laplacian(left, -rightExpr.Scale, rightExpr.Field);
}

template<int rank>
auto DIVGRAD(tScaledFieldExpr<rank> expr) {
  return tFieldLaplacianExpr<rank>(expr);
}
template<int rank>
auto DIVGRAD(tField<rank> field) {
  return tFieldLaplacianExpr<rank>(1.0, field);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

template<int rank>
struct tConvectionExpr : public _tFieldExpr {
  static constexpr bool cHasAssignEval = true;
  static constexpr bool cHasScaledLikeBody = true;
  double Scale;
  tField<rank> Field;
  tField<1> VelocityField;

  explicit tConvectionExpr(double scale, 
                           tField<rank> field, tField<1> velocityField) :
    Scale(scale), Field(field),VelocityField(velocityField) {
  }
}; // struct tConvectionExpr

template<int rank>
void operator-=(tField<rank> left, 
                tConvectionExpr<rank> rightExpr) {
  FDM_Convection(
    left, +rightExpr.Scale, rightExpr.Field, rightExpr.VelocityField);
}
template<int rank>
void operator+=(tField<rank> left, 
                tConvectionExpr<rank> rightExpr) {
  FDM_Convection(
    left, -rightExpr.Scale, rightExpr.Field, rightExpr.VelocityField);
}

template<int rank>
auto CONV(tField<rank> field, tField<1> velocityField) {
  return tConvectionExpr<rank>(1.0, field, velocityField);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

// Initiate tmp field allocation.
template<int rank, _TI_TYPENAME(tLeftExpr), _TI_TYPENAME(tRightExpr)>
auto operator+(tLeftExpr<rank> left, tRightExpr<rank> right) {
  auto tmp = AllocateField<rank>();
  return std::make_pair(tmp, left) + right;
}
template<int rank, _TI_TYPENAME(tLeftExpr), _TI_TYPENAME(tRightExpr)>
auto operator-(tLeftExpr<rank> left, tRightExpr<rank> right) {
  auto tmp = AllocateField<rank>();
  return std::make_pair(tmp, left) - right;
}

template<int rank, _TI_TYPENAME(tLeftExpr), typename tRightExpr>
auto operator+(std::pair<tField<rank>, tLeftExpr<rank>> leftExpr, tRightExpr rightExpr) {
  auto [tmp, left] = leftExpr;
  tmp << left;
  if constexpr (tRightExpr::cHasAssignEval) {
    tmp += rightExpr;
    return std::make_pair(tmp, tmp);
  } else {
    return std::make_pair(tmp, tmp + rightExpr);
  }
}
template<int rank, _TI_TYPENAME(tLeftExpr), typename tRightExpr>
auto operator-(std::pair<tField<rank>, tLeftExpr<rank>> leftExpr, tRightExpr rightExpr) {
  auto [tmp, left] = leftExpr;
  tmp << left;
  if constexpr (tRightExpr::cHasAssignEval) {
    tmp -= rightExpr;
    return std::make_pair(tmp, tmp);
  } else {
    return std::make_pair(tmp, tmp - rightExpr);
  }
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

template<int rank, typename tMeshOperator>
void SOLVE_BiCGSTAB(tMeshOperator&& meshOperator, 
                    tField<rank> solution, tField<rank> rhs) {
  Solve_BiCGStab(solution, rhs, std::forward<tMeshOperator>(meshOperator));
}

} // namespace StormRuler
