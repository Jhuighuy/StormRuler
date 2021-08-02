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
    
/**
 * Scalar/vector/tensor field class.
 */
template<int rank>
class tField {
private:
  void* mData;
public:
  explicit tField(void* data) :
    mData(data) {
  }
  void* Data() {
    return mData;
  }
}; // class tField

using tMFunc = void(*)(
  int* shape, double* in, double* out, void* env);
using tSMFunc = void(*)(
  int dim, double* x, int* shape, double* in, double* out, void* env);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

// Scaled field wrapper.
template<int rank>
struct tScaledField {
  double Scale;
  tField<rank> Field;
  explicit tScaledField(double scale, tField<rank>& field) :
    Scale(scale), Field(field) {
  }
}; // struct tScaledField

// Scale a field.
template<int rank>
auto operator*(double scale, tField<rank> field) {
  return tScaledField<rank>(scale, field);
} 
template<int rank>
auto operator*(tField<rank> field, double scale) {
  return tScaledField<rank>(scale, field);
}
template<int rank>
auto operator/(tField<rank> field, double scale) {
  assert(scale != 0.0);
  return tScaledField<rank>(1.0/scale, field);
}

// Rescale a field.
template<int rank>
auto operator*(double scale, tScaledField<rank> field) {
  return tScaledField<rank>(scale*field.Scale, field.Field);
} 
template<int rank>
auto operator*(tScaledField<rank> field, double scale) {
  return tScaledField<rank>(scale*field.Scale, field.Field);
}
template<int rank>
auto operator/(tScaledField<rank> field, double scale) {
  return tScaledField<rank>(field.Scale/scale, field.Field);
}

// Negate a field.

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

// Field sum wrapper.
template<int rank>
struct tFieldSum {
  tScaledField<rank> LeftField;
  tScaledField<rank> RightField;
  explicit tFieldSum(tScaledField<rank> leftField,
                     tScaledField<rank> rightField) :
    LeftField(leftField), RightField(rightField) {
  }
}; // struct tFieldSum

// Add two fields.
template<int rank>
auto operator+(tScaledField<rank> leftField,
               tScaledField<rank> rightField) {
  return tFieldSum<rank>(leftField, rightField);
}
template<int rank>
auto operator+(tField<rank> leftField, 
               tScaledField<rank> rightField) {
  return tFieldSum<rank>(tScaledField<rank>(1.0, leftField), rightField);
}
template<int rank>
auto operator+(tScaledField<rank> leftField, 
               tField<rank> rightField) {
  return tFieldSum<rank>(leftField, tScaledField<rank>(1.0, rightField));
}
template<int rank>
auto operator+(tField<rank> leftField, tField<rank> rightField) {
  return tFieldSum<rank>(leftField, rightField);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

// Field difference wrapper.
template<int rank>
struct tFieldDiff {
  tScaledField<rank> LeftField;
  tScaledField<rank> RightField;
  explicit tFieldDiff(tScaledField<rank> leftField,
                      tScaledField<rank> rightField) :
    LeftField(leftField), RightField(rightField) {
  }
}; // struct tFieldDiff

// Add two fields.
template<int rank>
auto operator-(tScaledField<rank> leftField,
               tScaledField<rank> rightField) {
  return tFieldDiff<rank>(leftField, rightField);
}
template<int rank>
auto operator-(tField<rank> leftField, 
               tScaledField<rank> rightField) {
  return tFieldDiff<rank>(tScaledField<rank>(1.0, leftField), rightField);
}
template<int rank>
auto operator-(tScaledField<rank> leftField, 
               tField<rank> rightField) {
  return tFieldDiff<rank>(leftField, tScaledField<rank>(1.0, rightField));
}
template<int rank>
auto operator-(tField<rank> leftField, tField<rank> rightField) {
  return tFieldDiff<rank>(leftField, rightField);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //



} // namespace StormRuler
