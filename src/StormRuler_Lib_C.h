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

#pragma once
#ifndef StormRuler_Lib_C_INCLUDED_
#define StormRuler_Lib_C_INCLUDED_

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/**
 * Tensor element type.
 */
typedef enum {
  StormRuler_Type_Real = 'd',    /** Real type, @c double. */
  StormRuler_Type_Complex = 'z', /** Real type, @c complex. */
  StormRuler_Type_Symbol = 's',  /** Symbolic element type, @c StormRuler_Symbol. */
} StormRuler_Type;

/**
 * Tensor function pointer (typeless).
 *
 * @param[in] type_x Type of the tensor arguments.
 * @param[in] rank_x Rank of the tensor arguments. 
 * @param[in] shape_x Shape of the tensor arguments, size is @p rank_x.
 * @param[out] Mx Pointer to the output tensor, size is product of @p shape_x elements.
 * @param[in] x Pointer to the input tensor, size is product of @p shape_x elements.
 * 
 * @param[in] env User environment.  
 */
typedef void(*StormRuler_MapFuncPtr)(StormRuler_Type type_x,
                                     int rank_x, const int* pShape_x, 
                                     void* pMx, const void* pX, const void* env);

/**
 * Spatial tensor function pointer (typeless).
 * 
 * @param[in] dim_r Number of spatial dimensions.
 * @param[in] pR Pointer to the spatial coordinate vector, size is @p dim_r.
 * 
 * @param[in] type_x Type of the tensor arguments.
 * @param[in] rank_x Rank of the tensor arguments. 
 * @param[in] shape_x Shape of the tensor arguments, size is @p rank_x.
 * @param[out] Mx Pointer to the output tensor, size is product of @p shape_x elements.
 * @param[in] x Pointer to the input tensor, size is product of @p shape_x elements.
 * 
 * @param[in] env User environment.  
 */
typedef void(*StormRuler_SMapFuncPtr)(int dim_r, const double* pR,
                                      StormRuler_Type type_x, 
                                      int rank_x, const int* pShape_x, 
                                      void* pMx, const void* pX, const void* env);

/**
 * @brief Field domain handle.
 * 
 * Field domain is responsible for field allocation,
 * running computational kernels, ...
 */
typedef struct StormRuler_Domain StormRuler_Domain;

/**
 * Tensor array handle (typeless).
 */
typedef struct StormRuler_Field StormRuler_Field;

/**
 * Matrix-vector product function (typeless).
 * 
 * @param[inout] domain Domain handle.
 * @param[inout] Ax Result field.
 * @param[in] x Input field.
 */
typedef void(*StormRuler_MatVecFuncPtr)(StormRuler_Domain domain,
                                        StormRuler_Field Ax, StormRuler_Field x);

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#endif // StormRuler_Lib_C_INCLUDED_
