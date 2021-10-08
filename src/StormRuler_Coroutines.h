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

#include "StormRuler_Params.h"

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

/**
 * Couroutine opaque handle.
 */
SR_OPAQUE_STRUCT(SR_tCoroutine);

/**
 * Coroutine function pointer.
 * 
 * @param[in] co Couroutine handle.
 * @param[inout] co_env Couroutine environment.
 * 
 * @returns Coroutine return value.
 */
typedef SR_INTEGER(*SR_tCoFunc)(SR_tCoroutine co, void* co_env);

/**
 * Create a new coroutine with the specified
 * entry point and environment.
 * 
 * Couroutine is created in the suspended state.
 * 
 * @param[in] co_func Couroutine function.
 * @param[inout] co_env Couroutine environment.
 * 
 * @returns Couroutine handle.
 */
SR_API SR_tCoroutine SR_Co_Create(SR_tCoFunc co_func, void* co_env);

/**
 * Delete the couroutine.
 * 
 * @param[in] co Couroutine handle.
 */
SR_API void SR_Co_Free(SR_tCoroutine co);

/**
 * Switch context to the couroutine and wait for it to yield.
 * 
 * @param[in] co Couroutine handle.
 * 
 * @returns Couroutine yeild value.
 */
SR_API SR_INTEGER SR_Co_Await(SR_tCoroutine co);

/**
 * Yeild with a value to the caller context.
 *  
 * @param[in] co Couroutine handle.
 * @param[in] value Value passed to the caller context.
 */
SR_API void SR_Co_Yield(SR_tCoroutine co, SR_INTEGER value);
