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

#include "StormRuler_API.h"
#include "StormRuler_Coroutines.h"

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

typedef struct {
  char type;
  SR_tMesh mesh;
  SR_tFieldA x, b;
  SR_eSolver Solver; 
  SR_ePrecond Precond;
  SR_tFieldA Ay, y;
  SR_tCoroutine coroutine;
} SR_RCI_tLinSolveA_Env;

static void SR_RCI_MatVecA_Co(SR_tMesh mesh,
    SR_tFieldA Ay, SR_tFieldA y, void* co_env_) {
  SR_RCI_tLinSolveA_Env* co_env = (SR_RCI_tLinSolveA_Env*)co_env_;

  puts("hui!!!"); 

  // Pass field through the environment
  // and yield to the main coroutine to compute the matrix-vector product.
  co_env->Ay = Ay, co_env->y = y;
  SR_Co_Yield(co_env->coroutine, SR_Request_MatVec);
} // SR_RCI_MatVecA_Co

static void SR_RCI_LinSolveA_Co(void* co_env_) {
  SR_RCI_tLinSolveA_Env* co_env = (SR_RCI_tLinSolveA_Env*)co_env_;

  SR_LinSolveR(co_env->mesh, 
    co_env->x.R, co_env->b.R, 
    (SR_tMatVecFuncR)SR_RCI_MatVecA_Co, co_env, 
    co_env->Solver, co_env->Precond, NULL, NULL);

  co_env->Ay = SR_NULL_A, co_env->y = SR_NULL_A;
  SR_Co_Return(co_env->coroutine, SR_Done);
} // SR_RCI_LinSolveR_Co

SR_API SR_eRequest SR_RCI_LinSolveR(SR_tMesh mesh,
    SR_tFieldR x, SR_tFieldR b, 
    SR_eSolver solver, SR_ePrecond precond,
    SR_tFieldR* pAy, SR_tFieldR* pY) {

  static SR_RCI_tLinSolveA_Env co_env = {};
  if (co_env.coroutine == NULL) {
    // Set the environment and create the compute coroutine.
    co_env = (SR_RCI_tLinSolveA_Env){'R', mesh, x, b, solver, precond};
    co_env.coroutine = SR_Co_Create(SR_RCI_LinSolveA_Co, &co_env); 
  }

  // Switch to the compute coroutine 
  // and wait for it to yield or return.
  SR_Co_Await(co_env.coroutine);
  *pAy = co_env.Ay.R, *pY = co_env.y.R;
  if (*pAy == NULL) {
    SR_Co_Free(co_env.coroutine), co_env.coroutine = NULL;
    return SR_Done;
  }
  return SR_Request_MatVec; 
} // SR_RCI_LinSolveR
