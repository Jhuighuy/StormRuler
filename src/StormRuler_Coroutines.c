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

#include "StormRuler_Coroutines.h"

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#define MINICORO_IMPL
#include "minicoro.h"

void SR_Co_EntryPoint(mco_coro* co) {
  void** userData = (void**)mco_get_user_data(co);
  SR_tCoFunc co_func = (SR_tCoFunc)userData[0];
  void* co_env = userData[1];
  co_func(co_env);
} // void SR_Co_EntryPoint

SR_tCoroutine SR_Co_Create(SR_tCoFunc co_func, void* co_env) {

  void** userData = (void**)calloc(2, sizeof(*userData));
  userData[0] = co_func, userData[1] = co_env;

  mco_desc desc = mco_desc_init(SR_Co_EntryPoint, 0);
  desc.user_data = userData;

  mco_coro* co;
  mco_result res = mco_create(&co, &desc);
  assert(res == MCO_SUCCESS);

  return (SR_tCoroutine)co;
} // SR_Co_Create

void SR_Co_Free(SR_tCoroutine co_) {
  mco_coro* co = (mco_coro*)co_;
  mco_result res = mco_destroy(co);  
  assert(res == MCO_SUCCESS);
}

SR_INTEGER SR_Co_Await(SR_tCoroutine co_) {
  mco_coro* co = (mco_coro*)co_;
  mco_result res = mco_resume(co);
  assert(res == MCO_SUCCESS);
  return 0;
}

void SR_Co_Yield(SR_tCoroutine co_, SR_INTEGER value) {
  mco_coro* co = (mco_coro*)co_;
  mco_result res = mco_yield(co);
  assert(res == MCO_SUCCESS);
}

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#if 0

#include <stdint.h>
#include <stdlib.h>
#include <ucontext.h>

typedef struct thread_t_ thread_t;

static const size_t SR_Co_StackSize = 20*1024*1024; // 20MB.

typedef enum {
  SR_Co_New,
  SR_Co_Running,
  SR_Co_Finished
} SR_eCoState;

#ifdef __x86_64__
typedef union {
  void* ptr;
  uint32_t part[sizeof(void*)/sizeof(uint32_t)];
} SR_uPtrSplitter;
#endif

struct thread_t_ {
  struct {
    ucontext_t callee, caller;
  } coro;
};

struct SR_tCoroutine_t {
  SR_eCoState state;
  SR_tCoFunc function;
  thread_t* thread;
  ucontext_t context;
  char* stack;
  int yield_value;
  void* env;
}; // struct SR_tCoroutine_t

void SR_Co_Yield(SR_tCoroutine coro, int value);

#ifdef __x86_64__
static void SR_Co_EntryPoint(uint32_t part0, uint32_t part1) {
  SR_uPtrSplitter p;
  p.part[0] = part0;
  p.part[1] = part1;
  SR_tCoroutine coro = p.ptr;
  int return_value = coro->function(coro->env);
  coro->state = SR_Co_Finished;
  SR_Co_Yield(coro, return_value);
}
#else
static void SR_Co_EntryPoint(coro_t* coro) {
  int return_value = coro->function(coro->env);
  coro->state = SR_Co_Finished;
  SR_Co_Yield(coro, return_value);
}
#endif

SR_tCoroutine SR_Co_Create(SR_tCoFunc function, void* env) {
  thread_t* thread = calloc(1, sizeof(*thread));

  SR_tCoroutine coroutine = calloc(1, sizeof(*coroutine));

  coroutine->state = SR_Co_New;
  coroutine->stack = calloc(1, SR_Co_StackSize);
  coroutine->thread = thread;
  coroutine->function = function;
  coroutine->env = env;

  getcontext(&coroutine->context);
  coroutine->context.uc_stack.ss_sp = coroutine->stack;
  coroutine->context.uc_stack.ss_size = SR_Co_StackSize;
  coroutine->context.uc_link = 0;

#ifdef __x86_64__
  SR_uPtrSplitter p;
  p.ptr = coroutine;
  makecontext(&coroutine->context,
    (void(*)())SR_Co_EntryPoint, 2, p.part[0], p.part[1]);
#else
  makecontext(&coro->context, 
    (void(*)())SR_Co_EntryPoint, 1, coro);
#endif

  return coroutine;
}

void SR_Co_Free(SR_tCoroutine coroutine) {
  free(coroutine->thread);
  free(coroutine->stack);
  free(coroutine);
} // SR_Co_Free

int SR_Co_Await(SR_tCoroutine coro) {
  if (coro->state == SR_Co_New) {
    coro->state = SR_Co_Running;
  } else if (coro->state == SR_Co_Finished) {
    return 0;
  }

  ucontext_t old_context = coro->thread->coro.caller;
  swapcontext(&coro->thread->coro.caller, &coro->context);
  coro->context = coro->thread->coro.callee;
  coro->thread->coro.caller = old_context;

  return coro->yield_value;
} // SR_Co_Await

void SR_Co_Yield(SR_tCoroutine coroutine, SR_INTEGER value) {
  coroutine->yield_value = value;
  swapcontext(&coroutine->thread->coro.callee, &coroutine->thread->coro.caller);
} // SR_Co_Await

#endif
