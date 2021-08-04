// 
// THIS IS THE AUTO-GENERATED FILE, DO NOT EDIT MANUALLY 
// 
#define EXTERN extern "C" 
namespace StormRuler { 

// 
// EXPORT LIB MESH & FIELDS ALLOCATION & FIELDS IO 
//
EXTERN void Lib_InitializeMesh();

template<int rank> 
tField<rank> AllocateField();
EXTERN void _Lib_AllocateField0(tFieldBase**); 
template<> 
tField<0> AllocateField<0>() { 
  tFieldBase* pData; 
  _Lib_AllocateField0(&pData); 
  return tField<0>(pData); 
}
EXTERN void _Lib_AllocateField1(tFieldBase**); 
template<> 
tField<1> AllocateField<1>() { 
  tFieldBase* pData; 
  _Lib_AllocateField1(&pData); 
  return tField<1>(pData); 
}
EXTERN void _Lib_AllocateField2(tFieldBase**); 
template<> 
tField<2> AllocateField<2>() { 
  tFieldBase* pData; 
  _Lib_AllocateField2(&pData); 
  return tField<2>(pData); 
}

template<int rank> 
void DeallocateField(tFieldBase*);
EXTERN void _Lib_DeallocateField0(tFieldBase*); 
template<> 
void DeallocateField<0>(tFieldBase* pData) { 
  _Lib_DeallocateField0(pData); 
}
EXTERN void _Lib_DeallocateField1(tFieldBase*); 
template<> 
void DeallocateField<1>(tFieldBase* pData) { 
  _Lib_DeallocateField1(pData); 
}
EXTERN void _Lib_DeallocateField2(tFieldBase*); 
template<> 
void DeallocateField<2>(tFieldBase* pData) { 
  _Lib_DeallocateField2(pData); 
}

EXTERN void _Lib_IO_Begin();

EXTERN void _Lib_IO_Add0(tFieldBase*, const char*, int); 
void _Lib_IO_Add(tField<0> u, const std::string& name) { 
  _Lib_IO_Add0(u.Data(), name.c_str(), name.size()); 
}
EXTERN void _Lib_IO_Add1(tFieldBase*, const char*, int); 
void _Lib_IO_Add(tField<1> u, const std::string& name) { 
  _Lib_IO_Add1(u.Data(), name.c_str(), name.size()); 
}
EXTERN void _Lib_IO_Add2(tFieldBase*, const char*, int); 
void _Lib_IO_Add(tField<2> u, const std::string& name) { 
  _Lib_IO_Add2(u.Data(), name.c_str(), name.size()); 
}

EXTERN void _Lib_IO_End();


// 
// EXPORT LIB BLAS 
//
EXTERN void _Lib_BLAS_Fill0(tFieldBase*, double); 
void BLAS_Fill(tField<0> u, double alpha) { 
  _Lib_BLAS_Fill0(u.Data(), alpha); 
}
EXTERN void _Lib_BLAS_Fill1(tFieldBase*, double); 
void BLAS_Fill(tField<1> u, double alpha) { 
  _Lib_BLAS_Fill1(u.Data(), alpha); 
}
EXTERN void _Lib_BLAS_Fill2(tFieldBase*, double); 
void BLAS_Fill(tField<2> u, double alpha) { 
  _Lib_BLAS_Fill2(u.Data(), alpha); 
}

EXTERN void _Lib_BLAS_Set0(tFieldBase*, tFieldBase*); 
void BLAS_Set(tField<0> u, tField<0> v) { 
  _Lib_BLAS_Set0(u.Data(), v.Data()); 
}
EXTERN void _Lib_BLAS_Set1(tFieldBase*, tFieldBase*); 
void BLAS_Set(tField<1> u, tField<1> v) { 
  _Lib_BLAS_Set1(u.Data(), v.Data()); 
}
EXTERN void _Lib_BLAS_Set2(tFieldBase*, tFieldBase*); 
void BLAS_Set(tField<2> u, tField<2> v) { 
  _Lib_BLAS_Set2(u.Data(), v.Data()); 
}

EXTERN void _Lib_BLAS_Add0( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Add(tField<0> u, tField<0> v, 
              tField<0> w, double alpha, double beta) { 
  _Lib_BLAS_Add0(u.Data(), v.Data(), w.Data(), alpha, beta); 
}
EXTERN void _Lib_BLAS_Add1( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Add(tField<1> u, tField<1> v, 
              tField<1> w, double alpha, double beta) { 
  _Lib_BLAS_Add1(u.Data(), v.Data(), w.Data(), alpha, beta); 
}
EXTERN void _Lib_BLAS_Add2( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Add(tField<2> u, tField<2> v, 
              tField<2> w, double alpha, double beta) { 
  _Lib_BLAS_Add2(u.Data(), v.Data(), w.Data(), alpha, beta); 
}

EXTERN void _Lib_BLAS_Sub0( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Sub(tField<0> u, tField<0> v, 
              tField<0> w, double alpha, double beta) { 
  _Lib_BLAS_Sub0(u.Data(), v.Data(), w.Data(), alpha, beta); 
}
EXTERN void _Lib_BLAS_Sub1( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Sub(tField<1> u, tField<1> v, 
              tField<1> w, double alpha, double beta) { 
  _Lib_BLAS_Sub1(u.Data(), v.Data(), w.Data(), alpha, beta); 
}
EXTERN void _Lib_BLAS_Sub2( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Sub(tField<2> u, tField<2> v, 
              tField<2> w, double alpha, double beta) { 
  _Lib_BLAS_Sub2(u.Data(), v.Data(), w.Data(), alpha, beta); 
}

EXTERN void _Lib_BLAS_Mul0(tFieldBase*, tFieldBase*, tFieldBase*); 
void BLAS_Mul(tField<0> u, 
              tField<0> v, tField<0> w) { 
  _Lib_BLAS_Mul0(u.Data(), v.Data(), w.Data()); 
}
EXTERN void _Lib_BLAS_Mul1(tFieldBase*, tFieldBase*, tFieldBase*); 
void BLAS_Mul(tField<1> u, 
              tField<0> v, tField<1> w) { 
  _Lib_BLAS_Mul1(u.Data(), v.Data(), w.Data()); 
}
EXTERN void _Lib_BLAS_Mul2(tFieldBase*, tFieldBase*, tFieldBase*); 
void BLAS_Mul(tField<2> u, 
              tField<0> v, tField<2> w) { 
  _Lib_BLAS_Mul2(u.Data(), v.Data(), w.Data()); 
}

EXTERN void _Lib_BLAS_FuncProd0(tFieldBase*, tFieldBase*, tMFuncPtr, void*); 
template<typename tMFunc> 
void BLAS_FuncProd(tField<0> v, 
                   tField<0> u, tMFunc&& func) { 
  _Lib_BLAS_FuncProd0(v.Data(), u.Data(), 
    [](int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tMFunc*>(env); 
      *out = func(*in);
    }, &func); 
}
EXTERN void _Lib_BLAS_FuncProd1(tFieldBase*, tFieldBase*, tMFuncPtr, void*); 
template<typename tMFunc> 
void BLAS_FuncProd(tField<1> v, 
                   tField<1> u, tMFunc&& func) { 
  _Lib_BLAS_FuncProd1(v.Data(), u.Data(), 
    [](int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tMFunc*>(env); 
      func(in, out);
    }, &func); 
}
EXTERN void _Lib_BLAS_FuncProd2(tFieldBase*, tFieldBase*, tMFuncPtr, void*); 
template<typename tMFunc> 
void BLAS_FuncProd(tField<2> v, 
                   tField<2> u, tMFunc&& func) { 
  _Lib_BLAS_FuncProd2(v.Data(), u.Data(), 
    [](int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tMFunc*>(env); 
      func(in, out);
    }, &func); 
}

EXTERN void _Lib_BLAS_SFuncProd0(tFieldBase*, tFieldBase*, tSMFuncPtr, void*); 
template<typename tSMFunc> 
void BLAS_SFuncProd(tField<0> v, 
                    tField<0> u, tSMFunc&& func) { 
  _Lib_BLAS_SFuncProd0(v.Data(), u.Data(), 
    [](int dim, double* x, int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tSMFunc*>(env); 
      func(x, in, out); 
    }, &func); 
}
EXTERN void _Lib_BLAS_SFuncProd1(tFieldBase*, tFieldBase*, tSMFuncPtr, void*); 
template<typename tSMFunc> 
void BLAS_SFuncProd(tField<1> v, 
                    tField<1> u, tSMFunc&& func) { 
  _Lib_BLAS_SFuncProd1(v.Data(), u.Data(), 
    [](int dim, double* x, int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tSMFunc*>(env); 
      func(x, in, out); 
    }, &func); 
}
EXTERN void _Lib_BLAS_SFuncProd2(tFieldBase*, tFieldBase*, tSMFuncPtr, void*); 
template<typename tSMFunc> 
void BLAS_SFuncProd(tField<2> v, 
                    tField<2> u, tSMFunc&& func) { 
  _Lib_BLAS_SFuncProd2(v.Data(), u.Data(), 
    [](int dim, double* x, int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tSMFunc*>(env); 
      func(x, in, out); 
    }, &func); 
}


// 
// EXPORT LIB FDM OPERATORS 
//
EXTERN void _Lib_FDM_Gradient0( 
    tFieldBase*, double, tFieldBase*, char); 
void FDM_Gradient(tField<1> vBar, 
                  double lambda, tField<0> u, char dir) { 
  _Lib_FDM_Gradient0(vBar.Data(), lambda, u.Data(), dir); 
}
EXTERN void _Lib_FDM_Gradient1( 
    tFieldBase*, double, tFieldBase*, char); 
void FDM_Gradient(tField<2> vBar, 
                  double lambda, tField<1> u, char dir) { 
  _Lib_FDM_Gradient1(vBar.Data(), lambda, u.Data(), dir); 
}

EXTERN void _Lib_FDM_Divergence0(tFieldBase*, double, tFieldBase*, char); 
void FDM_Divergence(tField<0> v, 
                    double lambda, tField<1> uBar, char dir) { 
  _Lib_FDM_Divergence0(v.Data(), lambda, uBar.Data(), dir); 
}
EXTERN void _Lib_FDM_Divergence1(tFieldBase*, double, tFieldBase*, char); 
void FDM_Divergence(tField<1> v, 
                    double lambda, tField<2> uBar, char dir) { 
  _Lib_FDM_Divergence1(v.Data(), lambda, uBar.Data(), dir); 
}

EXTERN void _Lib_FDM_Laplacian0( 
    tFieldBase*, double, tFieldBase*); 
void FDM_Laplacian(tField<0> v, 
                   double lambda, tField<0> u) { 
  _Lib_FDM_Laplacian0(v.Data(), lambda, u.Data()); 
}
EXTERN void _Lib_FDM_Laplacian1( 
    tFieldBase*, double, tFieldBase*); 
void FDM_Laplacian(tField<1> v, 
                   double lambda, tField<1> u) { 
  _Lib_FDM_Laplacian1(v.Data(), lambda, u.Data()); 
}
EXTERN void _Lib_FDM_Laplacian2( 
    tFieldBase*, double, tFieldBase*); 
void FDM_Laplacian(tField<2> v, 
                   double lambda, tField<2> u) { 
  _Lib_FDM_Laplacian2(v.Data(), lambda, u.Data()); 
}


// 
// EXPORT LIB FDM CONVECTION 
//
EXTERN void _Lib_FDM_Convection0( 
    tFieldBase*, double, tFieldBase*, tFieldBase*); 
void FDM_Convection(tField<0> v, 
                    double lambda, tField<0> u, tField<1> wBar) { 
  _Lib_FDM_Convection0(v.Data(), lambda, u.Data(), wBar.Data()); 
}
EXTERN void _Lib_FDM_Convection1( 
    tFieldBase*, double, tFieldBase*, tFieldBase*); 
void FDM_Convection(tField<1> v, 
                    double lambda, tField<1> u, tField<1> wBar) { 
  _Lib_FDM_Convection1(v.Data(), lambda, u.Data(), wBar.Data()); 
}


// 
// EXPORT LIB KRYLOV SOLVERS 
//
EXTERN void _Lib_Solve_BiCGStab0(tFieldBase*, tFieldBase*, tMeshOperatorPtr, void*); 
template<typename tMeshOperator> 
void Solve_BiCGStab(tField<0> u, 
                    tField<0> b, tMeshOperator&& meshOperator) { 
  _Lib_Solve_BiCGStab0(u.Data(), b.Data(), 
    [](tFieldBase* out, tFieldBase* in, void* env) { 
      auto& meshOperator = *reinterpret_cast<tMeshOperator*>(env); 
      meshOperator(tField<0>(in, nullptr), 
                   tField<0>(out, nullptr)); 
    }, &meshOperator); 
}
EXTERN void _Lib_Solve_BiCGStab1(tFieldBase*, tFieldBase*, tMeshOperatorPtr, void*); 
template<typename tMeshOperator> 
void Solve_BiCGStab(tField<1> u, 
                    tField<1> b, tMeshOperator&& meshOperator) { 
  _Lib_Solve_BiCGStab1(u.Data(), b.Data(), 
    [](tFieldBase* out, tFieldBase* in, void* env) { 
      auto& meshOperator = *reinterpret_cast<tMeshOperator*>(env); 
      meshOperator(tField<1>(in, nullptr), 
                   tField<1>(out, nullptr)); 
    }, &meshOperator); 
}
EXTERN void _Lib_Solve_BiCGStab2(tFieldBase*, tFieldBase*, tMeshOperatorPtr, void*); 
template<typename tMeshOperator> 
void Solve_BiCGStab(tField<2> u, 
                    tField<2> b, tMeshOperator&& meshOperator) { 
  _Lib_Solve_BiCGStab2(u.Data(), b.Data(), 
    [](tFieldBase* out, tFieldBase* in, void* env) { 
      auto& meshOperator = *reinterpret_cast<tMeshOperator*>(env); 
      meshOperator(tField<2>(in, nullptr), 
                   tField<2>(out, nullptr)); 
    }, &meshOperator); 
}

} // namespace StormRuler 
#undef EXTERN 
// 
// END OF THE AUTO-GENERATED FILE 
// 

