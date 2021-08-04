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
tField<rank> AllocField();
EXTERN void _Lib_AllocField0(tFieldBase**); 
template<> 
tField<0> AllocField<0>() { 
  tFieldBase* pBase; 
  _Lib_AllocField0(&pBase); 
  return tField<0>(pBase); 
}
EXTERN void _Lib_AllocField1(tFieldBase**); 
template<> 
tField<1> AllocField<1>() { 
  tFieldBase* pBase; 
  _Lib_AllocField1(&pBase); 
  return tField<1>(pBase); 
}
EXTERN void _Lib_AllocField2(tFieldBase**); 
template<> 
tField<2> AllocField<2>() { 
  tFieldBase* pBase; 
  _Lib_AllocField2(&pBase); 
  return tField<2>(pBase); 
}

template<int rank> 
void FreeField(tFieldBase*);
EXTERN void _Lib_FreeField0(tFieldBase*); 
template<> 
void FreeField<0>(tFieldBase* pData) { 
  if (pData != nullptr) _Lib_FreeField0(pData); 
}
EXTERN void _Lib_FreeField1(tFieldBase*); 
template<> 
void FreeField<1>(tFieldBase* pData) { 
  if (pData != nullptr) _Lib_FreeField1(pData); 
}
EXTERN void _Lib_FreeField2(tFieldBase*); 
template<> 
void FreeField<2>(tFieldBase* pData) { 
  if (pData != nullptr) _Lib_FreeField2(pData); 
}

EXTERN void _Lib_GetFieldData0(tFieldBase*, double**, int*); 
void GetFieldData(tField<0> field, 
                  double** ppFieldData, int* pFieldDataSize) { 
  _Lib_GetFieldData0(field.Base(), ppFieldData, pFieldDataSize); 
}
EXTERN void _Lib_GetFieldData1(tFieldBase*, double**, int*); 
void GetFieldData(tField<1> field, 
                  double** ppFieldData, int* pFieldDataSize) { 
  _Lib_GetFieldData1(field.Base(), ppFieldData, pFieldDataSize); 
}
EXTERN void _Lib_GetFieldData2(tFieldBase*, double**, int*); 
void GetFieldData(tField<2> field, 
                  double** ppFieldData, int* pFieldDataSize) { 
  _Lib_GetFieldData2(field.Base(), ppFieldData, pFieldDataSize); 
}

EXTERN void _Lib_IO_Begin(...);

EXTERN void _Lib_IO_Add0(tFieldBase*, const char*, int); 
void _Lib_IO_Add(tField<0> field, const std::string& name) { 
  _Lib_IO_Add0(field.Base(), name.c_str(), name.size()); 
}
EXTERN void _Lib_IO_Add1(tFieldBase*, const char*, int); 
void _Lib_IO_Add(tField<1> field, const std::string& name) { 
  _Lib_IO_Add1(field.Base(), name.c_str(), name.size()); 
}
EXTERN void _Lib_IO_Add2(tFieldBase*, const char*, int); 
void _Lib_IO_Add(tField<2> field, const std::string& name) { 
  _Lib_IO_Add2(field.Base(), name.c_str(), name.size()); 
}

EXTERN void _Lib_IO_End();


// 
// EXPORT LIB BLAS 
//
EXTERN void _Lib_BLAS_Fill0(tFieldBase*, double); 
void BLAS_Fill(tField<0> u, double alpha) { 
  _Lib_BLAS_Fill0(u.Base(), alpha); 
}
EXTERN void _Lib_BLAS_Fill1(tFieldBase*, double); 
void BLAS_Fill(tField<1> u, double alpha) { 
  _Lib_BLAS_Fill1(u.Base(), alpha); 
}
EXTERN void _Lib_BLAS_Fill2(tFieldBase*, double); 
void BLAS_Fill(tField<2> u, double alpha) { 
  _Lib_BLAS_Fill2(u.Base(), alpha); 
}

EXTERN void _Lib_BLAS_Set0(tFieldBase*, tFieldBase*); 
void BLAS_Set(tField<0> u, tField<0> v) { 
  _Lib_BLAS_Set0(u.Base(), v.Base()); 
}
EXTERN void _Lib_BLAS_Set1(tFieldBase*, tFieldBase*); 
void BLAS_Set(tField<1> u, tField<1> v) { 
  _Lib_BLAS_Set1(u.Base(), v.Base()); 
}
EXTERN void _Lib_BLAS_Set2(tFieldBase*, tFieldBase*); 
void BLAS_Set(tField<2> u, tField<2> v) { 
  _Lib_BLAS_Set2(u.Base(), v.Base()); 
}

EXTERN void _Lib_BLAS_Add0( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Add(tField<0> u, tField<0> v, 
              tField<0> w, double alpha, double beta) { 
  _Lib_BLAS_Add0(u.Base(), v.Base(), w.Base(), alpha, beta); 
}
EXTERN void _Lib_BLAS_Add1( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Add(tField<1> u, tField<1> v, 
              tField<1> w, double alpha, double beta) { 
  _Lib_BLAS_Add1(u.Base(), v.Base(), w.Base(), alpha, beta); 
}
EXTERN void _Lib_BLAS_Add2( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Add(tField<2> u, tField<2> v, 
              tField<2> w, double alpha, double beta) { 
  _Lib_BLAS_Add2(u.Base(), v.Base(), w.Base(), alpha, beta); 
}

EXTERN void _Lib_BLAS_Sub0( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Sub(tField<0> u, tField<0> v, 
              tField<0> w, double alpha, double beta) { 
  _Lib_BLAS_Sub0(u.Base(), v.Base(), w.Base(), alpha, beta); 
}
EXTERN void _Lib_BLAS_Sub1( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Sub(tField<1> u, tField<1> v, 
              tField<1> w, double alpha, double beta) { 
  _Lib_BLAS_Sub1(u.Base(), v.Base(), w.Base(), alpha, beta); 
}
EXTERN void _Lib_BLAS_Sub2( 
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); 
void BLAS_Sub(tField<2> u, tField<2> v, 
              tField<2> w, double alpha, double beta) { 
  _Lib_BLAS_Sub2(u.Base(), v.Base(), w.Base(), alpha, beta); 
}

EXTERN void _Lib_BLAS_Mul0( 
    tFieldBase*, tFieldBase*, tFieldBase*, int); 
void BLAS_Mul(tField<0> u, 
              tField<0> v, tField<0> w, int power = 1) { 
  _Lib_BLAS_Mul0(u.Base(), v.Base(), w.Base(), power); 
}
EXTERN void _Lib_BLAS_Mul1( 
    tFieldBase*, tFieldBase*, tFieldBase*, int); 
void BLAS_Mul(tField<1> u, 
              tField<0> v, tField<1> w, int power = 1) { 
  _Lib_BLAS_Mul1(u.Base(), v.Base(), w.Base(), power); 
}
EXTERN void _Lib_BLAS_Mul2( 
    tFieldBase*, tFieldBase*, tFieldBase*, int); 
void BLAS_Mul(tField<2> u, 
              tField<0> v, tField<2> w, int power = 1) { 
  _Lib_BLAS_Mul2(u.Base(), v.Base(), w.Base(), power); 
}

EXTERN void _Lib_BLAS_FuncProd0( 
    tFieldBase*, tFieldBase*, tMFuncPtr, void*); 
template<typename tMFunc> 
void BLAS_FuncProd(tField<0> v, 
                   tField<0> u, tMFunc&& func) { 
  _Lib_BLAS_FuncProd0(v.Base(), u.Base(), 
    [](int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tMFunc*>(env); 
      *out = func(*in);
    }, &func); 
}
EXTERN void _Lib_BLAS_FuncProd1( 
    tFieldBase*, tFieldBase*, tMFuncPtr, void*); 
template<typename tMFunc> 
void BLAS_FuncProd(tField<1> v, 
                   tField<1> u, tMFunc&& func) { 
  _Lib_BLAS_FuncProd1(v.Base(), u.Base(), 
    [](int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tMFunc*>(env); 
      func(in, out);
    }, &func); 
}
EXTERN void _Lib_BLAS_FuncProd2( 
    tFieldBase*, tFieldBase*, tMFuncPtr, void*); 
template<typename tMFunc> 
void BLAS_FuncProd(tField<2> v, 
                   tField<2> u, tMFunc&& func) { 
  _Lib_BLAS_FuncProd2(v.Base(), u.Base(), 
    [](int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tMFunc*>(env); 
      func(in, out);
    }, &func); 
}

EXTERN void _Lib_BLAS_SFuncProd0( 
    tFieldBase*, tFieldBase*, tSMFuncPtr, void*); 
template<typename tSMFunc> 
void BLAS_SFuncProd(tField<0> v, 
                    tField<0> u, tSMFunc&& func) { 
  _Lib_BLAS_SFuncProd0(v.Base(), u.Base(), 
    [](int dim, double* x, int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tSMFunc*>(env); 
      func(x, in, out); 
    }, &func); 
}
EXTERN void _Lib_BLAS_SFuncProd1( 
    tFieldBase*, tFieldBase*, tSMFuncPtr, void*); 
template<typename tSMFunc> 
void BLAS_SFuncProd(tField<1> v, 
                    tField<1> u, tSMFunc&& func) { 
  _Lib_BLAS_SFuncProd1(v.Base(), u.Base(), 
    [](int dim, double* x, int* shape, double* in, double* out, void* env) { 
      auto& func = *reinterpret_cast<tSMFunc*>(env); 
      func(x, in, out); 
    }, &func); 
}
EXTERN void _Lib_BLAS_SFuncProd2( 
    tFieldBase*, tFieldBase*, tSMFuncPtr, void*); 
template<typename tSMFunc> 
void BLAS_SFuncProd(tField<2> v, 
                    tField<2> u, tSMFunc&& func) { 
  _Lib_BLAS_SFuncProd2(v.Base(), u.Base(), 
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
void FDM_Gradient(tField<1> v_bar, 
                  double lambda, tField<0> u, char dir) { 
  _Lib_FDM_Gradient0(v_bar.Base(), lambda, u.Base(), dir); 
}
EXTERN void _Lib_FDM_Gradient1( 
    tFieldBase*, double, tFieldBase*, char); 
void FDM_Gradient(tField<2> v_bar, 
                  double lambda, tField<1> u, char dir) { 
  _Lib_FDM_Gradient1(v_bar.Base(), lambda, u.Base(), dir); 
}

EXTERN void _Lib_FDM_Divergence0( 
    tFieldBase*, double, tFieldBase*, char); 
void FDM_Divergence(tField<0> v, 
                    double lambda, tField<1> u_bar, char dir) { 
  _Lib_FDM_Divergence0(v.Base(), lambda, u_bar.Base(), dir); 
}
EXTERN void _Lib_FDM_Divergence1( 
    tFieldBase*, double, tFieldBase*, char); 
void FDM_Divergence(tField<1> v, 
                    double lambda, tField<2> u_bar, char dir) { 
  _Lib_FDM_Divergence1(v.Base(), lambda, u_bar.Base(), dir); 
}

EXTERN void _Lib_FDM_Laplacian0( 
    tFieldBase*, double, tFieldBase*); 
void FDM_Laplacian(tField<0> v, 
                   double lambda, tField<0> u) { 
  _Lib_FDM_Laplacian0(v.Base(), lambda, u.Base()); 
}
EXTERN void _Lib_FDM_Laplacian1( 
    tFieldBase*, double, tFieldBase*); 
void FDM_Laplacian(tField<1> v, 
                   double lambda, tField<1> u) { 
  _Lib_FDM_Laplacian1(v.Base(), lambda, u.Base()); 
}
EXTERN void _Lib_FDM_Laplacian2( 
    tFieldBase*, double, tFieldBase*); 
void FDM_Laplacian(tField<2> v, 
                   double lambda, tField<2> u) { 
  _Lib_FDM_Laplacian2(v.Base(), lambda, u.Base()); 
}


// 
// EXPORT LIB FDM CONVECTION 
//
EXTERN void _Lib_FDM_Convection0( 
    tFieldBase*, double, tFieldBase*, tFieldBase*); 
void FDM_Convection(tField<0> v, 
                    double lambda, tField<0> u, tField<1> w_bar) { 
  _Lib_FDM_Convection0(v.Base(), lambda, u.Base(), w_bar.Base()); 
}
EXTERN void _Lib_FDM_Convection1( 
    tFieldBase*, double, tFieldBase*, tFieldBase*); 
void FDM_Convection(tField<1> v, 
                    double lambda, tField<1> u, tField<1> w_bar) { 
  _Lib_FDM_Convection1(v.Base(), lambda, u.Base(), w_bar.Base()); 
}


// 
// EXPORT LIB KRYLOV SOLVERS 
//
EXTERN void _Lib_Solve_BiCGStab0( 
    tFieldBase*, tFieldBase*, tMeshOperatorPtr, void*); 
template<typename tMeshOperator> 
void Solve_BiCGStab(tField<0> u, 
                    tField<0> b, tMeshOperator&& meshOperator) { 
  _Lib_Solve_BiCGStab0(u.Base(), b.Base(), 
    [](tFieldBase* out, tFieldBase* in, void* env) { 
      auto& meshOperator = *reinterpret_cast<tMeshOperator*>(env); 
      meshOperator(tField<0>(in, nullptr), 
                   tField<0>(out, nullptr)); 
    }, &meshOperator); 
}
EXTERN void _Lib_Solve_BiCGStab1( 
    tFieldBase*, tFieldBase*, tMeshOperatorPtr, void*); 
template<typename tMeshOperator> 
void Solve_BiCGStab(tField<1> u, 
                    tField<1> b, tMeshOperator&& meshOperator) { 
  _Lib_Solve_BiCGStab1(u.Base(), b.Base(), 
    [](tFieldBase* out, tFieldBase* in, void* env) { 
      auto& meshOperator = *reinterpret_cast<tMeshOperator*>(env); 
      meshOperator(tField<1>(in, nullptr), 
                   tField<1>(out, nullptr)); 
    }, &meshOperator); 
}
EXTERN void _Lib_Solve_BiCGStab2( 
    tFieldBase*, tFieldBase*, tMeshOperatorPtr, void*); 
template<typename tMeshOperator> 
void Solve_BiCGStab(tField<2> u, 
                    tField<2> b, tMeshOperator&& meshOperator) { 
  _Lib_Solve_BiCGStab2(u.Base(), b.Base(), 
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

