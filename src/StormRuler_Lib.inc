// 
// THIS IS THE AUTO-GENERATED FILE, DO NOT EDIT MANUALLY 
// 
#define EXTERN extern "C" 
namespace StormRuler { 

// 
// EXPORT LIB MESH & FIELDS ALLOCATION  
//
EXTERN void Lib_InitializeMesh();

template<int rank> 
tField<rank> AllocateField();
EXTERN void Lib_AllocateField0(void**); 
template<> 
tField<0> AllocateField<0>() { 
  void* pData; 
  Lib_AllocateField0(&pData); 
  return tField<0>(pData); 
}
EXTERN void Lib_AllocateField1(void**); 
template<> 
tField<1> AllocateField<1>() { 
  void* pData; 
  Lib_AllocateField1(&pData); 
  return tField<1>(pData); 
}
EXTERN void Lib_AllocateField2(void**); 
template<> 
tField<2> AllocateField<2>() { 
  void* pData; 
  Lib_AllocateField2(&pData); 
  return tField<2>(pData); 
}

template<int rank> 
void DeallocateField(tField<rank>);
EXTERN void Lib_DeallocateField0(void*); 
template<> 
void DeallocateField<0>(tField<0> field) { 
  Lib_DeallocateField0(field.Data()); 
}
EXTERN void Lib_DeallocateField1(void*); 
template<> 
void DeallocateField<1>(tField<1> field) { 
  Lib_DeallocateField1(field.Data()); 
}
EXTERN void Lib_DeallocateField2(void*); 
template<> 
void DeallocateField<2>(tField<2> field) { 
  Lib_DeallocateField2(field.Data()); 
}


// 
// EXPORT LIB BLAS 
//
EXTERN void Lib_BLAS_Fill0(void*, double); 
void BLAS_Fill(tField<0> u, double alpha) { 
  Lib_BLAS_Fill0(u.Data(), alpha); 
}
EXTERN void Lib_BLAS_Fill1(void*, double); 
void BLAS_Fill(tField<1> u, double alpha) { 
  Lib_BLAS_Fill1(u.Data(), alpha); 
}
EXTERN void Lib_BLAS_Fill2(void*, double); 
void BLAS_Fill(tField<2> u, double alpha) { 
  Lib_BLAS_Fill2(u.Data(), alpha); 
}

EXTERN void Lib_BLAS_Set0(void*, void*); 
void BLAS_Set(tField<0> u, tField<0> v) { 
  Lib_BLAS_Set0(u.Data(), v.Data()); 
}
EXTERN void Lib_BLAS_Set1(void*, void*); 
void BLAS_Set(tField<1> u, tField<1> v) { 
  Lib_BLAS_Set1(u.Data(), v.Data()); 
}
EXTERN void Lib_BLAS_Set2(void*, void*); 
void BLAS_Set(tField<2> u, tField<2> v) { 
  Lib_BLAS_Set2(u.Data(), v.Data()); 
}

EXTERN void Lib_BLAS_Add0(void*, void*, void*, double, double); 
void BLAS_Add(tField<0> u, tField<0> v, 
              tField<0> w, double alpha, double beta) { 
  Lib_BLAS_Add0(u.Data(), v.Data(), w.Data(), alpha, beta); 
}
EXTERN void Lib_BLAS_Add1(void*, void*, void*, double, double); 
void BLAS_Add(tField<1> u, tField<1> v, 
              tField<1> w, double alpha, double beta) { 
  Lib_BLAS_Add1(u.Data(), v.Data(), w.Data(), alpha, beta); 
}
EXTERN void Lib_BLAS_Add2(void*, void*, void*, double, double); 
void BLAS_Add(tField<2> u, tField<2> v, 
              tField<2> w, double alpha, double beta) { 
  Lib_BLAS_Add2(u.Data(), v.Data(), w.Data(), alpha, beta); 
}

EXTERN void Lib_BLAS_Sub0(void*, void*, void*, double, double); 
void BLAS_Sub(tField<0> u, tField<0> v, 
              tField<0> w, double alpha, double beta) { 
  Lib_BLAS_Sub0(u.Data(), v.Data(), w.Data(), alpha, beta); 
}
EXTERN void Lib_BLAS_Sub1(void*, void*, void*, double, double); 
void BLAS_Sub(tField<1> u, tField<1> v, 
              tField<1> w, double alpha, double beta) { 
  Lib_BLAS_Sub1(u.Data(), v.Data(), w.Data(), alpha, beta); 
}
EXTERN void Lib_BLAS_Sub2(void*, void*, void*, double, double); 
void BLAS_Sub(tField<2> u, tField<2> v, 
              tField<2> w, double alpha, double beta) { 
  Lib_BLAS_Sub2(u.Data(), v.Data(), w.Data(), alpha, beta); 
}

EXTERN void Lib_BLAS_FuncProd0(void*, void*, tMFunc, void*); 
template<typename tFunc> 
void BLAS_FuncProd(tField<0> v, 
                   tField<0> u, tFunc&& func) { 
  Lib_BLAS_FuncProd0(v.Data(), u.Data(), 
                          [](int* shape, double* in, double* out, 
                             void* env) { 
    (*reinterpret_cast<tFunc*>(env))(in, out); 
  }, &func); 
}
EXTERN void Lib_BLAS_FuncProd1(void*, void*, tMFunc, void*); 
template<typename tFunc> 
void BLAS_FuncProd(tField<1> v, 
                   tField<1> u, tFunc&& func) { 
  Lib_BLAS_FuncProd1(v.Data(), u.Data(), 
                          [](int* shape, double* in, double* out, 
                             void* env) { 
    (*reinterpret_cast<tFunc*>(env))(in, out); 
  }, &func); 
}
EXTERN void Lib_BLAS_FuncProd2(void*, void*, tMFunc, void*); 
template<typename tFunc> 
void BLAS_FuncProd(tField<2> v, 
                   tField<2> u, tFunc&& func) { 
  Lib_BLAS_FuncProd2(v.Data(), u.Data(), 
                          [](int* shape, double* in, double* out, 
                             void* env) { 
    (*reinterpret_cast<tFunc*>(env))(in, out); 
  }, &func); 
}

EXTERN void Lib_BLAS_SFuncProd0(void*, void*, tSMFunc, void*); 
template<typename tFunc> 
void BLAS_SFuncProd(tField<0> v, 
                    tField<0> u, tFunc&& func) { 
  Lib_BLAS_SFuncProd0(v.Data(), u.Data(), 
                           [](int dim, double* x, 
                              int* shape, double* in, double* out, 
                              void* env) { 
    (*reinterpret_cast<tFunc*>(env))(x, in, out); 
  }, &func); 
}
EXTERN void Lib_BLAS_SFuncProd1(void*, void*, tSMFunc, void*); 
template<typename tFunc> 
void BLAS_SFuncProd(tField<1> v, 
                    tField<1> u, tFunc&& func) { 
  Lib_BLAS_SFuncProd1(v.Data(), u.Data(), 
                           [](int dim, double* x, 
                              int* shape, double* in, double* out, 
                              void* env) { 
    (*reinterpret_cast<tFunc*>(env))(x, in, out); 
  }, &func); 
}
EXTERN void Lib_BLAS_SFuncProd2(void*, void*, tSMFunc, void*); 
template<typename tFunc> 
void BLAS_SFuncProd(tField<2> v, 
                    tField<2> u, tFunc&& func) { 
  Lib_BLAS_SFuncProd2(v.Data(), u.Data(), 
                           [](int dim, double* x, 
                              int* shape, double* in, double* out, 
                              void* env) { 
    (*reinterpret_cast<tFunc*>(env))(x, in, out); 
  }, &func); 
}


// 
// EXPORT LIB FDM OPERATORS 
//
EXTERN void Lib_FDM_Gradient0(void*, double, void*, char); 
void FDM_Gradient(tField<1> vBar, 
                  double lambda, tField<0> u, char dir) { 
  Lib_FDM_Gradient0(vBar.Data(), lambda, u.Data(), dir); 
}
EXTERN void Lib_FDM_Gradient1(void*, double, void*, char); 
void FDM_Gradient(tField<2> vBar, 
                  double lambda, tField<1> u, char dir) { 
  Lib_FDM_Gradient1(vBar.Data(), lambda, u.Data(), dir); 
}

EXTERN void Lib_FDM_Divergence0(void*, double, void*, char); 
void FDM_Divergence(tField<0> v, 
                    double lambda, tField<1> uBar, char dir) { 
  Lib_FDM_Divergence0(v.Data(), lambda, uBar.Data(), dir); 
}
EXTERN void Lib_FDM_Divergence1(void*, double, void*, char); 
void FDM_Divergence(tField<1> v, 
                    double lambda, tField<2> uBar, char dir) { 
  Lib_FDM_Divergence1(v.Data(), lambda, uBar.Data(), dir); 
}

EXTERN void Lib_FDM_Laplacian0(void*, double, void*); 
void FDM_Laplacian(tField<0> v, 
                   double lambda, tField<0> u) { 
  Lib_FDM_Laplacian0(v.Data(), lambda, u.Data()); 
}
EXTERN void Lib_FDM_Laplacian1(void*, double, void*); 
void FDM_Laplacian(tField<1> v, 
                   double lambda, tField<1> u) { 
  Lib_FDM_Laplacian1(v.Data(), lambda, u.Data()); 
}
EXTERN void Lib_FDM_Laplacian2(void*, double, void*); 
void FDM_Laplacian(tField<2> v, 
                   double lambda, tField<2> u) { 
  Lib_FDM_Laplacian2(v.Data(), lambda, u.Data()); 
}


 
} // namespace StormRuler 
#undef EXTERN 
// 
// END OF THE AUTO-GENERATED FILE 
// 

