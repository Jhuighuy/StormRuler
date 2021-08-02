!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! Copyright (C) 2021 Oleg Butakov
!! 
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights  to use, 
!! copy, modify, merge, publish, distribute, sublicense, and/or
!! sell copies of the Software, and to permit persons to whom the  
!! Software is furnished to do so, subject to the following 
!! conditions:
!! 
!! The above copyright notice and this permission notice shall be 
!! included in all copies or substantial portions of the Software.
!! 
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!
module StormRuler_Lib

#$use 'StormRuler_Parameters.f90'

use StormRuler_Parameters, only: dp
use StormRuler_Helpers
use StormRuler_IO
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: &
  & Fill,Set,Dot,Add,Sub,FuncProd,SFuncProd
use StormRuler_FDM_Operators, only: &
  & FDM_Gradient_Central,FDM_Divergence_Central, &
  & FDM_Gradient_Forward,FDM_Divergence_Backward, &
  & FDM_Laplacian_Central

use, intrinsic :: iso_c_binding
use, intrinsic :: iso_fortran_env, only: error_unit

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

! Helper function to replace call results with empty sring.
#$let void(x) = ''

! Generated '.inc' file IO.
#$let incFile = open(__FILE__[:-4] + '.inc', mode='w')
#$let writeIncLine(line) = void(incFile.write(line + '\n'))

! Field wrapper.
#$do rank = 0, NUM_RANKS
type :: tField$rank
  real(dp), allocatable :: Data(@:,:)
end type !tField$rank
#$end do

! Math function wrapper.
abstract interface
  pure subroutine tLibMFunc$0(shape,tU,tV,env)
    import :: dp,c_int,c_ptr
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: tU
    real(dp), intent(out) :: tV
    type(c_ptr), intent(in), value :: env
  end subroutine tLibMFunc$0
#$do rank = 1, NUM_RANKS
  pure subroutine tLibMFunc$rank(shape,tU,tV,env)
    import :: dp,c_int,c_ptr
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: tU(*)
    real(dp), intent(out) :: tV(*)
    type(c_ptr), intent(in), value :: env
  end subroutine tLibMFunc$rank
#$end do
end interface

! Spatial math function wrapper.
abstract interface
  pure subroutine tLibSMFunc$0(dim,x,shape,tU,tV,env)
    import :: dp,c_int,c_ptr
    integer(c_int), intent(in), value :: dim
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: x(*),tU
    real(dp), intent(out) :: tV
    type(c_ptr), intent(in), value :: env
  end subroutine tLibSMFunc$0
#$do rank = 1, NUM_RANKS
  pure subroutine tLibSMFunc$rank(dim,x,shape,tU,tV,env)
    import :: dp,c_int,c_ptr
    integer(c_int), intent(in), value :: dim
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: x(*),tU(*)
    real(dp), intent(out) :: tV(*)
    type(c_ptr), intent(in), value :: env
  end subroutine tLibSMFunc$rank
#$end do
end interface

class(tMesh), allocatable :: gMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

${writeIncLine(fr'''// &
// THIS IS THE AUTO-GENERATED FILE, DO NOT EDIT MANUALLY &
// &
#define EXTERN extern "C" &
namespace StormRuler {{ &
''')}$

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

${writeIncLine(fr'''// &
// EXPORT LIB MESH & FIELDS ALLOCATION  &
//''')}$

subroutine Lib_InitializeMesh() &
  & bind(c,name='Lib_InitializeMesh')
  ! <<<<<<<<<<<<<<<<<<<<<<
  ! >>>>>>>>>>>>>>>>>>>>>>

#$if False
  real(8), parameter :: pi = 4*atan(1.0D0)
  Integer, Parameter :: Nx = 100, Ny = 100  
  Real(8), Parameter :: Dx = 2*pi/Nx, Dy = 2*pi/Ny, Dt = Dx*Dx
  allocate(gMesh)
  gMesh%dt = dt
  call gMesh%InitRect(dx,nx,.true.,dy,ny,.true.,20)
#$else
  integer, allocatable :: pixels(:,:)
  integer, allocatable :: colorToBCM(:)
  allocate(gMesh)
  call Load_PPM('test/Domain-10.ppm',pixels)
  colorToBCM = [PixelToInt([255,255,255]),PixelToInt([255,0,0])]
  call gMesh%InitFromImage2D(pixels,0,colorToBCM,3)
  call gMesh%PrintTo_Neato('test/c2c.dot')
  call gMesh%PrintTo_LegacyVTK('test/c2c.vtk')
#$endif
end subroutine Lib_InitializeMesh
${writeIncLine('EXTERN void Lib_InitializeMesh();')}$
${writeIncLine('')}$

${writeIncLine( &
fr'''template<int rank> &
tField<rank> AllocateField();''')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_AllocateField$rank(pField) &
  & bind(c,name='Lib_AllocateField$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(out) :: pField
  ! >>>>>>>>>>>>>>>>>>>>>>
  type(tField$rank), pointer :: field
  allocate(field)
  allocate(field%Data(@{gMesh%Dim}@,gMesh%NumAllCells)); field%Data(@:,:) = 1488.0_dp
  pField = c_loc(field)
end subroutine Lib_AllocateField$rank
${writeIncLine( &
fr'''EXTERN void Lib_AllocateField{rank}(void**); &
template<> &
tField<{rank}> AllocateField<{rank}>() {{ &
  void* pData; &
  Lib_AllocateField{rank}(&pData); &
  return tField<{rank}>(pData); &
}}''')}$
#$end do
${writeIncLine('')}$

${writeIncLine( &
fr'''template<int rank> &
void DeallocateField(tField<rank>);''')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_DeallocateField$rank(pField) &
  & bind(c,name='Lib_DeallocateField$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pField
  ! >>>>>>>>>>>>>>>>>>>>>>
  type(tField$rank), pointer :: field
  call c_f_pointer(cptr=pField, fptr=field)
  deallocate(field%Data)
  deallocate(field)
end subroutine Lib_DeallocateField$rank
${writeIncLine( &
fr'''EXTERN void Lib_DeallocateField{rank}(void*); &
template<> &
void DeallocateField<{rank}>(tField<{rank}> field) {{ &
  Lib_DeallocateField{rank}(field.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$

${writeIncLine('')}$

!! ----------------------------------------------------------------- !!
!! Dereference a C field pointer.
!! ----------------------------------------------------------------- !!
#$do rank = 0, NUM_RANKS
function DEREF$rank(pField) result(u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pField
  real(dp), pointer :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  type(tField$rank), pointer :: field
  call c_f_pointer(cptr=pField, fptr=field)
  u=>field%Data
end function DEREF$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

${writeIncLine(fr'''// &
// EXPORT LIB BLAS &
//''')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Fill$rank(pU,alpha) &
  & bind(c,name='Lib_BLAS_Fill$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: alpha
  type(c_ptr), intent(in), value :: pU
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:)
  u=>DEREF$rank(pU)
  ! ----------------------
  call Fill(gMesh,u,alpha)
end subroutine Lib_BLAS_Fill$rank
${writeIncLine( &
fr'''EXTERN void Lib_BLAS_Fill{rank}(void*, double); &
void BLAS_Fill(tField<{rank}> u, double alpha) {{ &
  Lib_BLAS_Fill{rank}(u.Data(), alpha); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Set$rank(pU,pV) &
  & bind(c,name='Lib_BLAS_Set$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU,pV
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:),v(@:,:)
  u=>DEREF$rank(pU); v=>DEREF$rank(pV) 
  ! ----------------------
  call Set(gMesh,u,v)
end subroutine Lib_BLAS_Set$rank
${writeIncLine( &
fr'''EXTERN void Lib_BLAS_Set{rank}(void*, void*); &
void BLAS_Set(tField<{rank}> u, tField<{rank}> v) {{ &
  Lib_BLAS_Set{rank}(u.Data(), v.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Add$rank(pU,pV,pW,alpha,beta) &
  & bind(c,name='Lib_BLAS_Add$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: alpha,beta
  type(c_ptr), intent(in), value :: pU,pV,pW
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:),v(@:,:),w(@:,:)
  u=>DEREF$rank(pU); v=>DEREF$rank(pV); w=>DEREF$rank(pW)  
  ! ----------------------
  call Add(gMesh,u,v,w,alpha,beta)
end subroutine Lib_BLAS_Add$rank
${writeIncLine( &
fr'''EXTERN void Lib_BLAS_Add{rank}(void*, void*, void*, double, double); &
void BLAS_Add(tField<{rank}> u, tField<{rank}> v, &
              tField<{rank}> w, double alpha, double beta) {{ &
  Lib_BLAS_Add{rank}(u.Data(), v.Data(), w.Data(), alpha, beta); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Sub$rank(pU,pV,pW,alpha,beta) &
  & bind(c,name='Lib_BLAS_Sub$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: alpha,beta
  type(c_ptr), intent(in), value :: pU,pV,pW
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:),v(@:,:),w(@:,:)
  u=>DEREF$rank(pU); v=>DEREF$rank(pV); w=>DEREF$rank(pW)  
  ! ----------------------
  call Sub(gMesh,u,v,w,alpha,beta)
end subroutine Lib_BLAS_Sub$rank
${writeIncLine( &
fr'''EXTERN void Lib_BLAS_Sub{rank}(void*, void*, void*, double, double); &
void BLAS_Sub(tField<{rank}> u, tField<{rank}> v, &
              tField<{rank}> w, double alpha, double beta) {{ &
  Lib_BLAS_Sub{rank}(u.Data(), v.Data(), w.Data(), alpha, beta); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_FuncProd$rank(pV,pU,pF,env) &
  & bind(c,name='Lib_BLAS_FuncProd$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU,pV
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: env
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:),v(@:,:)
  procedure(tLibMFunc$rank), pointer :: f
  u=>DEREF$rank(pU); v=>DEREF$rank(pV)
  call c_f_procpointer(cptr=pF,fptr=f)
  ! ----------------------
  call FuncProd(gMesh,v,u,wF)
contains
  pure function wF(u) result(v)
#$if rank == 0
    real(dp), intent(in) :: u
    real(dp) :: v
#$else
    real(dp), intent(in) :: u(@:)
    real(dp) :: v(@{size(u,dim=$$)}@)
#$endif
    call f(shape(u),u,v,env)
  end function wF 
end subroutine Lib_BLAS_FuncProd$rank
${writeIncLine( &
fr'''EXTERN void Lib_BLAS_FuncProd{rank}(void*, void*, tMFunc, void*); &
template<typename tFunc> &
void BLAS_FuncProd(tField<{rank}> v, &
                   tField<{rank}> u, tFunc&& func) {{ &
  Lib_BLAS_FuncProd{rank}(v.Data(), u.Data(), &
                          [](int* shape, double* in, double* out, &
                             void* env) {{ &
    (*reinterpret_cast<tFunc*>(env))(in, out); &
  }}, &func); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_SFuncProd$rank(pV,pU,pF,env) &
  & bind(c,name='Lib_BLAS_SFuncProd$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU,pV
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: env
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:),v(@:,:)
  procedure(tLibSMFunc$rank), pointer :: f
  u=>DEREF$rank(pU); v=>DEREF$rank(pV)
  call c_f_procpointer(cptr=pF,fptr=f)
  ! ----------------------
  call SFuncProd(gMesh,v,u,wF)
contains
  pure function wF(x,u) result(v)
#$if rank == 0
    real(dp), intent(in) :: x(:),u
    real(dp) :: v
#$else
    real(dp), intent(in) :: x(:),u(@:)
    real(dp) :: v(@{size(u,dim=$$)}@)
#$endif
    call f(size(x),x,shape(u),u,v,env)
  end function wF 
end subroutine Lib_BLAS_SFuncProd$rank
${writeIncLine( &
fr'''EXTERN void Lib_BLAS_SFuncProd{rank}(void*, void*, tSMFunc, void*); &
template<typename tFunc> &
void BLAS_SFuncProd(tField<{rank}> v, &
                    tField<{rank}> u, tFunc&& func) {{ &
  Lib_BLAS_SFuncProd{rank}(v.Data(), u.Data(), &
                           [](int dim, double* x, &
                              int* shape, double* in, double* out, &
                              void* env) {{ &
    (*reinterpret_cast<tFunc*>(env))(x, in, out); &
  }}, &func); &
}}''')}$
#$end do
${writeIncLine('')}$

${writeIncLine('')}$

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

${writeIncLine(fr'''// &
// EXPORT LIB FDM OPERATORS &
//''')}$

#$do rank = 0, NUM_RANKS-1
subroutine Lib_FDM_Gradient$rank(pVBar,lambda,pU,dir) &
  & bind(c,name='Lib_FDM_Gradient$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU,pVBar
  character(c_char), intent(in), value :: dir
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:),vBar(:,@:,:)
  u=>DEREF$rank(pU); vBar=>DEREF${rank+1}$(pVBar)
  ! ----------------------
  select case(dir)
    case('c')
      call FDM_Gradient_Central(gMesh,vBar,lambda,u)
    case('f')
      call FDM_Gradient_Forward(gMesh,vBar,lambda,u,dirAll=1_1)
    case('b')
      call FDM_Gradient_Forward(gMesh,vBar,lambda,u,dirAll=-1_1)
  end select
end subroutine Lib_FDM_Gradient$rank
${writeIncLine( &
fr'''EXTERN void Lib_FDM_Gradient{rank}(void*, double, void*, char); &
void FDM_Gradient(tField<{rank+1}> vBar, &
                  double lambda, tField<{rank}> u, char dir) {{ &
  Lib_FDM_Gradient{rank}(vBar.Data(), lambda, u.Data(), dir); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS-1
subroutine Lib_FDM_Divergence$rank(pV,lambda,pUBar,dir) &
  & bind(c,name='Lib_FDM_Divergence$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pUBar,pV
  character(c_char), intent(in), value :: dir
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: uBar(:,@:,:),v(@:,:)
  uBar=>DEREF${rank+1}$(pUBar); v=>DEREF$rank(pV)
  ! ----------------------
  select case(dir)
    case('c')
      call FDM_Divergence_Central(gMesh,v,lambda,uBar)
    case('f')
      call FDM_Divergence_Backward(gMesh,v,lambda,uBar,dirAll=1_1)
    case('b')
      call FDM_Divergence_Backward(gMesh,v,lambda,uBar,dirAll=-1_1)
  end select
end subroutine Lib_FDM_Divergence$rank
${writeIncLine( &
fr'''EXTERN void Lib_FDM_Divergence{rank}(void*, double, void*, char); &
void FDM_Divergence(tField<{rank}> v, &
                    double lambda, tField<{rank+1}> uBar, char dir) {{ &
  Lib_FDM_Divergence{rank}(v.Data(), lambda, uBar.Data(), dir); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_FDM_Laplacian$rank(pV,lambda,pU) &
  & bind(c,name='Lib_FDM_Laplacian$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU,pV
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:),v(@:,:)
  u=>DEREF$rank(pU); v=>DEREF$rank(pV)
  ! ----------------------
  call FDM_Laplacian_Central(gMesh,v,lambda,u)
end subroutine Lib_FDM_Laplacian$rank
${writeIncLine( &
fr'''EXTERN void Lib_FDM_Laplacian{rank}(void*, double, void*); &
void FDM_Laplacian(tField<{rank}> v, &
                   double lambda, tField<{rank}> u) {{ &
  Lib_FDM_Laplacian{rank}(v.Data(), lambda, u.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$

${writeIncLine('')}$

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

${writeIncLine(fr''' &
}} // namespace StormRuler &
#undef EXTERN &
// &
// END OF THE AUTO-GENERATED FILE &
// &
''')}$

end module StormRuler_Lib

${void(incFile.close())}$
