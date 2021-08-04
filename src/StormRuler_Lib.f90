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
  & Fill, Set, Dot, Add, Sub, Mul, &
  & Mul_Inner, Mul_Outer, FuncProd, SFuncProd
use StormRuler_FDM_Operators, only: &
  & FDM_Gradient_Central, FDM_Divergence_Central, &
  & FDM_Gradient_Forward, FDM_Divergence_Backward, &
  & FDM_Laplacian_Central
use StormRuler_FDM_Operators, only: &
  & FDM_Convection_Central ! TODO: should be StormRuler_FDM_Convection
use StormRuler_ConvParams, only: tConvParams
use StormRuler_KrylovSolvers, only: &
  & Solve_CG, Solve_BiCGStab

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
  real(dp), pointer :: Data(@:,:)
end type !tField$rank
#$end do

! Math function wrapper.
abstract interface
  pure subroutine tLibMFunc$0(shape, tU, tV, env)
    import :: dp, c_int, c_ptr
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: tU
    real(dp), intent(out) :: tV
    type(c_ptr), intent(in), value :: env
  end subroutine tLibMFunc$0
#$do rank = 1, NUM_RANKS
  pure subroutine tLibMFunc$rank(shape, tU, tV, env)
    import :: dp, c_int, c_ptr
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: tU(*)
    real(dp), intent(out) :: tV(*)
    type(c_ptr), intent(in), value :: env
  end subroutine tLibMFunc$rank
#$end do
end interface

! Spatial math function wrapper.
abstract interface
  pure subroutine tLibSMFunc$0(dim, x, shape, tU, tV, env)
    import :: dp, c_int, c_ptr
    integer(c_int), intent(in), value :: dim
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: x(*), tU
    real(dp), intent(out) :: tV
    type(c_ptr), intent(in), value :: env
  end subroutine tLibSMFunc$0
#$do rank = 1, NUM_RANKS
  pure subroutine tLibSMFunc$rank(dim, x, shape, tU, tV, env)
    import :: dp, c_int, c_ptr
    integer(c_int), intent(in), value :: dim
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: x(*), tU(*)
    real(dp), intent(out) :: tV(*)
    type(c_ptr), intent(in), value :: env
  end subroutine tLibSMFunc$rank
#$end do
end interface

! Mesh operator function wrapper.
abstract interface
#$do rank = 0, NUM_RANKS
  subroutine tLibMeshOperator$rank(pV, pW, env)
    import :: dp, c_int, c_ptr
    type(c_ptr), intent(in), value :: pV, pW
    type(c_ptr), intent(in), value :: env
  end subroutine tLibMeshOperator$rank
#$end do
end interface

! Global mesh object.
class(tMesh), allocatable :: gMesh

! Global IO list object.
class(IOList), allocatable :: gIOList
integer :: L, M

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
// EXPORT LIB MESH & FIELDS ALLOCATION & FIELDS IO &
//''')}$

subroutine Lib_InitializeMesh() &
  & bind(c, name='Lib_InitializeMesh')
  ! <<<<<<<<<<<<<<<<<<<<<<
  ! >>>>>>>>>>>>>>>>>>>>>>
#$if True
  real(8), parameter :: pi = 4*atan(1.0D0)
  Integer, Parameter :: Nx = 100, Ny = 100
  Real(8), Parameter :: Dx = 2*pi/Nx, Dy = 2*pi/Ny, Dt = Dx*Dx
  allocate(gMesh)
  gMesh%dt = dt
  call gMesh%InitRect(dx, nx, .true., dy, ny, .true., 20)
#$else
  integer, allocatable :: pixels(:,:)
  integer, allocatable :: colorToBCM(:)
  allocate(gMesh)
  call Load_PPM('test/Domain-10.ppm', pixels)
  colorToBCM = [PixelToInt([255, 255, 255]), PixelToInt([255, 0, 0])]
  call gMesh%InitFromImage2D(pixels, 0, colorToBCM, 3)
  call gMesh%PrintTo_Neato('test/c2c.dot')
  call gMesh%PrintTo_LegacyVTK('test/c2c.vtk')
#$endif
  L = 0
end subroutine Lib_InitializeMesh
${writeIncLine('EXTERN void Lib_InitializeMesh();')}$
${writeIncLine('')}$

${writeIncLine( &
fr'''template<int rank> &
tField<rank> AllocateField();''')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_AllocateField$rank(pField) &
  & bind(c, name='_Lib_AllocateField$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(out) :: pField
  ! >>>>>>>>>>>>>>>>>>>>>>
  type(tField$rank), pointer :: field
  allocate(field)
  allocate(field%Data(@{gMesh%Dim}@, gMesh%NumAllCells)); field%Data(@:,:) = 2.0_dp
  pField = c_loc(field)
end subroutine Lib_AllocateField$rank
${writeIncLine( &
fr'''EXTERN void _Lib_AllocateField{rank}(tFieldBase**); &
template<> &
tField<{rank}> AllocateField<{rank}>() {{ &
  tFieldBase* pData; &
  _Lib_AllocateField{rank}(&pData); &
  return tField<{rank}>(pData); &
}}''')}$
#$end do
${writeIncLine('')}$

${writeIncLine( &
fr'''template<int rank> &
void DeallocateField(tFieldBase*);''')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_DeallocateField$rank(pField) &
  & bind(c, name='_Lib_DeallocateField$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pField
  ! >>>>>>>>>>>>>>>>>>>>>>
  type(tField$rank), pointer :: field
  call c_f_pointer(cptr=pField, fptr=field)
  deallocate(field%Data)
  deallocate(field)
end subroutine Lib_DeallocateField$rank
${writeIncLine( &
fr'''EXTERN void _Lib_DeallocateField{rank}(tFieldBase*); &
template<> &
void DeallocateField<{rank}>(tFieldBase* pData) {{ &
  _Lib_DeallocateField{rank}(pData); &
}}''')}$
#$end do
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
  u => field%Data
end function DEREF$rank
#$end do

subroutine Lib_IO_Begin() &
  & bind(c, name='_Lib_IO_Begin')
  ! <<<<<<<<<<<<<<<<<<<<<<
  ! >>>>>>>>>>>>>>>>>>>>>>
  allocate(gIOList)
  M = 1
end subroutine Lib_IO_Begin
${writeIncLine( &
fr'''EXTERN void _Lib_IO_Begin();''')}$
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_IO_Add$rank(pU, pName, nameLen) &
  & bind(c, name='_Lib_IO_Add$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pName
  integer(c_int), intent(in), value :: nameLen
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:)
  character(len=nameLen), pointer :: name
  u => DEREF$rank(pU)
  call c_f_pointer(cptr=pName, fptr=name)
  ! ----------------------
  call gIOList%Add(name, u)
end subroutine Lib_IO_Add$rank
${writeIncLine( &
fr'''EXTERN void _Lib_IO_Add{rank}(tFieldBase*, const char*, int); &
void _Lib_IO_Add(tField<{rank}> u, const std::string& name) {{ &
  _Lib_IO_Add{rank}(u.Data(), name.c_str(), name.size()); &
}}''')}$
#$end do
${writeIncLine('')}$

subroutine Lib_IO_End() &
  & bind(c, name='_Lib_IO_End')
  ! <<<<<<<<<<<<<<<<<<<<<<
  ! >>>>>>>>>>>>>>>>>>>>>>
  call gMesh%PrintTo_LegacyVTK('out/fld-'//I2S(L)//'.vtk', gIOList)
  deallocate(gIOList)
  L = L + 1
end subroutine Lib_IO_End
${writeIncLine( &
fr'''EXTERN void _Lib_IO_End();''')}$
${writeIncLine('')}$

${writeIncLine('')}$

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

${writeIncLine(fr'''// &
// EXPORT LIB BLAS &
//''')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Fill$rank(pU, alpha) &
  & bind(c, name='_Lib_BLAS_Fill$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: alpha
  type(c_ptr), intent(in), value :: pU
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:)
  u => DEREF$rank(pU)
  ! ----------------------
  call Fill(gMesh, u, alpha)
end subroutine Lib_BLAS_Fill$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Fill{rank}(tFieldBase*, double); &
void BLAS_Fill(tField<{rank}> u, double alpha) {{ &
  _Lib_BLAS_Fill{rank}(u.Data(), alpha); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Set$rank(pU, pV) &
  & bind(c, name='_Lib_BLAS_Set$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pV
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:)
  u => DEREF$rank(pU); v => DEREF$rank(pV)
  ! ----------------------
  call Set(gMesh, u, v)
end subroutine Lib_BLAS_Set$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Set{rank}(tFieldBase*, tFieldBase*); &
void BLAS_Set(tField<{rank}> u, tField<{rank}> v) {{ &
  _Lib_BLAS_Set{rank}(u.Data(), v.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Add$rank(pU, pV, pW, alpha, beta) &
  & bind(c, name='_Lib_BLAS_Add$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: alpha, beta
  type(c_ptr), intent(in), value :: pU, pV, pW
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:), w(@:,:)
  u => DEREF$rank(pU); v => DEREF$rank(pV); w => DEREF$rank(pW)
  ! ----------------------
  call Add(gMesh, u, v, w, alpha, beta)
end subroutine Lib_BLAS_Add$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Add{rank}( &
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); &
void BLAS_Add(tField<{rank}> u, tField<{rank}> v, &
              tField<{rank}> w, double alpha, double beta) {{ &
  _Lib_BLAS_Add{rank}(u.Data(), v.Data(), w.Data(), alpha, beta); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Sub$rank(pU, pV, pW, alpha, beta) &
  & bind(c, name='_Lib_BLAS_Sub$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: alpha, beta
  type(c_ptr), intent(in), value :: pU, pV, pW
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:), w(@:,:)
  u => DEREF$rank(pU); v => DEREF$rank(pV); w => DEREF$rank(pW)
  ! ----------------------
  call Sub(gMesh, u, v, w, alpha, beta)
end subroutine Lib_BLAS_Sub$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Sub{rank}( &
  tFieldBase*, tFieldBase*, tFieldBase*, double, double); &
void BLAS_Sub(tField<{rank}> u, tField<{rank}> v, &
              tField<{rank}> w, double alpha, double beta) {{ &
  _Lib_BLAS_Sub{rank}(u.Data(), v.Data(), w.Data(), alpha, beta); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Mul$rank(pU, pV, pW) &
  & bind(c, name='_Lib_BLAS_Mul$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pV, pW
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(:), w(@:,:)
  u => DEREF$rank(pU); v => DEREF$0(pV); w => DEREF$rank(pW)
  ! ----------------------
  call Mul(gMesh, u, v, w)
end subroutine Lib_BLAS_Mul$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Mul{rank}(tFieldBase*, tFieldBase*, tFieldBase*); &
void BLAS_Mul(tField<{rank}> u, &
              tField<0> v, tField<{rank}> w) {{ &
  _Lib_BLAS_Mul{rank}(u.Data(), v.Data(), w.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$

#$if False
#$do rank = 0, NUM_RANKS-1
subroutine Lib_BLAS_Mul_Inner$rank(pU, pVBar, pWBar) &
  & bind(c, name='_Lib_BLAS_Mul_Inner$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pVBar, pWBar
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), vBar(:,:), wBar(:,@:,:)
  u => DEREF$rank(pU); vBar => DEREF$1(pVBar); wBar => DEREF${rank+1}$(pWBar)
  ! ----------------------
  call Mul_Inner(gMesh, u, vBar, wBar)
end subroutine Lib_BLAS_Mul_Inner$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Mul_Inner{rank}(tFieldBase*, tFieldBase*, tFieldBase*); &
void BLAS_Mul_Inner(tField<{rank}> u, &
                    tField<1> vBar, tField<{rank+1}> wBar) {{ &
  _Lib_BLAS_Mul_Inner{rank}(u.Data(), vBar.Data(), wBar.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS-1
subroutine Lib_BLAS_Mul_Outer$rank(pUHat, pVBar, pWBar) &
  & bind(c, name='_Lib_BLAS_Mul_Outer$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pUHat, pVBar, pWBar
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: uHat(:,@:,:), vBar(:,:), wBar(@:,:)
  uHat => DEREF${rank+1}$(pUHat); vBar => DEREF$1(pVBar); wBar => DEREF$rank(pWBar)
  ! ----------------------
  call Mul_Outer(gMesh, uHat, vBar, wBar)
end subroutine Lib_BLAS_Mul_Outer$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Mul_Outer{rank}(tFieldBase*, tFieldBase*, tFieldBase*); &
void BLAS_Mul_Outer(tField<{rank+1}> uHat, &
                    tField<1> vBar, tField<{rank}> wBar) {{ &
  _Lib_BLAS_Mul_Outer{rank}(uHat.Data(), vBar.Data(), wBar.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$
#$end if

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_FuncProd$rank(pV, pU, pF, env) &
  & bind(c, name='_Lib_BLAS_FuncProd$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pV
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: env
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:)
  procedure(tLibMFunc$rank), pointer :: f
  u => DEREF$rank(pU); v => DEREF$rank(pV)
  call c_f_procpointer(cptr=pF, fptr=f)
  ! ----------------------
  call FuncProd(gMesh, v, u, wF)
contains
  pure function wF(u) result(v)
#$if rank == 0
    real(dp), intent(in) :: u
    real(dp) :: v
#$else
    real(dp), intent(in) :: u(@:)
    real(dp) :: v(@{size(u, dim=$$)}@)
#$endif
    call f(shape(u), u, v, env)
  end function wF
end subroutine Lib_BLAS_FuncProd$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_FuncProd{rank}(tFieldBase*, tFieldBase*, tMFuncPtr, void*); &
template<typename tMFunc> &
void BLAS_FuncProd(tField<{rank}> v, &
                   tField<{rank}> u, tMFunc&& func) {{ &
  _Lib_BLAS_FuncProd{rank}(v.Data(), u.Data(), &
    [](int* shape, double* in, double* out, void* env) {{ &
      auto& func = *reinterpret_cast<tMFunc*>(env); ''')}$
#$if rank == 0
${writeIncLine('      *out = func(*in);')}$;
#$else
${writeIncLine('      func(in, out);')}$;
#$endif
${writeIncLine( &
fr'''    }}, &func); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_SFuncProd$rank(pV, pU, pF, env) &
  & bind(c, name='_Lib_BLAS_SFuncProd$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pV
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: env
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:)
  procedure(tLibSMFunc$rank), pointer :: f
  u => DEREF$rank(pU); v => DEREF$rank(pV)
  call c_f_procpointer(cptr=pF, fptr=f)
  ! ----------------------
  call SFuncProd(gMesh, v, u, wF)
contains
  pure function wF(x, u) result(v)
#$if rank == 0
    real(dp), intent(in) :: x(:), u
    real(dp) :: v
#$else
    real(dp), intent(in) :: x(:), u(@:)
    real(dp) :: v(@{size(u, dim=$$)}@)
#$endif
    call f(size(x), x, shape(u), u, v, env)
  end function wF
end subroutine Lib_BLAS_SFuncProd$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_SFuncProd{rank}(tFieldBase*, tFieldBase*, tSMFuncPtr, void*); &
template<typename tSMFunc> &
void BLAS_SFuncProd(tField<{rank}> v, &
                    tField<{rank}> u, tSMFunc&& func) {{ &
  _Lib_BLAS_SFuncProd{rank}(v.Data(), u.Data(), &
    [](int dim, double* x, int* shape, double* in, double* out, void* env) {{ &
      auto& func = *reinterpret_cast<tSMFunc*>(env); & 
      func(x, in, out); &
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
subroutine Lib_FDM_Gradient$rank(pVBar, lambda, pU, dir) &
  & bind(c, name='_Lib_FDM_Gradient$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pVBar
  character(c_char), intent(in), value :: dir
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), vBar(:,@:,:)
  u => DEREF$rank(pU); vBar => DEREF${rank+1}$(pVBar)
  ! ----------------------
  select case(dir)
    case('c')
      call FDM_Gradient_Central(gMesh, vBar, lambda, u)
    case('f')
      call FDM_Gradient_Forward(gMesh, vBar, lambda, u, dirAll=1_1)
    case('b')
      call FDM_Gradient_Forward(gMesh, vBar, lambda, u, dirAll=-1_1)
  end select
end subroutine Lib_FDM_Gradient$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FDM_Gradient{rank}( &
    tFieldBase*, double, tFieldBase*, char); &
void FDM_Gradient(tField<{rank+1}> vBar, &
                  double lambda, tField<{rank}> u, char dir) {{ &
  _Lib_FDM_Gradient{rank}(vBar.Data(), lambda, u.Data(), dir); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS-1
subroutine Lib_FDM_Divergence$rank(pV, lambda, pUBar, dir) &
  & bind(c, name='_Lib_FDM_Divergence$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pUBar, pV
  character(c_char), intent(in), value :: dir
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: uBar(:,@:,:), v(@:,:)
  uBar => DEREF${rank+1}$(pUBar); v => DEREF$rank(pV)
  ! ----------------------
  select case(dir)
    case('c')
      call FDM_Divergence_Central(gMesh, v, lambda, uBar)
    case('f')
      call FDM_Divergence_Backward(gMesh, v, lambda, uBar, dirAll=1_1)
    case('b')
      call FDM_Divergence_Backward(gMesh, v, lambda, uBar, dirAll=-1_1)
  end select
end subroutine Lib_FDM_Divergence$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FDM_Divergence{rank}(tFieldBase*, double, tFieldBase*, char); &
void FDM_Divergence(tField<{rank}> v, &
                    double lambda, tField<{rank+1}> uBar, char dir) {{ &
  _Lib_FDM_Divergence{rank}(v.Data(), lambda, uBar.Data(), dir); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_FDM_Laplacian$rank(pV, lambda, pU) &
  & bind(c, name='_Lib_FDM_Laplacian$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:)
  u => DEREF$rank(pU); v => DEREF$rank(pV)
  ! ----------------------
  call FDM_Laplacian_Central(gMesh, v, lambda, u)
end subroutine Lib_FDM_Laplacian$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FDM_Laplacian{rank}( &
    tFieldBase*, double, tFieldBase*); &
void FDM_Laplacian(tField<{rank}> v, &
                   double lambda, tField<{rank}> u) {{ &
  _Lib_FDM_Laplacian{rank}(v.Data(), lambda, u.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$

${writeIncLine('')}$

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

${writeIncLine(fr'''// &
// EXPORT LIB FDM CONVECTION &
//''')}$

#$do rank = 0, NUM_RANKS-1
subroutine Lib_FDM_Convection$rank(pV, lambda, pU, pWBar) &
  & bind(c, name='_Lib_FDM_Convection$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV, pWBar
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:), wBar(:,:)
  u => DEREF$rank(pU); v => DEREF$rank(pV); wBar => DEREF$1(pWBar)
  ! ----------------------
  call FDM_Convection_Central(gMesh, v, lambda, u, wBar)
end subroutine Lib_FDM_Convection$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FDM_Convection{rank}( &
    tFieldBase*, double, tFieldBase*, tFieldBase*); &
void FDM_Convection(tField<{rank}> v, &
                    double lambda, tField<{rank}> u, tField<1> wBar) {{ &
  _Lib_FDM_Convection{rank}(v.Data(), lambda, u.Data(), wBar.Data()); &
}}''')}$
#$end do
${writeIncLine('')}$

${writeIncLine('')}$

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

${writeIncLine(fr'''// &
// EXPORT LIB KRYLOV SOLVERS &
//''')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_Solve_BiCGStab$rank(pU, pB, pA, env) &
  & bind(c, name='_Lib_Solve_BiCGStab$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pB
  type(c_funptr), intent(in), value :: pA
  type(c_ptr), intent(in), value :: env
  ! >>>>>>>>>>>>>>>>>>>>>>
  class(tConvParams), allocatable :: Params
  real(dp), pointer :: u(@:,:), b(@:,:)
  procedure(tLibMeshOperator$rank), pointer :: A
  u => DEREF$rank(pU); b => DEREF$rank(pB)
  call c_f_procpointer(cptr=pA, fptr=A)
  ! ----------------------
  allocate(Params)
  call Params%Init(gMesh%Dl(1)*gMesh%Dl(2)*1.0D-8, &
    &              gMesh%Dl(1)*gMesh%Dl(1)*1.0D-8, 100000)
  call Solve_CG(gMesh, u, b, wA, Params, Params)
contains
  subroutine wA(mesh, v, w, opParams)
    class(tMesh), intent(in) :: mesh
    real(dp), intent(in), pointer :: v(@:,:), w(@:,:)
    class(*), intent(in) :: opParams
    type(tField$rank), target :: fV, fW
    type(c_ptr) :: pV, pW
    fV%Data => v; fW%Data => w
    pV = c_loc(fV); pW = c_loc(fW)
    call A(pV, pW, env)
  end subroutine wA
end subroutine Lib_Solve_BiCGStab$rank
${writeIncLine( &
fr'''EXTERN void _Lib_Solve_BiCGStab{rank}(tFieldBase*, tFieldBase*, tMeshOperatorPtr, void*); &
template<typename tMeshOperator> &
void Solve_BiCGStab(tField<{rank}> u, &
                    tField<{rank}> b, tMeshOperator&& meshOperator) {{ &
  _Lib_Solve_BiCGStab{rank}(u.Data(), b.Data(), &
    [](tFieldBase* out, tFieldBase* in, void* env) {{ &
      auto& meshOperator = *reinterpret_cast<tMeshOperator*>(env); &
      meshOperator(tField<{rank}>(in, nullptr), &
                   tField<{rank}>(out, nullptr)); & 
    }}, &meshOperator); &
}}''')}$
#$end do

${writeIncLine('')}$

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

${writeIncLine(fr'''}} // namespace StormRuler &
#undef EXTERN &
// &
// END OF THE AUTO-GENERATED FILE &
// &
''')}$

end module StormRuler_Lib

${void(incFile.close())}$
