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
  & FDM_Laplacian_Central, FDM_DivWGrad_Central
use StormRuler_FDM_Operators, only: &
  & FDM_Convection_Central ! TODO: should be StormRuler_FDM_Convection
use StormRuler_ConvParams, only: tConvParams
use StormRuler_KrylovSolvers, only: &
  & Solve_CG, Solve_BiCGStab
#$if HAS_MKL
use StormRuler_KrylovSolvers_MKL, only: &
  & Solve_CG_MKL, Solve_FGMRES_MKL
#$end if

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: &
  & c_char, c_int, c_size_t, c_ptr, c_funptr, &
  & c_loc, c_f_pointer, c_f_procpointer   

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
type :: tLibFieldBase$rank
  real(dp), pointer :: Values(@:,:)
end type !tLibFieldBase$rank
#$end do

! Math function wrapper.
abstract interface
  pure subroutine tLibMapFunc$0(shape, tU, tV, env) bind(c)
    import :: dp, c_int, c_ptr
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: tU
    real(dp), intent(out) :: tV
    type(c_ptr), intent(in), value :: env
  end subroutine tLibMapFunc$0
#$do rank = 1, NUM_RANKS
  pure subroutine tLibMapFunc$rank(shape, tU, tV, env) bind(c)
    import :: dp, c_int, c_ptr
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: tU(*)
    real(dp), intent(out) :: tV(*)
    type(c_ptr), intent(in), value :: env
  end subroutine tLibMapFunc$rank
#$end do
end interface

! Spatial math function wrapper.
abstract interface
  pure subroutine tLibSMapFunc$0(dim, x, shape, tU, tV, env) bind(c)
    import :: dp, c_int, c_ptr
    integer(c_int), intent(in), value :: dim
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: x(*), tU
    real(dp), intent(out) :: tV
    type(c_ptr), intent(in), value :: env
  end subroutine tLibSMapFunc$0
#$do rank = 1, NUM_RANKS
  pure subroutine tLibSMapFunc$rank(dim, x, shape, tU, tV, env) bind(c)
    import :: dp, c_int, c_ptr
    integer(c_int), intent(in), value :: dim
    integer(c_int), intent(in) :: shape(*)
    real(dp), intent(in) :: x(*), tU(*)
    real(dp), intent(out) :: tV(*)
    type(c_ptr), intent(in), value :: env
  end subroutine tLibSMapFunc$rank
#$end do
end interface

! Mesh operator function wrapper.
abstract interface
#$do rank = 0, NUM_RANKS
  subroutine tLibMeshOperator$rank(pV, pW, env) bind(c)
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
integer(ip) :: L

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
  integer(ip), Parameter :: Nx = 100, Ny = 100
  Real(8), Parameter :: Dx = 2*pi/Nx, Dy = 2*pi/Ny, Dt = Dx*Dx
  allocate(gMesh)
  gMesh%dt = dt
  call gMesh%InitRect(dx, nx, .true., dy, ny, .true., 20)
#$else
  integer(ip), allocatable :: pixels(:,:)
  integer(ip), allocatable :: colorToBCM(:)
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
tField<rank> AllocField();''')}$
#$do rank = 0, NUM_RANKS
subroutine Lib_AllocField$rank(ppField) &
  & bind(c, name='_Lib_AllocField$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(out) :: ppField
  ! >>>>>>>>>>>>>>>>>>>>>>
  type(tLibFieldBase$rank), pointer :: pFieldBase
  allocate(pFieldBase)
  allocate(pFieldBase%Values(@{gMesh%Dim}@, gMesh%NumAllCells))
  ppField = c_loc(pFieldBase)
end subroutine Lib_AllocField$rank
${writeIncLine( &
fr'''EXTERN void _Lib_AllocField{rank}(tFieldBase**); &
template<> &
tField<{rank}> AllocField<{rank}>() {{ &
  tFieldBase* pBase; &
  _Lib_AllocField{rank}(&pBase); &
  return tField<{rank}>(pBase); &
}}''')}$
#$end do
${writeIncLine('')}$

${writeIncLine( &
fr'''template<int rank> &
void FreeField(tFieldBase*);''')}$
#$do rank = 0, NUM_RANKS
subroutine Lib_FreeField$rank(pField) &
  & bind(c, name='_Lib_FreeField$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pField
  ! >>>>>>>>>>>>>>>>>>>>>>
  type(tLibFieldBase$rank), pointer :: pFieldBase
  call c_f_pointer(cptr=pField, fptr=pFieldBase)
  deallocate(pFieldBase%Values)
  deallocate(pFieldBase)
end subroutine Lib_FreeField$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FreeField{rank}(tFieldBase*); &
template<> &
void FreeField<{rank}>(tFieldBase* pData) {{ &
  if (pData != nullptr) _Lib_FreeField{rank}(pData); &
}}''')}$
#$end do
${writeIncLine('')}$

!! ----------------------------------------------------------------- !!
!! Dereference a C field pointer.
!! ----------------------------------------------------------------- !!
#$do rank = 0, NUM_RANKS
function DEREF$rank(pField) result(field)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pField
  real(dp), pointer :: field(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  type(tLibFieldBase$rank), pointer :: pFieldBase
  call c_f_pointer(cptr=pField, fptr=pFieldBase)
  field => pFieldBase%Values
end function DEREF$rank
#$end do

#$do rank = 0, NUM_RANKS
subroutine Lib_GetFieldData$rank(pField, ppFieldData, pFieldDataSize) &
  & bind(c, name='_Lib_GetFieldData$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pField
  type(c_ptr), intent(out) :: ppFieldData
  integer(c_size_t), intent(out) :: pFieldDataSize
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: field(@:,:)
  field => DEREF$rank(pField)
  ppFieldData = c_loc(field); pFieldDataSize = size(field)
end subroutine Lib_GetFieldData$rank
${writeIncLine( &
fr'''EXTERN void _Lib_GetFieldData{rank}(tFieldBase*, double**, int*); &
void GetFieldData(tField<{rank}> field, &
                  double** ppFieldData, int* pFieldDataSize) {{ &
  _Lib_GetFieldData{rank}(field.Base(), ppFieldData, pFieldDataSize); &
}}''')}$
#$end do
${writeIncLine('')}$

!!
!! TODO: add IO parameters:
!! 'Begin' should allocate IO object of the specific type,
!! 'End' should envoke the 'Write' action and deallocate it.
!! 'Add' looks fine (for now).
!!

subroutine Lib_IO_Begin() &
  & bind(c, name='_Lib_IO_Begin')
  ! <<<<<<<<<<<<<<<<<<<<<<
  ! >>>>>>>>>>>>>>>>>>>>>>
  allocate(gIOList)
end subroutine Lib_IO_Begin
${writeIncLine( &
fr'''EXTERN void _Lib_IO_Begin(...);''')}$
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_IO_Add$rank(pField, pName, nameLen) &
  & bind(c, name='_Lib_IO_Add$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pField, pName
  integer(c_int), intent(in), value :: nameLen
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: field(@:,:)
  character(len=nameLen), pointer :: name
  field => DEREF$rank(pField)
  call c_f_pointer(cptr=pName, fptr=name)
  ! ----------------------
  call gIOList%Add(name, field)
end subroutine Lib_IO_Add$rank
${writeIncLine( &
fr'''EXTERN void _Lib_IO_Add{rank}(tFieldBase*, const char*, int); &
void _Lib_IO_Add(tField<{rank}> field, const std::string& name) {{ &
  _Lib_IO_Add{rank}(field.Base(), name.c_str(), name.size()); &
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
  _Lib_BLAS_Fill{rank}(u.Base(), alpha); &
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
  _Lib_BLAS_Set{rank}(u.Base(), v.Base()); &
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
  _Lib_BLAS_Add{rank}(u.Base(), v.Base(), w.Base(), alpha, beta); &
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
  _Lib_BLAS_Sub{rank}(u.Base(), v.Base(), w.Base(), alpha, beta); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_BLAS_Mul$rank(pU, pV, pW, power) &
  & bind(c, name='_Lib_BLAS_Mul$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pV, pW
  integer(c_int), intent(in), value :: power
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(:), w(@:,:)
  u => DEREF$rank(pU); v => DEREF$0(pV); w => DEREF$rank(pW)
  ! ----------------------
  call Mul(gMesh, u, v, w, int(power))
end subroutine Lib_BLAS_Mul$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Mul{rank}( &
    tFieldBase*, tFieldBase*, tFieldBase*, int); &
void BLAS_Mul(tField<{rank}> u, &
              tField<0> v, tField<{rank}> w, int power = 1) {{ &
  _Lib_BLAS_Mul{rank}(u.Base(), v.Base(), w.Base(), power); &
}}''')}$
#$end do
${writeIncLine('')}$

#$if False
#$do rank = 0, NUM_RANKS-1
subroutine Lib_BLAS_Mul_Inner$rank(pU, pV_bar, pW_bar) &
  & bind(c, name='_Lib_BLAS_Mul_Inner$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pU, pV_bar, pW_bar
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v_bar(:,:), w_bar(:,@:,:)
  u => DEREF$rank(pU); v_bar => DEREF$1(pV_bar); w_bar => DEREF${rank+1}$(pW_bar)
  ! ----------------------
  call Mul_Inner(gMesh, u, v_bar, w_bar)
end subroutine Lib_BLAS_Mul_Inner$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Mul_Inner{rank}(tFieldBase*, tFieldBase*, tFieldBase*); &
void BLAS_Mul_Inner(tField<{rank}> u, &
                    tField<1> v_bar, tField<{rank+1}> w_bar) {{ &
  _Lib_BLAS_Mul_Inner{rank}(u.Base(), v_bar.Base(), w_bar.Base()); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS-1
subroutine Lib_BLAS_Mul_Outer$rank(pUHat, pV_bar, pW_bar) &
  & bind(c, name='_Lib_BLAS_Mul_Outer$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pUHat, pV_bar, pW_bar
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: uHat(:,@:,:), v_bar(:,:), w_bar(@:,:)
  uHat => DEREF${rank+1}$(pUHat); v_bar => DEREF$1(pV_bar); w_bar => DEREF$rank(pW_bar)
  ! ----------------------
  call Mul_Outer(gMesh, uHat, v_bar, w_bar)
end subroutine Lib_BLAS_Mul_Outer$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_Mul_Outer{rank}(tFieldBase*, tFieldBase*, tFieldBase*); &
void BLAS_Mul_Outer(tField<{rank+1}> uHat, &
                    tField<1> v_bar, tField<{rank}> w_bar) {{ &
  _Lib_BLAS_Mul_Outer{rank}(uHat.Base(), v_bar.Base(), w_bar.Base()); &
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
  procedure(tLibMapFunc$rank), pointer :: f
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
    call f(int(shape(u), kind=c_int), u, v, env)
  end function wF
end subroutine Lib_BLAS_FuncProd$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_FuncProd{rank}( &
    tFieldBase*, tFieldBase*, tMFuncPtr, void*); &
template<typename tMFunc> &
void BLAS_FuncProd(tField<{rank}> v, &
                   tField<{rank}> u, tMFunc&& func) {{ &
  _Lib_BLAS_FuncProd{rank}(v.Base(), u.Base(), &
    [](int* shape, double* in, double* out, void* env) {{ &
      auto& func = *reinterpret_cast<tMFunc*>(env); ''')}$
#$if rank == 0
${writeIncLine('      *out = func(*in);')}$
#$else
${writeIncLine('      func(in, out);')}$
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
  procedure(tLibSMapFunc$rank), pointer :: f
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
    call f(int(size(x), kind=c_int), x, int(shape(u), kind=c_int), u, v, env)
  end function wF
end subroutine Lib_BLAS_SFuncProd$rank
${writeIncLine( &
fr'''EXTERN void _Lib_BLAS_SFuncProd{rank}( &
    tFieldBase*, tFieldBase*, tSMFuncPtr, void*); &
template<typename tSMFunc> &
void BLAS_SFuncProd(tField<{rank}> v, &
                    tField<{rank}> u, tSMFunc&& func) {{ &
  _Lib_BLAS_SFuncProd{rank}(v.Base(), u.Base(), &
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
subroutine Lib_FDM_Gradient$rank(pV_bar, lambda, pU, dir) &
  & bind(c, name='_Lib_FDM_Gradient$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV_bar
  character(c_char), intent(in), value :: dir
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v_bar(:,@:,:)
  u => DEREF$rank(pU); v_bar => DEREF${rank+1}$(pV_bar)
  ! ----------------------
  !select case(dir)
  ! case('c')
      call FDM_Gradient_Central(gMesh, v_bar, lambda, u)
  !  case('f')
  !    call FDM_Gradient_Forward(gMesh, v_bar, lambda, u, dirAll=1_1)
  !  case('b')
  !    call FDM_Gradient_Forward(gMesh, v_bar, lambda, u, dirAll=-1_1)
  !end select
end subroutine Lib_FDM_Gradient$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FDM_Gradient{rank}( &
    tFieldBase*, double, tFieldBase*, char); &
void FDM_Gradient(tField<{rank+1}> v_bar, &
                  double lambda, tField<{rank}> u, char dir) {{ &
  _Lib_FDM_Gradient{rank}(v_bar.Base(), lambda, u.Base(), dir); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS-1
subroutine Lib_FDM_Divergence$rank(pV, lambda, pU_bar, dir) &
  & bind(c, name='_Lib_FDM_Divergence$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU_bar, pV
  character(c_char), intent(in), value :: dir
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u_bar(:,@:,:), v(@:,:)
  u_bar => DEREF${rank+1}$(pU_bar); v => DEREF$rank(pV)
  ! ----------------------
  !select case(dir)
  !  case('c')
      call FDM_Divergence_Central(gMesh, v, lambda, u_bar)
  !  case('f')
  !    call FDM_Divergence_Backward(gMesh, v, lambda, u_bar, dirAll=1_1)
  !  case('b')
  !    call FDM_Divergence_Backward(gMesh, v, lambda, u_bar, dirAll=-1_1)
  !end select
end subroutine Lib_FDM_Divergence$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FDM_Divergence{rank}( &
    tFieldBase*, double, tFieldBase*, char); &
void FDM_Divergence(tField<{rank}> v, &
                    double lambda, tField<{rank+1}> u_bar, char dir) {{ &
  _Lib_FDM_Divergence{rank}(v.Base(), lambda, u_bar.Base(), dir); &
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
  _Lib_FDM_Laplacian{rank}(v.Base(), lambda, u.Base()); &
}}''')}$
#$end do
${writeIncLine('')}$

#$do rank = 0, NUM_RANKS
subroutine Lib_FDM_DivWGrad$rank(pV, lambda, pW, pU) &
  & bind(c, name='_Lib_FDM_DivWGrad$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV, pW
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:), w(@:,:)
  u => DEREF$rank(pU); v => DEREF$rank(pV); w => DEREF$rank(pW)
  ! ----------------------
  call FDM_DivWGrad_Central(gMesh, v, lambda, w, u)
end subroutine Lib_FDM_DivWGrad$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FDM_DivWGrad{rank}( &
    tFieldBase*, double, tFieldBase*, tFieldBase*); &
void FDM_DivWGrad(tField<{rank}> v, &
                  double lambda, tField<{rank}> w, tField<{rank}> u) {{ &
  _Lib_FDM_DivWGrad{rank}(v.Base(), lambda, w.Base(), u.Base()); &
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
subroutine Lib_FDM_Convection$rank(pV, lambda, pU, pW_bar) &
  & bind(c, name='_Lib_FDM_Convection$rank')
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV, pW_bar
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), pointer :: u(@:,:), v(@:,:), w_bar(:,:)
  u => DEREF$rank(pU); v => DEREF$rank(pV); w_bar => DEREF$1(pW_bar)
  ! ----------------------
  call FDM_Convection_Central(gMesh, v, lambda, u, w_bar)
end subroutine Lib_FDM_Convection$rank
${writeIncLine( &
fr'''EXTERN void _Lib_FDM_Convection{rank}( &
    tFieldBase*, double, tFieldBase*, tFieldBase*); &
void FDM_Convection(tField<{rank}> v, &
                    double lambda, tField<{rank}> u, tField<1> w_bar) {{ &
  _Lib_FDM_Convection{rank}(v.Base(), lambda, u.Base(), w_bar.Base()); &
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
  call Params%Init(1.0D-8, &
    &              1.0D-8, 100000)
  call Solve_CG(gMesh, u, b, wA, Params, Params)
contains
  subroutine wA(mesh, v, w, opParams)
    class(tMesh), intent(in) :: mesh
    real(dp), intent(in), target :: w(@:,:)
    real(dp), intent(inout), target :: v(@:,:)
    class(*), intent(in) :: opParams
    type(tLibFieldBase$rank), target :: fV, fW
    type(c_ptr) :: pV, pW
    fV%Values => v; fW%Values => w
    pV = c_loc(fV); pW = c_loc(fW)
    call A(pV, pW, env)
  end subroutine wA
end subroutine Lib_Solve_BiCGStab$rank
${writeIncLine( &
fr'''EXTERN void _Lib_Solve_BiCGStab{rank}( &
    tFieldBase*, tFieldBase*, tMeshOperatorPtr, void*); &
template<typename tMeshOperator> &
void Solve_BiCGStab(tField<{rank}> u, &
                    tField<{rank}> b, tMeshOperator&& meshOperator) {{ &
  _Lib_Solve_BiCGStab{rank}(u.Base(), b.Base(), &
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
