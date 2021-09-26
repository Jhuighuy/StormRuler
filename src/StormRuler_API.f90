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
module StormRuler_API

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8
use StormRuler_Helpers
use StormRuler_IO, only: tIOList => IOList ! TODO:
use StormRuler_IO
use StormRuler_Mesh, only: tMesh

use StormRuler_BLAS, only: Fill, Fill_Random, Set, &
  & Scale, Add, Sub, FuncProd, SFuncProd

#$for type_, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$type_
use StormRuler_Solvers_Precond, only: tPrecondFunc$type_
#$end for

use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_CG, only: Solve_CG, Solve_BiCGStab
use StormRuler_Solvers_Chebyshev, only: Solve_Chebyshev
use StormRuler_Solvers_MINRES, only: Solve_MINRES

use StormRuler_FDM_Operators, only: &
  & FDM_Convection_Central ! TODO: should be StormRuler_FDM_Convection
use StormRuler_FDM_Operators, only: &
  & FDM_Gradient_Central, FDM_Divergence_Central, &
  & FDM_Laplacian_Central, FDM_DivWGrad_Central
use StormRuler_FDM_BCs, only: &
  & FDM_ApplyBCs

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: &
  & c_char, c_int, c_size_t, c_ptr, c_funptr, &
  & c_loc, c_f_pointer, c_f_procpointer   

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

#$let CLASSES = ['tMesh', 'tIOList']

!! ----------------------------------------------------------------- !!
!! Class wrapper struct.
!! ----------------------------------------------------------------- !!
#$for klass in CLASSES
type :: ${klass}$_C
  class($klass), pointer :: mObj
end type !${klass}$_C
#$end for

!! ----------------------------------------------------------------- !!
!! Field wrapper class.
!! ----------------------------------------------------------------- !!
#$for type, typename in SCALAR_TYPES
type :: tField_C$type
  character :: mType = '$type'
  $typename, pointer :: mData(:,:) => null()
end type !tField_C$type
#$end for

interface Wrap
#$for klass in CLASSES
  module procedure Wrap$klass
#$end for
end interface

interface Unwrap
  module procedure UnwrapString
#$for klass in CLASSES
  module procedure Unwrap$klass
#$end for
#$for type, typename in SCALAR_TYPES
  module procedure Unwrap2D$type
#$end for
end interface Unwrap

enum, bind(C)
  enumerator :: SR_Mat_SymmDefinite = 100
  enumerator :: SR_Mat_SymmSemiDefinite, SR_Mat_Symm
  enumerator :: SR_Mat_General, SR_Mat_General_Singular
end enum

enum, bind(C)
  enumerator :: SR_Auto = 200
  enumerator :: SR_CG, SR_BiCGStab
  enumerator :: SR_Cheby, SR_ChebyCG
  enumerator :: SR_MINRES, SR_GMRES
end enum

enum, bind(C)
  enumerator :: SR_Precond_None = 300
  enumerator :: SR_Precond_Jacobi, SR_Precond_LU_SGS
end enum

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Unwrap a class pointer.
!! ----------------------------------------------------------------- !!
#$for klass in CLASSES
function Wrap$klass(obj) result(pObj)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class($klass), intent(in), target :: obj
  type(c_ptr) :: pObj
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(${klass}$_C), pointer :: pObj_C

  allocate(pObj_C)
  pObj_C%mObj => obj
  pObj = c_loc(pObj_C)

end function Wrap$klass
#$end for

!! ----------------------------------------------------------------- !!
!! Unwrap a C-string pointer.
!! ----------------------------------------------------------------- !!
subroutine UnwrapString(string, pString)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pString
  character(len=:), intent(out), pointer :: string
  ! >>>>>>>>>>>>>>>>>>>>>>

  interface
    pure function strlen(pString) bind(C, name='strlen')
      import :: c_ptr, c_size_t
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(c_ptr), intent(in), value :: pString
      integer(c_size_t) :: strlen
      ! >>>>>>>>>>>>>>>>>>>>>>
    end function strlen
  end interface

  integer(c_size_t) :: length
  
  ! Compute C-string length and allocate Fortran string.
  length = strlen(pString)
  allocate( character(len=length) :: string )

  ! Cast C-pointer and copy characters.
  block
    character(kind=c_char, len=length), pointer :: string_C
    call c_f_pointer(cptr=pString, fptr=string_C)
    string(:) = string_C(:)
  end block
end subroutine UnwrapString

!! ----------------------------------------------------------------- !!
!! Unwrap a class pointer.
!! ----------------------------------------------------------------- !!
#$for klass in CLASSES
subroutine Unwrap$klass(obj, pObj, free)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pObj
  class($klass), intent(out), pointer :: obj
  logical, intent(in), optional :: free
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(${klass}$_C), pointer :: pObj_C

  call c_f_pointer(cptr=pObj, fptr=pObj_C)
  obj => pObj_C%mObj

  if (present(free)) then
    if (free) deallocate(pObj_C)
  end if

end subroutine Unwrap$klass
#$end for

#$for type, typename in SCALAR_TYPES
subroutine Unwrap2D$type(x, pX, free)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  $typename, intent(out), pointer :: x(:,:)
  logical, intent(in), optional :: free
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(tField_C$type), pointer :: pX_C

  call c_f_pointer(cptr=pX, fptr=pX_C)
  x => pX_C%mData
  !x => Reshape2D(pX_C%mData)

  if (present(free)) then
    if (free) deallocate(pX_C)
  end if

end subroutine Unwrap2D$type
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

function InitMesh() result(pMesh) &
    & bind(C, name='SR_InitMesh')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr) :: pMesh
  ! >>>>>>>>>>>>>>>>>>>>>>
  class(tMesh), pointer :: gMesh
  type(tMesh_C), pointer :: pMesh_C
  real(8), parameter :: pi = 4*atan(1.0D0)
  integer(ip), Parameter :: Nx = 100, Ny = 100
  Real(8), Parameter :: Dx = 2*pi/Nx, Dy = 2*pi/Ny, Dt = Dx*Dx
#$if False
  allocate(gMesh)
  gMesh%dt = dt
  call gMesh%InitRect(dx, nx, .true., dy, ny, .true., 20)
#$else
  integer(ip), allocatable :: pixels(:,:)
  integer(ip), allocatable :: colorToBCM(:)
  allocate(gMesh)
  call Load_PPM('test/Domain-100.ppm', pixels)
  colorToBCM = [PixelToInt([255, 255, 255]), PixelToInt([255, 0, 0])]
  call gMesh%InitFromImage2D(pixels, 0, colorToBCM, 10)
  allocate(gMesh%dr(1:2, 1:4))
  gMesh%dl = [Dx,Dx,Dy,Dy]
  gMesh%dr(:,1) = [gMesh%dl(1), 0.0_dp]
  gMesh%dr(:,2) = [gMesh%dl(1), 0.0_dp]
  gMesh%dr(:,3) = [0.0_dp, gMesh%dl(3)]
  gMesh%dr(:,4) = [0.0_dp, gMesh%dl(3)]
  call gMesh%PrintTo_Neato('test/c2c.dot')
  call gMesh%PrintTo_LegacyVTK('test/c2c.vtk')
#$endif
  call gMesh%SetRange()
  allocate(pMesh_C)
  pMesh_C%mObj => gMesh
  pMesh = c_loc(pMesh_C)
  !L = 0
end function InitMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
function Alloc$type(pMesh, numVars, rank) result(pY) &
    & bind(C, name='SR_Alloc$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  integer(c_int), intent(in), value :: numVars, rank
  type(c_ptr) :: pY
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  type(tField_C$type), pointer :: pY_C

  call Unwrap(mesh, pMesh)
  allocate(pY_C)
  allocate(pY_C%mData( (mesh%Dim**rank)*numVars, mesh%NumAllCells ))
  pY = c_loc(pY_C)

end function Alloc$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
function Alloc_Mold$type(pX) result(pY) &
    & bind(C, name='SR_Alloc_Mold$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  type(c_ptr) :: pY
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(tField_C$type), pointer :: pY_C
  $typename, pointer :: x(:,:)

  call Unwrap(x, pX)
  allocate(pY_C)
  allocate(pY_C%mData, mold=x)
  pY = c_loc(pY_C)

end function Alloc_Mold$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
subroutine Free$type(pX) &
    & bind(C, name='SR_Free$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  ! >>>>>>>>>>>>>>>>>>>>>>

  $typename, pointer :: x(:,:)

  call Unwrap(x, pX, free=.true.)
  deallocate(x)

end subroutine Free$type
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
function IO_Begin() result(pIOList) &
    & bind(C, name='SR_IO_Begin') 
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr) :: pIOList
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(tIOList_C), pointer :: pIOList_C

  allocate(pIOList_C)
  allocate(pIOList_C%mObj)
  pIOList = c_loc(pIOList_C)

end function IO_Begin

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
subroutine IO_Add(pIOList, pX, pName) &
  & bind(C, name='SR_IO_Add') 
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pIOList
  type(c_ptr), intent(in), value :: pX, pName
  ! >>>>>>>>>>>>>>>>>>>>>>

  character(len=:), pointer :: name

  call Unwrap(name, pName)

end subroutine IO_Add

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
subroutine Fill$type(pMesh, pX, alpha, beta) &
    & bind(C, name='SR_Fill$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX
  real(dp), intent(in), value :: alpha
  $typename, intent(in), value :: beta
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX)

  call Fill(mesh, x, alpha, beta)

end subroutine Fill$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
subroutine Fill_Random$type(pMesh, pX, alpha, beta) &
    & bind(C, name='SR_Fill_Random$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX
  real(dp), intent(in), value :: alpha, beta
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX)

  call Fill_Random(mesh, x, alpha, beta)

end subroutine Fill_Random$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
subroutine Set$type(pMesh, pY, pX) &
    & bind(C, name='SR_Set$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY)

  call Set(mesh, y, x)

end subroutine Set$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
subroutine Scale$type(pMesh, pY, pX, alpha) &
    & bind(C, name='SR_Scale$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  $typename, intent(in), value :: alpha
  type(c_ptr), intent(in), value :: pX, pY
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY)

  call Scale(mesh, y, x, alpha)

end subroutine Scale$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
subroutine Add$type(pMesh, pZ, pY, pX, alpha, beta) &
    & bind(C, name='SR_Add$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  $typename, intent(in), value :: alpha, beta
  type(c_ptr), intent(in), value :: pX, pY, pZ
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:), z(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY); call Unwrap(z, pZ)

  call Add(mesh, z, y, x, alpha, beta)

end subroutine Add$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
subroutine FuncProd$type(pMesh, pY, pX, pF, env) &
    & bind(C, name='SR_FuncProd$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: env
  ! >>>>>>>>>>>>>>>>>>>>>>

  abstract interface
    pure subroutine tMapFunc_C(size, Fx, x, env) bind(C)
      import c_int, c_ptr, dp
      ! <<<<<<<<<<<<<<<<<<<<<<
      integer(c_int), intent(in), value :: size
      $typename, intent(in) :: x(*)
      $typename, intent(inout) :: Fx(*)
      type(c_ptr), intent(in), value :: env
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine tMapFunc_C
  end interface

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:)
  procedure(tMapFunc_C), pointer :: f

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY)
  call c_f_procpointer(cptr=pF, fptr=f)

  call FuncProd(mesh, y, x, f_wrapper)

contains
  pure function f_wrapper(x) result(Fx)
    ! <<<<<<<<<<<<<<<<<<<<<<
    $typename, intent(in) :: x(:)
    $typename :: Fx(size(x))
    ! >>>>>>>>>>>>>>>>>>>>>>
    
    call f(int(size(x), kind=c_int), Fx, x, env)

  end function f_wrapper
end subroutine FuncProd$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in SCALAR_TYPES
subroutine SFuncProd$type(pMesh, pY, pX, pF, env) &
    & bind(C, name='SR_SFuncProd$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: env
  ! >>>>>>>>>>>>>>>>>>>>>>

  abstract interface
    pure subroutine tSMapFunc_C(dim, r, size, Fx, x, env) bind(C)
      import c_int, c_ptr, dp
      ! <<<<<<<<<<<<<<<<<<<<<<
      integer(c_int), intent(in), value :: dim, size
      real(dp), intent(in) :: r(*)
      $typename, intent(in) :: x(*)
      $typename, intent(inout) :: Fx(*)
      type(c_ptr), intent(in), value :: env
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine tSMapFunc_C
  end interface

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:)
  procedure(tSMapFunc_C), pointer :: f

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY)
  call c_f_procpointer(cptr=pF, fptr=f)

  call SFuncProd(mesh, y, x, f_wrapper)

contains
  pure function f_wrapper(r, x) result(Fx)
    ! <<<<<<<<<<<<<<<<<<<<<<
    real(dp), intent(in) :: r(:)
    $typename, intent(in) :: x(:)
    $typename :: Fx(size(x))
    ! >>>>>>>>>>>>>>>>>>>>>>
    
    call f(int(size(r), kind=c_int), r, &
      & int(size(x), kind=c_int), Fx, x, env)

  end function f_wrapper
end subroutine SFuncProd$type
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in [SCALAR_TYPES[0]]
subroutine LinSolve$type(pMesh, pX, pB, pMatVec, env, &
    & Solver, Precond, pMatVec_H, env_H) bind(C, name='SR_LinSolve$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pB
  type(c_funptr), intent(in), value :: pMatVec, pMatVec_H
  type(c_ptr), intent(in), value :: env, env_H
  integer(c_int), intent(in), value :: Solver, Precond
  ! >>>>>>>>>>>>>>>>>>>>>>

  abstract interface
    subroutine tMatVec_C(pMesh, pAx, pX, env) bind(C)
      import :: c_ptr, dp
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(c_ptr), intent(in), value :: pMesh
      type(c_ptr), intent(in), value :: pAx, pX
      type(c_ptr), intent(in), value :: env
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine tMatVec_C
  end interface

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), b(:,:)
  procedure(tMatVec_C), pointer :: MatVec
  type(tConvParams) :: Params

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(b, pB)
  call c_f_procpointer(cptr=pMatVec, fptr=MatVec)

  call Params%Init(1.0D-8, 1.0D-8, 2000, 'CG')
  call Solve_CG(mesh, x, b, MatVec_wrapper, Params, Params)

contains
  subroutine MatVec_wrapper(mesh_, Ax, x, env_)
    ! <<<<<<<<<<<<<<<<<<<<<<
    class(tMesh), intent(in) :: mesh_
    $typename, intent(in), target :: x(:,:)
    $typename, intent(inout), target :: Ax(:,:)
    class(*), intent(inout) :: env_
    ! >>>>>>>>>>>>>>>>>>>>>>

    !type(tMesh_C), target :: pMesh_C
    !type(c_ptr) :: pMesh,
    type(tField_C$type), target :: pX_C, pAx_C
    type(c_ptr) :: pX, pAx

    !pMesh_C%mMesh => mesh
    !pMesh = c_loc(pMesh_C)
    pX_C%mData => x; pAx_C%mData => Ax
    pX = c_loc(pX_C); pAx = c_loc(pAx_C)

    call MatVec(pMesh, pAx, pX, env)

  end subroutine MatVec_wrapper
end subroutine LinSolve$type
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in [SCALAR_TYPES[0]]
subroutine ApplyBCs$type(pMesh, pU, BCmask, alpha, beta, gamma) &
    & bind(C, name='SR_ApplyBCs$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  integer(c_int), intent(in), value :: BCmask
  $typename, intent(in), value :: alpha, beta, gamma
  type(c_ptr), intent(in), value :: pU
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: u(:,:)
  integer(ip) :: iBC, firstBC, lastBC

  call Unwrap(mesh, pMesh)
  call Unwrap(u, pU)

  if (BCmask /= 0) then
    firstBC = BCmask; lastBC = BCmask
  else
    firstBC = 1; lastBC = mesh%NumBCMs
  end if
  do iBC = firstBC, lastBC
    call FDM_ApplyBCs(mesh, iBC, u, alpha, beta, gamma)
  end do

end subroutine ApplyBCs$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in [SCALAR_TYPES[0]]
subroutine Conv$type(pMesh, pV, lambda, pU, pA) &
    & bind(C, name='SR_Conv$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  $typename, intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV, pA
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: u(:,:), v(:,:), a(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(u, pU); call Unwrap(v, pV); call Unwrap(a, pA)

  call FDM_Convection_Central(mesh, v, lambda, u, a)

end subroutine Conv$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in [SCALAR_TYPES[0]]
subroutine DivGrad$type(pMesh, pV, lambda, pU) &
    & bind(C, name='SR_DivGrad$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  $typename, intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: u(:,:), v(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(u, pU); call Unwrap(v, pV)

  call FDM_Laplacian_Central(mesh, v, lambda, u)

end subroutine DivGrad$type
#$end for

!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
!! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ !!
#$for type, typename in [SCALAR_TYPES[0]]
subroutine DivKGrad$type(pMesh, pV, lambda, pK, pU) &
    & bind(C, name='SR_DivKGrad$type')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  $typename, intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV, pK
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: u(:,:), v(:,:), k(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(u, pU); call Unwrap(v, pV); call Unwrap(k, pK)

  call FDM_DivWGrad_Central(mesh, v, lambda, k, u)

end subroutine DivKGrad$type
#$end for

end module StormRuler_API
