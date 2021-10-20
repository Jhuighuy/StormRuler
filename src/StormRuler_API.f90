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

use StormRuler_Parameters, only: dp, ip, i8, error_code
use StormRuler_Helpers
use StormRuler_IO, only: tIOList => IOList ! TODO:
use StormRuler_IO
use StormRuler_Mesh, only: tMesh

use StormRuler_BLAS, only: Norm_2, &
  & Fill, Fill_Random, Set, &
  & Scale, Add, Sub, Mul, FuncProd, SFuncProd
#$for type_, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$type_
#$end for

use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_CG, only: Solve_CG, Solve_BiCGStab
use StormRuler_Solvers_Chebyshev, only: Solve_Chebyshev
use StormRuler_Solvers_MINRES, only: Solve_MINRES
use StormRuler_Solvers_GMRES, only: Solve_GMRES
use StormRuler_Solvers_LSQR, only: Solve_LSQR, Solve_LSMR
#$for type_, _ in SCALAR_TYPES
use StormRuler_Solvers_Precond, only: tPrecondFunc$type_
#$end for
use StormRuler_Solvers_Precond, only: &
  & Precondition_Jacobi, Precondition_LU_SGS

use StormRuler_FDM_BCs, only: FDM_ApplyBCs
use StormRuler_FDM_Operators, only: &
  & FDM_Gradient_Central, FDM_Divergence_Central, &
  & FDM_Laplacian_Central, FDM_DivWGrad_Central
  use StormRuler_FDM_Operators, only: &
  & FDM_Convection_Central ! TODO: should be StormRuler_FDM_Convection

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: c_char, c_int, c_long_long, &
  & c_size_t, c_ptr, c_funptr, c_null_ptr, c_null_funptr, &
  & c_loc, c_funloc, c_f_pointer, c_f_procpointer, c_associated

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface
  subroutine PrintPointer(p) bind(C, name='SR_PrintPointer')
    import :: c_ptr
    ! <<<<<<<<<<<<<<<<<<<<<<
    type(c_ptr), intent(in), value :: p
    ! >>>>>>>>>>>>>>>>>>>>>>
  end subroutine PrintPointer
end interface

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
#$for T, typename in SCALAR_TYPES
type :: tField_C$T
  character :: mType = '$T'
  integer(ip) :: mRank = -1
  $typename, pointer :: mData(:,:) => null()
end type !tField_C$T
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
#$for T, typename in SCALAR_TYPES
  module procedure UnwrapField$T
  module procedure UnwrapVecField$T
#$end for
end interface Unwrap

type(c_ptr), parameter :: ceMatVec = transfer(10000, c_null_ptr)
type(c_ptr), parameter :: ceMatVec_H = transfer(20000, c_null_ptr)

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
!! Unwrap a string pointer.
!! ----------------------------------------------------------------- !!
subroutine UnwrapString(string, pString)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pString
  character(len=:), intent(out), pointer :: string
  ! >>>>>>>>>>>>>>>>>>>>>>

  interface
    pure function cStrlen(pString) bind(C, name='strlen')
      import :: c_ptr, c_size_t
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(c_ptr), intent(in), value :: pString
      integer(c_size_t) :: cStrlen
      ! >>>>>>>>>>>>>>>>>>>>>>
    end function cStrlen
  end interface

  integer(c_size_t) :: length
  
  length = cStrlen(pString)
  allocate( character(len=length) :: string )

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

!! ----------------------------------------------------------------- !!
!! ----------------------------------------------------------------- !!
#$for T, typename in SCALAR_TYPES
subroutine UnwrapField$T(x, pX, free, pX_C_out, rank)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  $typename, intent(out), pointer :: x(:,:)
  logical, intent(in), optional :: free
  type(tField_C$T), intent(out), pointer, optional :: pX_C_out
  integer(ip), intent(out), optional :: rank
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(tField_C$T), pointer :: pX_C

  call c_f_pointer(cptr=pX, fptr=pX_C)
  if (present(pX_C_out)) pX_C_out => pX_C
  x => pX_C%mData
  !x => Reshape2D(pX_C%mData)

  if (present(free)) then
    if (free) deallocate(pX_C)
  end if

  if (present(rank)) then
    rank = pX_C%mRank
  end if

end subroutine UnwrapField$T
#$end for

!! ----------------------------------------------------------------- !!
!! ----------------------------------------------------------------- !!
#$for T, typename in SCALAR_TYPES
subroutine UnwrapVecField$T(x, pX, dim, free, pX_C_out)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  $typename, intent(out), pointer :: x(:,:,:)
  integer(ip), intent(in) :: dim
  logical, intent(in), optional :: free
  type(tField_C$T), intent(out), pointer, optional :: pX_C_out
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(tField_C$T), pointer :: pX_C

  call c_f_pointer(cptr=pX, fptr=pX_C)
  if (present(pX_C_out)) pX_C_out => pX_C
  !x => pX_C%mData
#$if T == 'R'
  x => AsVecField(dim, pX_C%mData)
#$end if

  if (present(free)) then
    if (free) deallocate(pX_C)
  end if

end subroutine UnwrapVecField$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
function cInitMesh() result(pMesh) bind(C, name='SR_InitMesh')
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
  call Load_PPM('test/Domain-100-Tube.ppm', pixels)
  !colorToBCM = [PixelToInt([255, 255, 255]), PixelToInt([255, 0, 0])]
  colorToBCM = [PixelToInt([255, 255, 255]), PixelToInt([255, 0, 0]), &
    & PixelToInt([0, 255, 0]), PixelToInt([0, 0, 255])]
  call gMesh%InitFromImage2D(pixels, 0, colorToBCM, 1)
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
end function cInitMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
function cAlloc$T(pMesh, numVars, rank) result(pY)  bind(C, name='SR_Alloc$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  integer(c_int), intent(in), value :: numVars, rank
  type(c_ptr) :: pY
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  type(tField_C$T), pointer :: pY_C

  call Unwrap(mesh, pMesh)
  allocate(pY_C)
  allocate(pY_C%mData( (mesh%NumDims**rank)*numVars, mesh%NumAllCells ))
  pY_C%mRank = rank
  pY = c_loc(pY_C)

end function cAlloc$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
function cAlloc_Mold$T(pX) result(pY) bind(C, name='SR_Alloc_Mold$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  type(c_ptr) :: pY
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(tField_C$T), pointer :: pX_C, pY_C
  $typename, pointer :: x(:,:)

  call Unwrap(x, pX, pX_C_out=pX_C)
  allocate(pY_C)
  allocate(pY_C%mData, mold=x)
  pY_C%mRank = pX_C%mRank
  pY = c_loc(pY_C)

end function cAlloc_Mold$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
subroutine cFree$T(pX) bind(C, name='SR_Free$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  ! >>>>>>>>>>>>>>>>>>>>>>

  $typename, pointer :: x(:,:)

  call Unwrap(x, pX, free=.true.)
  deallocate(x)

end subroutine cFree$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
integer(c_int) function cIsVecField$T(pX) bind(C, name='SR_IsVecField$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: x(:,:)
  integer :: rank

  call Unwrap(x, pX, rank=rank)
  cIsVecField$T = merge(1, 0, rank == 2)

end function cIsVecField$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
integer(c_int) function cIsMatField$T(pX) bind(C, name='SR_IsMatField$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pX
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: x(:,:)
  integer :: rank

  call Unwrap(x, pX, rank=rank)
  cIsMatField$T = merge(1, 0, rank == 3)

end function cIsMatField$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
function cIO_Begin() result(pIOList) bind(C, name='SR_IO_Begin') 
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr) :: pIOList
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(tIOList_C), pointer :: pIOList_C

  allocate(pIOList_C)
  allocate(pIOList_C%mObj)
  pIOList = c_loc(pIOList_C)

end function cIO_Begin

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
subroutine cIO_Add(pIOList, pX, pName) bind(C, name='SR_IO_Add') 
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pIOList
  type(c_ptr), intent(in), value :: pX, pName
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tIOList), pointer :: ioList
  type(tField_CR), pointer :: pX_C
  real(dp), pointer :: x(:,:)
  character(len=:), pointer :: name

  call Unwrap(ioList, pIOList)
  call Unwrap(name, pName)
  call Unwrap(x, pX, pX_C_out=pX_C)

  if (pX_C%mRank == 0) then
    call ioList%Add(name, x(1,:))
  else
    call ioList%Add(name, x)
  end if

end subroutine cIO_Add

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
subroutine cIO_Flush(pIOList, pMesh, pFileName) bind(C, name='SR_IO_Flush') 
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pIOList
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pFileName
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tIOList), pointer :: ioList
  class(tMesh), pointer :: mesh
  character(len=:), pointer :: filename

  call Unwrap(ioList, pIOList, free=.true.)
  call Unwrap(mesh, pMesh)
  call Unwrap(filename, pFileName)

  call mesh%PrintTo_LegacyVTK(filename, ioList)
  deallocate(ioList)

end subroutine cIO_Flush

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
subroutine cFill$T(pMesh, pX, alpha, beta) bind(C, name='SR_Fill$T')
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

end subroutine cFill$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
subroutine cFill_Random$T(pMesh, pX, alpha, beta) bind(C, name='SR_Fill_Random$T')
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

end subroutine cFill_Random$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
subroutine cSet$T(pMesh, pY, pX) bind(C, name='SR_Set$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY)

  call Set(mesh, y, x)

end subroutine cSet$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
subroutine cScale$T(pMesh, pY, pX, alpha) bind(C, name='SR_Scale$T')
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

end subroutine cScale$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
subroutine cAdd$T(pMesh, pZ, pY, pX, alpha, beta) bind(C, name='SR_Add$T')
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

end subroutine cAdd$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in {SCALAR_TYPES[0]}
subroutine cMul$T(pMesh, pZ, pY, pX) bind(C, name='SR_Mul$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY, pZ
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:), z(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY); call Unwrap(z, pZ)

  call Mul(mesh, z, y(1,:), x)

end subroutine cMul$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
subroutine cFuncProd$T(pMesh, pY, pX, pF, pEnv) bind(C, name='SR_FuncProd$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: pEnv
  ! >>>>>>>>>>>>>>>>>>>>>>

  abstract interface
    pure subroutine ctMapFunc(size, Fx, x, pEnv) bind(C)
      import :: c_int, c_ptr, dp
      ! <<<<<<<<<<<<<<<<<<<<<<
      integer(c_int), intent(in), value :: size
      $typename, intent(in) :: x(*)
      $typename, intent(inout) :: Fx(*)
      type(c_ptr), intent(in), value :: pEnv
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine ctMapFunc
  end interface

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:)
  procedure(ctMapFunc), pointer :: f

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY)
  call c_f_procpointer(cptr=pF, fptr=f)

  call FuncProd(mesh, y, x, cF)

contains
  pure function cF(x) result(Fx)
    ! <<<<<<<<<<<<<<<<<<<<<<
    $typename, intent(in) :: x(:)
    $typename :: Fx(size(x))
    ! >>>>>>>>>>>>>>>>>>>>>>
    
    call f(int(size(x), kind=c_int), Fx, x, pEnv)

  end function cF
end subroutine cFuncProd$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in SCALAR_TYPES
subroutine cSFuncProd$T(pMesh, pY, pX, pF, pEnv) bind(C, name='SR_SFuncProd$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: pEnv
  ! >>>>>>>>>>>>>>>>>>>>>>

  abstract interface
    pure subroutine ctSMapFunc(dim, r, size, Fx, x, pEnv) bind(C)
      import :: c_int, c_ptr, dp
      ! <<<<<<<<<<<<<<<<<<<<<<
      integer(c_int), intent(in), value :: dim, size
      real(dp), intent(in) :: r(*)
      $typename, intent(in) :: x(*)
      $typename, intent(inout) :: Fx(*)
      type(c_ptr), intent(in), value :: pEnv
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine ctSMapFunc
  end interface

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), y(:,:)
  procedure(ctSMapFunc), pointer :: f

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(y, pY)
  call c_f_procpointer(cptr=pF, fptr=f)

  call SFuncProd(mesh, y, x, cF)

contains
  pure function cF(r, x) result(Fx)
    ! <<<<<<<<<<<<<<<<<<<<<<
    real(dp), intent(in) :: r(:)
    $typename, intent(in) :: x(:)
    $typename :: Fx(size(x))
    ! >>>>>>>>>>>>>>>>>>>>>>
    
    call f(int(size(r), kind=c_int), r, &
      & int(size(x), kind=c_int), Fx, x, pEnv)

  end function cF
end subroutine cSFuncProd$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine cLinSolve$T(pMesh, pX, pB, pMatVec, pEnv, &
    & solverType, precondType, pMatVec_H, pEnv_H) !bind(C, name='SR_LinSolve$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pB
  type(c_funptr), intent(in), value :: pMatVec, pMatVec_H
  type(c_ptr), intent(in), value :: pEnv, pEnv_H
  integer(c_int), intent(in), value :: solverType, precondType
  ! >>>>>>>>>>>>>>>>>>>>>>

  abstract interface
    subroutine ctMatVec(pMesh, pAx, pX, pEnv) !bind(C)
      import :: c_ptr, dp
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(c_ptr), intent(in), value :: pMesh
      type(c_ptr), intent(in), value :: pAx, pX
      type(c_ptr), intent(in), value :: pEnv
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine ctMatVec
  end interface

  enum, bind(C)
    enumerator :: eAuto = 200, eCG, eBiCGStab, &
      & eCheby, eChebyCG, eMINRES, eGMRES, eLSQR, eLSMR
  end enum

  enum, bind(C)
    enumerator :: ePrecond_None = 300, &
      & ePrecond_Jacobi, ePrecond_LU_SGS
  end enum

  class(tMesh), pointer :: mesh
  $typename, pointer :: x(:,:), b(:,:), z(:,:), Az(:,:), f(:,:)
  procedure(ctMatVec), pointer :: MatVec, MatVec_H
  type(tConvParams) :: params
  procedure(tPrecondFunc$T), pointer :: Precond, Precond_T

  ! Select the preconditioner.
  Precond => null()
  Precond_T => null()
  select case(precondType)
    case(ePrecond_None)
    case(ePrecond_Jacobi)
      Precond => Precondition_Jacobi
      Precond_T => Precondition_Jacobi
    case(ePrecond_LU_SGS)
      Precond => Precondition_LU_SGS
      Precond_T => Precondition_LU_SGS
    case default
      error stop error_code
  end select

  call Unwrap(mesh, pMesh)
  call Unwrap(x, pX); call Unwrap(b, pB)
  call c_f_procpointer(cptr=pMatVec, fptr=MatVec)
  if (.not.c_associated(pMatVec_H, c_null_funptr)) then
    call c_f_procpointer(cptr=pMatVec_H, fptr=MatVec_H)
  end if

  ! Fix the possible non-uniform operator.
  allocate(z, Az, f, mold=x)
  call Fill(mesh, z, 0.0_dp)
  call Fill(mesh, Az, 0.0_dp)
  call cMatVec(mesh, f, z, Params)
  call Set(mesh, Az, f)
  call Sub(mesh, f, b, Az)

  ! Launch the solver.
  call params%Init(1.0D-8, 1.0D-8, 2000)
  select case(solverType)
    case(eCG)
      params%Name = 'CG'
      call Solve_CG(mesh, x, f, cMatVec, params, params)
    case(eBiCGStab)
      params%Name = 'BiCGStab'
      call Solve_BiCGStab(mesh, x, f, cMatVec, params, params)
    case(eMINRES)
      params%Name = 'MINRES'
      call Solve_MINRES(mesh, x, f, cMatVec, params, params)
    case(eGMRES)
      params%Name = 'GMRES'
      call Solve_GMRES(mesh, x, f, cMatVec, params, params)
    case(eLSQR)
      params%Name = 'LSQR'
      call Solve_LSQR(mesh, x, f, cMatVec, params, params)
    case(eLSMR)
      params%Name = 'LSMR'
      call Solve_LSMR(mesh, x, f, cMatVec, params, params)
    case default
      error stop error_code
  end select

contains
  subroutine cMatVec(mesh_, Ax, x, env_)
    ! <<<<<<<<<<<<<<<<<<<<<<
    class(tMesh), intent(in) :: mesh_
    $typename, intent(in), target :: x(:,:)
    $typename, intent(inout), target :: Ax(:,:)
    class(*), intent(inout) :: env_
    ! >>>>>>>>>>>>>>>>>>>>>>

    type(tField_C$T), target :: pX_C, pAx_C
    type(c_ptr) :: pX, pAx

    pX_C%mData => x; pAx_C%mData => Ax
    pX = c_loc(pX_C); pAx = c_loc(pAx_C)

    call MatVec(pMesh, pAx, pX, pEnv)
    call Sub(mesh, Ax, Ax, Az)

  end subroutine cMatVec
end subroutine cLinSolve$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in [SCALAR_TYPES[0]]
function cRCI_LinSolve$T(pMesh, pX, pB, solverType, precondType, &
    & pAy, pY) result(request) bind(C, name='SR_RCI_LinSolve$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pB
  integer(c_int), intent(in), value :: solverType, precondType
  type(c_ptr), intent(inout), target :: pAy, pY
  integer(c_int) :: request
  ! >>>>>>>>>>>>>>>>>>>>>>

  enum, bind(C)
    enumerator :: eDone = 1000, eMatVec, eMatVec_H
  end enum

  !! --------------------------------------------------------------- !!
  !! PThread API.
  !! --------------------------------------------------------------- !!
  type, bind(C) :: ctPThread
    integer(c_long_long) :: hidden
  end type ctPThread
  interface
    subroutine cPThreadCreate(pThread, pAttr, pFunc, env) bind(C, name='pthread_create')
      import :: c_ptr, c_funptr, ctPThread
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(ctPThread), intent(in) :: pThread
      type(c_ptr), intent(in), value :: pAttr
      type(c_funptr), intent(in), value :: pFunc
      type(c_ptr), intent(in), value :: env
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cPThreadCreate
    subroutine cPThreadJoin(pThread, pRetval) bind(C, name='pthread_join')
      import :: c_ptr, ctPThread
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(ctPThread), intent(in), value :: pThread
      type(c_ptr), intent(in), value :: pRetval
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cPThreadJoin
  end interface

  !! --------------------------------------------------------------- !!
  !! Unix semaphores API.
  !! --------------------------------------------------------------- !!
  type, bind(C) :: ctSem
    integer(c_size_t) :: hidden(2)
  end type ctSem
  interface
    subroutine cSemInit(pSem, shared, value) bind(C, name='sem_init')
      import :: c_int, ctSem
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(ctSem), intent(in) :: pSem
      integer(c_int), intent(in), value :: shared, value
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cSemInit
    subroutine cSemClose(pSem) bind(C, name='sem_close')
      import :: ctSem
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(ctSem), intent(in) :: pSem
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cSemClose
    subroutine cSemPost(pSem) bind(C, name='sem_post')
      import :: ctSem
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(ctSem), intent(in) :: pSem
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cSemPost
    subroutine cSemWait(pSem) bind(C, name='sem_wait')
      import :: ctSem
      ! <<<<<<<<<<<<<<<<<<<<<<
      type(ctSem), intent(in) :: pSem
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cSemWait
  end interface
  
  type(c_ptr), save :: spMesh
  type(c_ptr), save :: spX, spB
  integer(c_int), save :: sSolverType, sPrecondType
  type(c_ptr), save :: spAy, spY
  integer(c_int), save :: sRequest

  logical, save :: sFirstCall = .true.
  type(ctPThread), save :: sThread
  type(ctSem), save :: sSemYield, sSemAwait

  if (sFirstCall) then
    ! Save the arguments.
    spMesh = pMesh
    spX = pX; spB = pB
    sSolverType = solverType; sPrecondType = precondType

    ! Create the semaphores and launch the compute thread.
    sFirstCall = .false.
    call cSemInit(sSemYield, 0, 0)
    call cSemInit(sSemAwait, 0, 0)
    call cPThreadCreate(sThread, c_null_ptr, c_funloc(cThreadFunc), c_null_ptr)
  end if

  ! ----------------------
  ! Await for the compute thread to yield.
  ! ----------------------
  call cSemPost(sSemYield)
  call cSemWait(sSemAwait)

  ! ----------------------
  ! Read the request.
  ! ----------------------
  request = sRequest
  pAy = spAy; pY = spY

  ! ----------------------
  ! Clean-up on the exit request.
  ! ----------------------
  if (request == eDone) then
    sFirstCall = .true.
    call cPThreadJoin(sThread, c_null_ptr)
    call cSemClose(sSemYield)
    call cSemClose(sSemAwait)
  end if

contains
  subroutine cThreadFunc(env) bind(C)
    ! <<<<<<<<<<<<<<<<<<<<<<
    type(c_ptr), intent(in), value :: env
    ! >>>>>>>>>>>>>>>>>>>>>>

    print *, 'Entering the RCI thread function..'

    ! ----------------------
    ! Wait for the main thread to resume.
    ! ----------------------
    call cSemWait(sSemYield)

    ! ----------------------
    ! Enter the computational routine.
    ! ----------------------
    call cLinSolve$T(spMesh, spX, spB, &
      & c_funloc(cMatVecFromThread), ceMatVec, &
      & sSolverType, sPrecondType, &
      & c_funloc(cMatVecFromThread), ceMatVec_H)

    ! ----------------------
    ! Pass the exit request to the main thread and leave.
    ! ----------------------
    sRequest = eDone
    call cSemPost(sSemAwait)

    print *, 'Leaving the RCI thread function..'

  end subroutine cThreadFunc
  subroutine cMatVecFromThread(pMesh, pAy, pY, pType) !bind(C)
    ! <<<<<<<<<<<<<<<<<<<<<<
    type(c_ptr), intent(in), value :: pMesh
    type(c_ptr), intent(in), value :: pAy, pY
    type(c_ptr), intent(in), value :: pType
    ! >>>>>>>>>>>>>>>>>>>>>>

    ! ----------------------
    ! Pass the 'MatVec' or 'MatVec_H' request 
    ! and field pointers to the main thread and yield.
    ! ----------------------
    if (c_associated(pType, ceMatVec)) then
      sRequest = eMatVec
    else if (c_associated(pType, ceMatVec_H)) then
      sRequest = eMatVec_H
    end if
    spAy = pAy; spY = pY

    call cSemPost(sSemAwait)
    call cSemWait(sSemYield)

  end subroutine cMatVecFromThread
end function cRCI_LinSolve$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine cApplyBCs$T(pMesh, pU, BCmask, &
      & alpha, beta, gamma) bind(C, name='SR_ApplyBCs$T')
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

end subroutine cApplyBCs$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine cGrad$T(pMesh, pVVec, lambda, pU) bind(C, name='SR_Grad$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  $typename, intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pVVec
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: u(:,:), vVec(:,:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(u, pU); call Unwrap(vVec, pVVec, mesh%NumDims)

  call FDM_Gradient_Central(mesh, vVec, lambda, u)

end subroutine cGrad$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine cDiv$T(pMesh, pV, lambda, pUVec) bind(C, name='SR_Div$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  $typename, intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pUVec, pV
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: uVec(:,:,:), v(:,:)

  call Unwrap(mesh, pMesh)
  call Unwrap(uVec, pUVec, mesh%NumDims); call Unwrap(v, pV)

  call FDM_Divergence_Central(mesh, v, lambda, uVec)

end subroutine cDiv$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine cConv$T(pMesh, pV, lambda, pU, pA) bind(C, name='SR_Conv$T')
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

end subroutine cConv$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine cDivGrad$T(pMesh, pV, lambda, pU) bind(C, name='SR_DivGrad$T')
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(c_ptr), intent(in), value :: pMesh
  $typename, intent(in), value :: lambda
  type(c_ptr), intent(in), value :: pU, pV
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tMesh), pointer :: mesh
  $typename, pointer :: u(:,:), v(:,:)
  $typename, pointer :: uVec(:,:,:), vVec(:,:,:)

  call Unwrap(mesh, pMesh)
  if (cIsVecField$T(pU) == 1) then
    call Unwrap(uVec, pU, mesh%NumDims)
    call Unwrap(vVec, pV, mesh%NumDims)

    call FDM_Laplacian_Central(mesh, vVec, lambda, uVec)
  else
    call Unwrap(u, pU); call Unwrap(v, pV)

    call FDM_Laplacian_Central(mesh, v, lambda, u)
  end if

end subroutine cDivGrad$T
#$end for

!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine cDivKGrad$T(pMesh, pV, lambda, pK, pU) bind(C, name='SR_DivKGrad$T')
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

end subroutine cDivKGrad$T
#$end for

end module StormRuler_API
