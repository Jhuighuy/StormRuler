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
use StormRuler_Helpers, only: PrintBanner, PixelToInt

use StormRuler_Mesh, only: tMesh
use StormRuler_Mesh_Ordering, only: Mesh_Ordering_Quality, &
  & Mesh_Ordering_Dump, Mesh_Ordering_HilbertCurve, Mesh_Ordering_METIS

use StormRuler_Array, only: tArray, AllocArray, FreeArray
use StormRuler_IO, only: tIOList => IOList

use StormRuler_BLAS, only: Norm_2, &
  & Fill, Fill_Random, Set, Scale, Add, Sub, Mul, &
  & Integrate, FuncProd, SpFuncProd
use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_Unified, only: LinSolve
use StormRuler_Solvers_Newton, only: Solve_JFNK

use StormRuler_FDM_BCs, only: &
  & FDM_ApplyBCs, FDM_ApplyBCs_SlipWall, FDM_ApplyBCs_InOutLet
use StormRuler_FDM_Operators, only: FDM_Gradient, FDM_Divergence, &
  & FDM_Laplacian_Central, FDM_DivWGrad_Central
use StormRuler_FDM_RhieChow, only: FDM_RhieChow_Correction
use StormRuler_FDM_Convection, only: FDM_Convection_Central

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: c_char, c_int, &
  & c_double, c_size_t, c_ptr, c_funptr, c_null_char, &
  & c_loc, c_f_pointer, c_f_procpointer, c_associated

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

#$let c_typename(T) = {'R': 'real(c_double)'}[T]

#$let KLASSES = ['tMesh', 'tIOList']

!! ----------------------------------------------------------------- !!
!! Class wrapper struct.
!! ----------------------------------------------------------------- !!
#$for klass in KLASSES
type :: c$klass
  class($klass), pointer :: mObject
end type !c$klass
#$end for

!! ----------------------------------------------------------------- !!
!! Field wrapper class.
!! ----------------------------------------------------------------- !!
type :: ctField
  character :: mType = ''
  integer(ip) :: mRank = -1
  real(dp), pointer :: mData(:,:) => null()
  type(tArray), pointer :: mArray => null()
end type ctField

interface Wrap
#$for klass in KLASSES
  module procedure Wrap$klass
#$end for
end interface

interface Unwrap
  module procedure UnwrapString
#$for klass in KLASSES
  module procedure Unwrap$klass
#$end for
  module procedure UnwrapArray
  module procedure UnwrapField
end interface Unwrap

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Unwrap a class pointer.
!! ----------------------------------------------------------------- !!
#$for klass in KLASSES
function Wrap$klass(obj) result(pObj)
  class($klass), intent(inout), pointer :: obj
  type(c_ptr) :: pObj

  type(c$klass), pointer :: cpObj

  allocate(cpObj)
  cpObj%mObject => obj
  pObj = c_loc(cpObj)

end function Wrap$klass
#$end for

!! ----------------------------------------------------------------- !!
!! Unwrap a string pointer.
!! ----------------------------------------------------------------- !!
subroutine UnwrapString(string, pString)
  character(c_char), intent(in) :: pString(*)
  character(len=:), intent(out), pointer :: string

  interface
    pure function cStrlen(pString) bind(C, name='strlen')
      import :: c_char, c_size_t
      character(c_char), intent(in) :: pString(*)
      integer(c_size_t) :: cStrlen
    end function cStrlen
    pure subroutine cStrcpy(pOutString, pString) bind(C, name='strcpy')
      import :: c_char, c_ptr
      type(c_ptr), intent(in), value:: pOutString
      character(c_char), intent(in) :: pString(*)
    end subroutine cStrcpy
  end interface

  allocate(character(len=cStrlen(pString)) :: string)
  call cStrcpy(c_loc(string), pString)

end subroutine UnwrapString

!! ----------------------------------------------------------------- !!
!! Unwrap a class pointer.
!! ----------------------------------------------------------------- !!
#$for klass in KLASSES
subroutine Unwrap$klass(object, pObject, free)
  type(c_ptr), intent(in), value :: pObject
  class($klass), intent(out), pointer :: object
  logical, intent(in), optional :: free

  type(c$klass), pointer :: cObject

  call c_f_pointer(cptr=pObject, fptr=cObject)
  object => cObject%mObject

  if (present(free)) then
    if (free) deallocate(cObject)
  end if

end subroutine Unwrap$klass
#$end for

!! ----------------------------------------------------------------- !!
!! ----------------------------------------------------------------- !!
subroutine UnwrapArray(x, pX)
  type(c_ptr), intent(in), value :: pX
  class(tArray), intent(out), pointer :: x

  type(ctField), pointer :: pX_C

  call c_f_pointer(cptr=pX, fptr=pX_C)
  x => pX_C%mArray

end subroutine UnwrapArray

!! ----------------------------------------------------------------- !!
!! ----------------------------------------------------------------- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine UnwrapField(x, pX, pX_C_out, rank)
  type(c_ptr), intent(in), value :: pX
  real(dp), intent(out), pointer :: x(:,:)
  type(ctField), intent(out), pointer, optional :: pX_C_out
  integer(ip), intent(out), optional :: rank

  type(ctField), pointer :: pX_C

  call c_f_pointer(cptr=pX, fptr=pX_C)
  if (present(pX_C_out)) pX_C_out => pX_C
  x => pX_C%mData
  !x => Reshape2D(pX_C%mData)

  if (present(rank)) then
    rank = pX_C%mRank
  end if

end subroutine UnwrapField
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
function cInitMesh() result(pMesh) bind(C, name='SR_InitMesh')
  use StormRuler_IO, only: Load_PPM
  type(c_ptr) :: pMesh

  class(tMesh), pointer :: gMesh

  call PrintBanner

#$if False
  block

    integer(ip), parameter :: nx = 256, ny = 256
    real(dp), parameter :: dx = 1.0_dp/Nx, dy = 1.0_dp/Ny
    allocate(gMesh)
    call gMesh%InitRect(dx, nx, .true., dy, ny, .true., 20)

  end block
#$else
  block

    real(dp), parameter :: dx = 0.01_dp, dy = 0.01_dp
    integer(ip), allocatable :: pixels(:,:)
    integer(ip), allocatable :: colorToBCM(:)

    allocate(gMesh)

    call Load_PPM('test/Domain-100-Tube.ppm', pixels)

    colorToBCM = &
      & [ PixelToInt([255, 255, 255]), PixelToInt([255, 0, 0]), &
      &   PixelToInt([  0, 255,   0]), PixelToInt([0, 0, 255]), &
      &   PixelToInt([255,   0, 255]) ]
    call gMesh%InitFromImage2D(pixels, 0, colorToBCM, 2)
    
    allocate(gMesh%dr(1:2, 1:4))
    
    gMesh%dl = [Dx,Dx,Dy,Dy]
    gMesh%dr(:,1) = [gMesh%dl(1), 0.0_dp]
    gMesh%dr(:,2) = [gMesh%dl(1), 0.0_dp]
    gMesh%dr(:,3) = [0.0_dp, gMesh%dl(3)]
    gMesh%dr(:,4) = [0.0_dp, gMesh%dl(3)]

    call gMesh%PrintTo_Neato('test/c2c.dot')
    
    call gMesh%PrintTo_LegacyVTK('test/c2c.vtk')

  end block
#$endif

  call gMesh%SetRange()

  block
    integer(ip), allocatable :: iperm(:)

    allocate(iperm(gMesh%NumCells))

    call Mesh_Ordering_Dump(gMesh, 'test/MO-O.txt')
    print *, 'quality(i) = ', Mesh_Ordering_Quality(gMesh)

    call Mesh_Ordering_HilbertCurve(gMesh, iperm)
    call Mesh_Ordering_Dump(gMesh, 'test/MO-H.txt', iperm)
    print *, 'quality(h) = ', Mesh_Ordering_Quality(gMesh, iperm)

#$if HAS_METIS
    call Mesh_Ordering_METIS(gMesh, iperm)
    call Mesh_Ordering_Dump(gMesh, 'test/MO-M.txt', iperm)
    print *, 'quality(m) = ', Mesh_Ordering_Quality(gMesh, iperm)
#$end if

    call gMesh%ApplyOrdering(iperm)

  end block

  pMesh = Wrap(gMesh)

end function cInitMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
#$for T, typename in [SCALAR_TYPES[0]]
function cAlloc(pMesh, numVars, rank) result(pY) bind(C, name='SR_Alloc')
  type(c_ptr), intent(in), value :: pMesh
  integer(c_int), intent(in), value :: numVars, rank
  type(c_ptr) :: pY

  integer(ip) :: shape(rank + 2)

  class(tMesh), pointer :: mesh
  type(ctField), pointer :: pY_C

  call Unwrap(mesh, pMesh)

  shape(:rank) = mesh%NumDims
  shape(rank + 1) = numVars
  shape(rank + 2) = mesh%NumAllCells

  allocate(pY_C)
#$if T == 'R'
  allocate(pY_C%mArray)
  call AllocArray(pY_C%mArray, shape=shape)
  call pY_C%mArray%Get(pY_C%mData)
#$end if
  pY_C%mRank = rank
  pY = c_loc(pY_C)

end function cAlloc
#$end for

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
#$for T, typename in [SCALAR_TYPES[0]]
function stormAllocLike(pX) result(pY) bind(C, name='stormAllocLike')
  type(c_ptr), intent(in), value :: pX
  type(c_ptr) :: pY

  type(ctField), pointer :: pX_C, pY_C
  real(dp), pointer :: x(:,:)

  call Unwrap(x, pX, pX_C_out=pX_C)
  allocate(pY_C)
#$if T == 'R'
  allocate(pY_C%mArray)
  call AllocArray(pY_C%mArray, mold=pX_C%mArray)
  call pY_C%mArray%Get(pY_C%mData)
#$end if
  pY_C%mRank = pX_C%mRank
  pY = c_loc(pY_C)

end function stormAllocLike
#$end for

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormFree(pX) bind(C, name='stormFree')
  type(c_ptr), intent(in), value :: pX

  type(ctField), pointer :: pX_C
  real(dp), pointer :: x(:,:)

  call Unwrap(x, pX, pX_C_out=pX_C)
  call pX_C%mArray%Free()
  deallocate(pX_C%mArray)

end subroutine stormFree

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
function cIO_Begin() result(pIOList) bind(C, name='SR_IO_Begin') 
  type(c_ptr) :: pIOList
  class(tIOList), pointer :: ioList

  allocate(ioList)
  pIOList = Wrap(ioList)

end function cIO_Begin

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine cIO_Add(pIOList, pX, pName) bind(C, name='SR_IO_Add') 
  type(c_ptr), intent(in), value :: pIOList
  type(c_ptr), intent(in), value :: pX
  character(c_char), intent(in) :: pName(*)

  class(tIOList), pointer :: ioList
  type(ctField), pointer :: pX_C
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

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine cIO_Flush(pIOList, pMesh, pFileName) bind(C, name='SR_IO_Flush') 
  type(c_ptr), intent(in), value :: pIOList
  type(c_ptr), intent(in), value :: pMesh
  character(c_char), intent(in) :: pFileName(*)

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

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormFill(pMesh, pX, alpha, beta) bind(C, name='stormFill')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX
  real(c_double), intent(in), value :: alpha
  real(c_double), intent(in), value :: beta

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX)

  call Fill(mesh, xArr, alpha, beta)

end subroutine stormFill

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormRandFill(pMesh, pX, a, b) bind(C, name='stormRandFill')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX
  real(c_double), intent(in), value :: a, b

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX)

  call Fill_Random(mesh, xArr, a, b)

end subroutine stormRandFill

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormSet(pMesh, pY, pX) bind(C, name='stormSet')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX); call Unwrap(yArr, pY)

  call Set(mesh, yArr, xArr)

end subroutine stormSet

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormScale(pMesh, pY, pX, alpha) bind(C, name='stormScale')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  real(c_double), intent(in), value :: alpha

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX); call Unwrap(yArr, pY)

  call Scale(mesh, yArr, xArr, alpha)

end subroutine stormScale

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormAdd(pMesh, pZ, pY, pX, alpha, beta) bind(C, name='stormAdd')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY, pZ
  real(c_double), intent(in), value :: alpha, beta

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr, zArr

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX); call Unwrap(yArr, pY); call Unwrap(zArr, pZ)

  call Add(mesh, zArr, yArr, xArr, alpha, beta)

end subroutine stormAdd

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormMul(pMesh, pZ, pY, pX) bind(C, name='stormMul')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY, pZ

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr, zArr

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX); call Unwrap(yArr, pY); call Unwrap(zArr, pZ)

  call Mul(mesh, zArr, yArr, xArr)

end subroutine stormMul

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
real(c_double) function stormIntegrate(pMesh, pX, pF, pEnv) &
    & bind(C, name='stormIntegrate')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX
  type(c_funptr), intent(in), value :: pF
  type(c_ptr), intent(in), value :: pEnv

  abstract interface
    pure subroutine ctMapFunc(size, Fx, x, pEnv) bind(C)
      import :: c_int, c_ptr, dp
      integer(c_int), intent(in), value :: size
      real(dp), intent(in) :: x(*)
      real(dp), intent(inout) :: Fx(*)
      type(c_ptr), intent(in), value :: pEnv
    end subroutine ctMapFunc
  end interface

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr
  procedure(ctMapFunc), pointer :: f

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX)
  call c_f_procpointer(cptr=pF, fptr=f)

  stormIntegrate = Integrate(mesh, xArr, cFunc)

contains
  pure function cFunc(x) result(Fx)
    real(dp), intent(in) :: x(:)
    real(dp) :: Fx(size(x))
    
    call f(int(size(x), kind=c_int), Fx, x, pEnv)

  end function cFunc
end function stormIntegrate

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormFuncProd(pMesh, pY, pX, pF, pEnv) bind(C, name='stormFuncProd')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  type(c_funptr), intent(in), value :: pF
  type(*), intent(in) :: pEnv

  abstract interface
    pure subroutine ctMapFunc(size, Fx, x, pEnv) bind(C)
      import :: c_int, c_double
      integer(c_int), intent(in), value :: size
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(inout) :: Fx(*)
      type(*), intent(in) :: pEnv
    end subroutine ctMapFunc
  end interface

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr
  procedure(ctMapFunc), pointer :: f

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX); call Unwrap(yArr, pY)
  call c_f_procpointer(cptr=pF, fptr=f)

  call FuncProd(mesh, yArr, xArr, cFunc)

contains
  pure function cFunc(x) result(Fx)
    real(c_double), intent(in) :: x(:)
    real(c_double) :: Fx(size(x))
    
    call f(int(size(x), kind=c_int), Fx, x, pEnv)

  end function cFunc
end subroutine stormFuncProd

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormSpFuncProd(pMesh, pY, pX, pF, pEnv) bind(C, name='stormSpFuncProd')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pY
  type(c_funptr), intent(in), value :: pF
  type(*), intent(in) :: pEnv

  abstract interface
    pure subroutine ctSpMapFunc(dim, r, size, Fx, x, pEnv) bind(C)
      import :: c_int, c_double
      integer(c_int), intent(in), value :: dim, size
      real(c_double), intent(in) :: r(*)
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(inout) :: Fx(*)
      type(*), intent(in) :: pEnv
    end subroutine ctSpMapFunc
  end interface

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr
  procedure(ctSpMapFunc), pointer :: f

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX); call Unwrap(yArr, pY)
  call c_f_procpointer(cptr=pF, fptr=f)

  call SpFuncProd(mesh, yArr, xArr, cFunc)

contains
  pure function cFunc(r, x) result(Fx)
    real(dp), intent(in) :: r(:)
    real(c_double), intent(in) :: x(:)
    real(c_double) :: Fx(size(x))
    
    call f(int(size(r), kind=c_int), r, int(size(x), kind=c_int), Fx, x, pEnv)

  end function cFunc
end subroutine stormSpFuncProd

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormLinSolve(pMesh, pMethod, pPreMethod, &
    & pX, pB, pMatVec, pEnv, pMatVec_H, pEnv_H) bind(C, name='stormLinSolve')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pB
  character(c_char), intent(in) :: pMethod(*), pPreMethod(*)
  type(c_funptr), intent(in), value :: pMatVec, pMatVec_H
  type(*), intent(in) :: pEnv, pEnv_H

  abstract interface
    subroutine ctMatVec(pMesh, pAx, pX, pEnv) bind(C)
      import :: c_ptr, dp
      type(c_ptr), intent(in), value :: pMesh
      type(c_ptr), intent(in), value :: pAx, pX
      type(*), intent(in) :: pEnv
    end subroutine ctMatVec
  end interface

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, bArr
  procedure(ctMatVec), pointer :: MatVec, MatVec_H
  character(len=:), pointer :: method, preMethod
  type(tConvParams) :: params

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX); call Unwrap(bArr, pB)
  call c_f_procpointer(cptr=pMatVec, fptr=MatVec)
  if (c_associated(pMatVec_H)) then
    call c_f_procpointer(cptr=pMatVec_H, fptr=MatVec_H)
  end if
  call Unwrap(method, pMethod); call Unwrap(preMethod, pPreMethod)

  call params%Init(1.0D-6, 1.0D-6, 2000)
  call LinSolve(mesh, method, preMethod, xArr, bArr, cMatVec, params)

contains
  subroutine cMatVec(mesh_, Ax, x)
    class(tMesh), intent(inout), target :: mesh_
    class(tArray), intent(inout), target :: x, Ax

    type(ctField), target :: pX_C, pAx_C
    type(c_ptr) :: pX, pAx

    pX_C%mArray => x
    pAx_C%mArray => Ax
    call pX_C%mArray%Get(pX_C%mData)
    call pAx_C%mArray%Get(pAx_C%mData)
    pX = c_loc(pX_C); pAx = c_loc(pAx_C)

    call MatVec(pMesh, pAx, pX, pEnv)

  end subroutine cMatVec
end subroutine stormLinSolve

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormNonlinSolve(pMesh, pMethod, &
    & pX, pB, pMatVec, pEnv) bind(C, name='stormNonlinSolve')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pX, pB
  character(c_char), intent(in) :: pMethod(*)
  type(c_funptr), intent(in), value :: pMatVec
  type(*), intent(in) :: pEnv

  abstract interface
    subroutine ctMatVec(pMesh, pAx, pX, pEnv) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: pMesh
      type(c_ptr), intent(in), value :: pAx, pX
      type(*), intent(in) :: pEnv
    end subroutine ctMatVec
  end interface

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, bArr
  procedure(ctMatVec), pointer :: MatVec, MatVec_H
  type(tConvParams) :: params

  call Unwrap(mesh, pMesh)
  call Unwrap(xArr, pX); call Unwrap(bArr, pB)
  call c_f_procpointer(cptr=pMatVec, fptr=MatVec)

  call params%Init(1.0D-4, 1.0D-4, 100, 'JFNK')
  call Solve_JFNK(mesh, cMatVec, xArr, bArr, params)

contains
  subroutine cMatVec(mesh_, Ax, x)
    class(tMesh), intent(inout), target :: mesh_
    class(tArray), intent(inout), target :: x, Ax

    type(ctField), target :: pX_C, pAx_C
    type(c_ptr) :: pX, pAx

    pX_C%mArray => x
    pAx_C%mArray => Ax
    call pX_C%mArray%Get(pX_C%mData)
    call pAx_C%mArray%Get(pAx_C%mData)
    pX = c_loc(pX_C); pAx = c_loc(pAx_C)

    call MatVec(pMesh, pAx, pX, pEnv)

  end subroutine cMatVec
end subroutine stormNonlinSolve

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormApplyBCs(pMesh, pU, BCmask, alpha, beta, gamma) bind(C, name='stormApplyBCs')
  type(c_ptr), intent(in), value :: pMesh
  integer(c_int), intent(in), value :: BCmask
  type(c_ptr), intent(in), value :: pU
  real(dp), intent(in), value :: alpha, beta, gamma

  class(tMesh), pointer :: mesh
  real(dp), pointer :: u(:,:)
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

end subroutine stormApplyBCs

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormApplyBCs_SlipWall(pMesh, pU, BCmask) bind(C, name='stormApplyBCs_SlipWall')
  type(c_ptr), intent(in), value :: pMesh
  integer(c_int), intent(in), value :: BCmask
  type(c_ptr), intent(in), value :: pU

  class(tMesh), pointer :: mesh
  real(dp), pointer :: u(:,:)
  integer(ip) :: iBC, firstBC, lastBC

  call Unwrap(mesh, pMesh)
  call Unwrap(u, pU)

  if (BCmask /= 0) then
    firstBC = BCmask; lastBC = BCmask
  else
    firstBC = 1; lastBC = mesh%NumBCMs
  end if

  do iBC = firstBC, lastBC
    call FDM_ApplyBCs_SlipWall(mesh, iBC, u)
  end do

end subroutine stormApplyBCs_SlipWall

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormApplyBCs_InOutLet(pMesh, pU, BCmask) bind(C, name='stormApplyBCs_InOutLet')
  type(c_ptr), intent(in), value :: pMesh
  integer(c_int), intent(in), value :: BCmask
  type(c_ptr), intent(in), value :: pU

  class(tMesh), pointer :: mesh
  real(dp), pointer :: u(:,:)
  integer(ip) :: iBC, firstBC, lastBC

  call Unwrap(mesh, pMesh)
  call Unwrap(u, pU)

  if (BCmask /= 0) then
    firstBC = BCmask; lastBC = BCmask
  else
    firstBC = 1; lastBC = mesh%NumBCMs
  end if

  do iBC = firstBC, lastBC
    call FDM_ApplyBCs_InOutLet(mesh, iBC, u)
  end do

end subroutine stormApplyBCs_InOutLet

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormGradient(pMesh, pVVec, lambda, pU) bind(C, name='stormGradient')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pU, pVVec
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr, vVecArr

  call Unwrap(mesh, pMesh)
  call Unwrap(uArr, pU); call Unwrap(vVecArr, pVVec)

  call FDM_Gradient(mesh, vVecArr, lambda, uArr)

end subroutine stormGradient

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormDivergence(pMesh, pV, lambda, pUVec) bind(C, name='stormDivergence')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pUVec, pV
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uVecArr, vArr

  call Unwrap(mesh, pMesh)
  call Unwrap(uVecArr, pUVec); call Unwrap(vArr, pV)

  call FDM_Divergence(mesh, vArr, lambda, uVecArr)

end subroutine stormDivergence

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormRhieChowCorrection(pMesh, pV, lambda, tau, &
    & pP, pRho) bind(C, name='stormRhieChowCorrection')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value ::  pV, pP, pRho
  real(c_double), intent(in), value :: lambda, tau

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uVecArr, vArr, pArr, rhoArr

  call Unwrap(mesh, pMesh)
  call Unwrap(vArr, pV); call Unwrap(pArr, pP); call Unwrap(rhoArr, pRho)

  call FDM_RhieChow_Correction(mesh, vArr, lambda, tau, pArr, rhoArr)

end subroutine stormRhieChowCorrection

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormConvection(pMesh, pV, lambda, pU, pA) bind(C, name='stormConvection')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pU, pV, pA
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr, vArr, aArr

  call Unwrap(mesh, pMesh)
  call Unwrap(uArr, pU); call Unwrap(vArr, pV); call Unwrap(aArr, pA)

  call FDM_Convection_Central(mesh, vArr, lambda, uArr, aArr)

end subroutine stormConvection

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormDivGrad(pMesh, pV, lambda, pU) bind(C, name='stormDivGrad')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pU, pV
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr, vArr

  call Unwrap(mesh, pMesh)
  call Unwrap(uArr, pU); call Unwrap(vArr, pV)

  call FDM_Laplacian_Central(mesh, vArr, lambda, uArr)

end subroutine stormDivGrad

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormDivWGrad(pMesh, pV, lambda, pW, pU) bind(C, name='stormDivWGrad')
  type(c_ptr), intent(in), value :: pMesh
  type(c_ptr), intent(in), value :: pU, pV, pW
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr, vArr, wArr

  call Unwrap(mesh, pMesh)
  call Unwrap(uArr, pU); call Unwrap(vArr, pV); call Unwrap(wArr, pW)

  call FDM_DivWGrad_Central(mesh, vArr, lambda, wArr, uArr)

end subroutine stormDivWGrad

end module StormRuler_API
