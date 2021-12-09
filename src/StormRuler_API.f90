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

use StormRuler_Consts, only: bp, ip, dp

use StormRuler_Helpers, only: PrintBanner, RgbToInt

use StormRuler_Mesh, only: tMesh, InitMeshStencil, InitMeshFromImage
use StormRuler_Mesh_Ordering, only: Mesh_Ordering_Quality, &
  & Mesh_Ordering_Dump, Mesh_Ordering_HilbertCurve, Mesh_Ordering_METIS

use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_IO, only: tIOList => IOList
use StormRuler_IO_VTK!, only: ...

use StormRuler_BLAS, only: Norm_2, &
  & Fill, Fill_Random, Set, Scale, Add, Sub, Mul, &
  & Integrate, FuncProd, SpFuncProd
use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers, only: LinSolve
use StormRuler_Solvers_Newton, only: Solve_JFNK

use StormRuler_FDM_BCs, only: &
  & FDM_ApplyBCs, FDM_ApplyBCs_SlipWall, FDM_ApplyBCs_InOutLet
use StormRuler_FDM_Operators, only: &
  & FDM_Gradient, FDM_Divergence, FDM_DivGrad, FDM_DivWGrad
use StormRuler_FDM_RhieChow, only: FDM_RhieChow_Correction
use StormRuler_FDM_Convection, only: FDM_Convection_Central

use StormRuler_LBM_BCs, only:
use StormRuler_LBM_Operators, only: &
  & LBM_Stream, LBM_Macroscopics, LBM_Collision_BGK

use, intrinsic :: iso_c_binding, only: c_char, c_int, &
  & c_double, c_size_t, c_ptr, c_funptr, c_null_char, &
  & c_associated, c_loc, c_f_pointer, c_f_procpointer

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

#$let KLASSES = ['tMesh', 'tArray', 'tIOList']

!! ----------------------------------------------------------------- !!
!! Class wrapper struct.
!! ----------------------------------------------------------------- !!
#$for klass in KLASSES
type :: tWrapped_$klass
  class($klass), pointer :: mPointer
end type !tWrapped_$klass
#$end for

interface Wrap
#$for klass in KLASSES
  module procedure Wrap_$klass
#$end for
end interface

interface Unwrap
#$for klass in KLASSES
  module procedure Unwrap_$klass
#$end for
  module procedure Unwrap_String
end interface Unwrap

interface Free
#$for klass in KLASSES
  module procedure Free_$klass
#$end for
end interface

abstract interface
  pure subroutine ctMapFunc(size, Fx, x, env) bind(C)
    import :: c_size_t, c_double
    integer(c_size_t), intent(in), value :: size
    real(c_double), intent(in) :: x(*)
    real(c_double), intent(inout) :: Fx(*)
    type(*), intent(in) :: env
  end subroutine ctMapFunc
  pure subroutine ctSpMapFunc(dim, r, size, Fx, x, env) bind(C)
    import :: c_size_t, c_double
    integer(c_size_t), intent(in), value :: dim, size
    real(c_double), intent(in) :: r(*), x(*)
    real(c_double), intent(inout) :: Fx(*)
    type(*), intent(in) :: env
  end subroutine ctSpMapFunc
end interface

abstract interface
  subroutine ctMatVecFunc(meshPtr, AxPtr, xPtr, env) bind(C)
    import :: c_ptr
    type(c_ptr), intent(in), value :: meshPtr
    type(c_ptr), intent(in), value :: AxPtr, xPtr
    type(*), intent(in) :: env
  end subroutine ctMatVecFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Unwrap a class pointer.
!! ----------------------------------------------------------------- !!
#$for klass in KLASSES
function Wrap_$klass(object) result(objectPtr)
  class($klass), intent(in), pointer :: object
  type(c_ptr) :: objectPtr

  type(tWrapped_$klass), pointer :: wrappedObject

  allocate(wrappedObject)
  wrappedObject%mPointer => object
  objectPtr = c_loc(wrappedObject)

end function Wrap_$klass
#$end for

!! ----------------------------------------------------------------- !!
!! Unwrap a class pointer.
!! ----------------------------------------------------------------- !!
#$for klass in KLASSES
subroutine Unwrap_$klass(objectPtr, object, free)
  type(c_ptr), intent(in), value :: objectPtr
  class($klass), intent(out), pointer :: object
  logical, intent(in), optional :: free

  type(tWrapped_$klass), pointer :: wrappedObject

  call c_f_pointer(cptr=objectPtr, fptr=wrappedObject)
  object => wrappedObject%mPointer

  if (present(free)) then
    if (free) deallocate(wrappedObject)
  end if

end subroutine Unwrap_$klass
#$end for

!! ----------------------------------------------------------------- !!
!! Unwrap a string pointer.
!! ----------------------------------------------------------------- !!
subroutine Unwrap_String(stringPtr, string)
  character(c_char), intent(in) :: stringPtr(*)
  character(len=:), intent(out), pointer :: string

  interface
    pure function Strlen(stringPtr) bind(C, name='strlen')
      import :: c_char, c_size_t
      character(c_char), intent(in) :: stringPtr(*)
      integer(c_size_t) :: Strlen
    end function Strlen
    pure subroutine Strncpy(pOutString, stringPtr, len) bind(C, name='strncpy')
      import :: c_char, c_size_t, c_ptr
      type(c_ptr), intent(in), value :: pOutString
      character(c_char), intent(in) :: stringPtr(*)
      integer(c_size_t), intent(in), value :: len
    end subroutine Strncpy
  end interface

  integer(c_size_t) :: len

  len = Strlen(stringPtr)
  allocate(character(len=len) :: string)
  call Strncpy(c_loc(string), stringPtr, len)

end subroutine Unwrap_String

#$for klass in KLASSES
subroutine Free_$klass(objectPtr, mold)
  type(c_ptr), intent(in), value :: objectPtr
  class($klass), intent(in) :: mold

  type(tWrapped_$klass), pointer :: wrappedObject

  call c_f_pointer(cptr=objectPtr, fptr=wrappedObject)
  deallocate(wrappedObject)

end subroutine Free_$klass
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
function cInitMesh() result(meshPtr) bind(C, name='SR_InitMesh')
  use StormRuler_IO, only: Load_PPM
  type(c_ptr) :: meshPtr

  class(tMesh), pointer :: gMesh

  call PrintBanner()

#$if False
  block

    integer(ip), parameter :: nx = 128, ny = 128
    real(dp), parameter :: dx = 1.0_dp/Nx, dy = 1.0_dp/Ny
    allocate(gMesh)
    call gMesh%InitRect(dx, nx, .true., dy, ny, .true., 20)

  end block
#$else
  block

    real(dp), parameter :: dx = 0.01_dp, dy = 0.01_dp
    integer(ip), allocatable :: pixels(:,:)
    integer(ip), allocatable :: colorToMark(:)


    colorToMark = &
      & [ RgbToInt([255, 255, 255]), RgbToInt([255, 0, 0]), &
      &   RgbToInt([  0, 255,   0]), RgbToInt([0, 0, 255]), &
      &   RgbToInt([255,   0, 255]) ]

    allocate(gMesh)

    call Load_PPM('test/Domain-100-Tube.ppm', pixels)
    call InitMeshStencil(gMesh, [Dx,Dy], 'D2Q4')
    call InitMeshFromImage(gMesh, pixels, 0, colorToMark, 2, .true.)

    !call Load_PPM('test/Domain-500-Cube.ppm', pixels)
    !call InitMeshStencil(gMesh, [Dx,Dy], 'D2Q9')
    !call InitMeshFromImage(gMesh, pixels, 0, colorToMark, 1, .false.)

    call gMesh%PrintTo_Neato('test/c2c.dot')

    call IO_WriteVtkImageData(gMesh, 'test/c2c.vts')

  end block
#$endif

#$if False
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
#$end if

  meshPtr = Wrap(gMesh)

end function cInitMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
function cAlloc(meshPtr, numVars, rank) result(yPtr) bind(C, name='SR_Alloc')
  type(c_ptr), intent(in), value :: meshPtr
  integer(c_int), intent(in), value :: numVars, rank
  type(c_ptr) :: yPtr

  integer(ip) :: shape(rank + 2)

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: yArr

  call Unwrap(meshPtr, mesh)

  shape(:rank) = mesh%NumDims
  shape(rank + 1) = numVars
  shape(rank + 2) = mesh%NumAllCells

  allocate(yArr)
  call AllocArray(yArr, shape=shape)

  yPtr = Wrap(yArr)

end function cAlloc

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
function stormAllocOnMesh(meshPtr, rank, &
    & shape) result(yPtr) bind(C, name='stormAllocOnMesh')
  type(c_ptr), intent(in), value :: meshPtr
  integer(c_size_t), intent(in), value :: rank
  integer(c_size_t), intent(in) :: shape(rank)
  type(c_ptr) :: yPtr

  integer(ip) :: trueShape(rank + 1)

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: yArr

  call Unwrap(meshPtr, mesh)

  trueShape(:rank) = shape
  trueShape(rank + 1) = mesh%NumAllCells

  allocate(yArr)
  call AllocArray(yArr, shape=trueShape)

  yPtr = Wrap(yArr)

end function stormAllocOnMesh

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
function stormAllocLike(xPtr) result(yPtr) bind(C, name='stormAllocLike')
  type(c_ptr), intent(in), value :: xPtr
  type(c_ptr) :: yPtr

  class(tArray), pointer :: xArr, yArr

  call Unwrap(xPtr, xArr)

  allocate(yArr)
  call AllocArray(yArr, mold=xArr)

  yPtr = Wrap(yArr)

end function stormAllocLike

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormFree(xPtr) bind(C, name='stormFree')
  type(c_ptr), intent(in), value :: xPtr

  class(tArray), pointer :: xArr

  call Unwrap(xPtr, xArr, free=.true.)
  call FreeArray(xArr)

end subroutine stormFree

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
function cIO_Begin() result(ioListPtr) bind(C, name='SR_IO_Begin') 
  type(c_ptr) :: ioListPtr
  class(tIOList), pointer :: io_List

  allocate(io_List)
  ioListPtr = Wrap(io_List)

end function cIO_Begin

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine cIO_Add(ioListPtr, xPtr, namePtr) bind(C, name='SR_IO_Add') 
  type(c_ptr), intent(in), value :: ioListPtr
  type(c_ptr), intent(in), value :: xPtr
  character(c_char), intent(in) :: namePtr(*)

  class(tIOList), pointer :: io_List
  class(tArray), pointer :: xArr
  real(dp), pointer :: x(:,:)
  character(len=:), pointer :: name

  call Unwrap(ioListPtr, io_List)
  call Unwrap(namePtr, name)
  call Unwrap(xPtr, xArr)

  call xArr%Get(x)

  !! TODO: shape of x?
  if (xArr%mShape(1) == 1) then
    call io_List%Add(name, x(1,:))
  else
    call io_List%Add(name, x)
  end if

end subroutine cIO_Add

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine cIO_Flush(ioListPtr, meshPtr, filenamePtr) bind(C, name='SR_IO_Flush') 
  type(c_ptr), intent(in), value :: ioListPtr
  type(c_ptr), intent(in), value :: meshPtr
  character(c_char), intent(in) :: filenamePtr(*)

  class(tIOList), pointer :: io_List
  class(tMesh), pointer :: mesh
  character(len=:), pointer :: filename

  call Unwrap(ioListPtr, io_List, free=.true.)
  call Unwrap(meshPtr, mesh)
  call Unwrap(filenamePtr, filename)

  call IO_WriteVtkImageData(mesh, filename, io_List)
  deallocate(io_List)

end subroutine cIO_Flush

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormFill(meshPtr, xPtr, alpha) bind(C, name='stormFill')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr
  real(c_double), intent(in), value :: alpha

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr)

  call Fill(mesh, xArr, alpha)

end subroutine stormFill

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormRandFill(meshPtr, xPtr, a, b) bind(C, name='stormRandFill')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr
  real(c_double), intent(in), value :: a, b

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr)

  call Fill_Random(mesh, xArr, a, b)

end subroutine stormRandFill

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormSet(meshPtr, yPtr, xPtr) bind(C, name='stormSet')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr, yPtr

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr); call Unwrap(yPtr, yArr)

  call Set(mesh, yArr, xArr)

end subroutine stormSet

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormScale(meshPtr, yPtr, xPtr, alpha) bind(C, name='stormScale')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr, yPtr
  real(c_double), intent(in), value :: alpha

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr); call Unwrap(yPtr, yArr)

  call Scale(mesh, yArr, xArr, alpha)

end subroutine stormScale

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormAdd(meshPtr, zPtr, yPtr, xPtr, alpha, beta) bind(C, name='stormAdd')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr, yPtr, zPtr
  real(c_double), intent(in), value :: alpha, beta

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr, zArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr); call Unwrap(yPtr, yArr); call Unwrap(zPtr, zArr)

  call Add(mesh, zArr, yArr, xArr, alpha, beta)

end subroutine stormAdd

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormMul(meshPtr, zPtr, yPtr, xPtr) bind(C, name='stormMul')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr, yPtr, zPtr

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr, zArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr); call Unwrap(yPtr, yArr); call Unwrap(zPtr, zArr)

  call Mul(mesh, zArr, yArr, xArr)

end subroutine stormMul

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

#$macro WrapMapFunc ^(?P<cFunc>\w+)\s+(?P<Func>\w+)\s+(?P<env>\w+)\s*$
pure function $Func(x) result(y)
  real(dp), intent(in) :: x(:)
  real(dp) :: y(size(x))

  y(:) = x(:)
  call $cFunc(size(x, kind=c_size_t), y, x, $env)

end function $Func
#$end macro

#$macro WrapSpMapFunc ^(?P<cSpFunc>\w+)\s+(?P<SpFunc>\w+)\s+(?P<env>\w+)\s*$
pure function $SpFunc(r, x) result(y)
  real(dp), intent(in) :: r(:), x(:)
  real(dp) :: y(size(x))

  y(:) = x(:)
  call $cSpFunc(size(r, kind=c_size_t), r, size(x, kind=c_size_t), y, x, $env)

end function $SpFunc
#$end macro

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
real(c_double) function stormIntegrate(meshPtr, xPtr, cFuncPtr, env) &
    & bind(C, name='stormIntegrate')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr
  type(c_funptr), intent(in), value :: cFuncPtr
  type(*), intent(in) :: env

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr

  procedure(ctMapFunc), pointer :: cFunc

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr)

  call c_f_procpointer(cptr=cFuncPtr, fptr=cFunc)

  stormIntegrate = Integrate(mesh, xArr, Func)

contains
  @WrapMapFunc cFunc Func env
end function stormIntegrate

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormFuncProd(meshPtr, yPtr, xPtr, cFuncPtr, env) &
    & bind(C, name='stormFuncProd')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr, yPtr
  type(c_funptr), intent(in), value :: cFuncPtr
  type(*), intent(in) :: env

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr

  procedure(ctMapFunc), pointer :: cFunc

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr); call Unwrap(yPtr, yArr)

  call c_f_procpointer(cptr=cFuncPtr, fptr=cFunc)

  call FuncProd(mesh, yArr, xArr, Func)

contains
  @WrapMapFunc cFunc Func env
end subroutine stormFuncProd

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormSpFuncProd(meshPtr, yPtr, xPtr, cSpFuncPtr, env) &
    & bind(C, name='stormSpFuncProd')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr, yPtr
  type(c_funptr), intent(in), value :: cSpFuncPtr
  type(*), intent(in) :: env

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, yArr
  procedure(ctSpMapFunc), pointer :: cSpFunc

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr); call Unwrap(yPtr, yArr)

  call c_f_procpointer(cptr=cSpFuncPtr, fptr=cSpFunc)

  call SpFuncProd(mesh, yArr, xArr, SpFunc)

contains
  @WrapSpMapFunc cSpFunc SpFunc env
end subroutine stormSpFuncProd

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

#$macro WrapMatVecFunc ^(?P<cMatVec>\w+)\s+(?P<MatVec>\w+)\s+(?P<env>\w+)\s*$
subroutine $MatVec(mesh, yArr, xArr)
  class(tMesh), intent(in), target :: mesh
  class(tArray), intent(inout), target :: xArr, yArr

  type(c_ptr) :: meshPtr, xPtr, yPtr

  meshPtr = Wrap(mesh)
  xPtr = Wrap(xArr); yPtr = Wrap(yArr)

  call $cMatVec(meshPtr, yPtr, xPtr, $env)

  call Free(meshPtr, mold=mesh)
  call Free(xPtr, mold=xArr); call Free(yPtr, mold=yArr)

end subroutine $MatVec
#$end macro

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormLinSolve(meshPtr, methodPtr, preMethodPtr, &
    & xPtr, bPtr, cMatVecPtr, env, cConjMatVecPtr, conjEnv) &
    & bind(C, name='stormLinSolve')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr, bPtr
  character(c_char), intent(in) :: methodPtr(*), preMethodPtr(*)
  type(c_funptr), intent(in), value :: cMatVecPtr, cConjMatVecPtr
  type(*), intent(in) :: env, conjEnv

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, bArr
  character(len=:), pointer :: method, preMethod
  procedure(ctMatVecFunc), pointer :: cMatVec, cConjMatVec
  type(tConvParams) :: params

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr); call Unwrap(bPtr, bArr)
  call Unwrap(methodPtr, method); call Unwrap(preMethodPtr, preMethod)

  call c_f_procpointer(cptr=cMatVecPtr, fptr=cMatVec)
  if (c_associated(cConjMatVecPtr)) then
    call c_f_procpointer(cptr=cConjMatVecPtr, fptr=cConjMatVec)
  end if

  !! TODO: convergence parameters.
  !! TODO: utilize `ConjMatVec` somehow.
  call params%Init(1.0D-6, 1.0D-6, 2000)
  call LinSolve(mesh, method, preMethod, xArr, bArr, MatVec, params)

contains
  @WrapMatVecFunc cMatVec MatVec env
  @WrapMatVecFunc cConjMatVec ConjMatVec conjEnv
end subroutine stormLinSolve

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormNonlinSolve(meshPtr, methodPtr, &
    & xPtr, bPtr, cMatVecPtr, env) bind(C, name='stormNonlinSolve')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: xPtr, bPtr
  character(c_char), intent(in) :: methodPtr(*)
  type(c_funptr), intent(in), value :: cMatVecPtr
  type(*), intent(in) :: env

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: xArr, bArr
  procedure(ctMatVecFunc), pointer :: cMatVec
  type(tConvParams) :: params

  call Unwrap(meshPtr, mesh)
  call Unwrap(xPtr, xArr); call Unwrap(bPtr, bArr)

  call c_f_procpointer(cptr=cMatVecPtr, fptr=cMatVec)

  !! TODO: convergence parameters.
  !! TODO: select method.
  call params%Init(1.0D-4, 1.0D-4, 100, 'JFNK')
  call Solve_JFNK(mesh, MatVec, xArr, bArr, params)

contains
  @WrapMatVecFunc cMatVec MatVec env
end subroutine stormNonlinSolve

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! TODO: some more generic `stormApplyBCs` function
!!       that includes the specific BCs. 

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormApplyBCs(meshPtr, uPtr, markMask, alpha, beta, gamma) &
    & bind(C, name='stormApplyBCs')
  type(c_ptr), intent(in), value :: meshPtr
  integer(c_int), intent(in), value :: markMask
  type(c_ptr), intent(in), value :: uPtr
  real(c_double), intent(in), value :: alpha, beta, gamma

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr
  integer(ip) :: mark, firstMark, lastMark

  call Unwrap(meshPtr, mesh)
  call Unwrap(uPtr, uArr)

  if (markMask /= 0) then
    firstMark = markMask; lastMark = markMask
  else
    firstMark = 1; lastMark = mesh%NumBndMarks
  end if

  do mark = firstMark, lastMark
    block
      !! TODO: fix me!
      real(dp), pointer :: u(:,:)
      call uArr%Get(u)
      call FDM_ApplyBCs(mesh, mark, u, alpha, beta, gamma)
    end block
  end do

end subroutine stormApplyBCs

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormApplyBCs_SlipWall(meshPtr, uPtr, markMask) &
    & bind(C, name='stormApplyBCs_SlipWall')
  type(c_ptr), intent(in), value :: meshPtr
  integer(c_int), intent(in), value :: markMask
  type(c_ptr), intent(in), value :: uPtr

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr
  integer(ip) :: mark, firstMark, lastMark

  call Unwrap(meshPtr, mesh)
  call Unwrap(uPtr, uArr)

  if (markMask /= 0) then
    firstMark = markMask; lastMark = markMask
  else
    firstMark = 1; lastMark = mesh%NumBndMarks
  end if

  do mark = firstMark, lastMark
    block
      !! TODO: fix me!
      real(dp), pointer :: u(:,:)
      call uArr%Get(u)
      call FDM_ApplyBCs_SlipWall(mesh, mark, u)
    end block
  end do

end subroutine stormApplyBCs_SlipWall

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormApplyBCs_InOutLet(meshPtr, uPtr, markMask) &
    & bind(C, name='stormApplyBCs_InOutLet')
  type(c_ptr), intent(in), value :: meshPtr
  integer(c_int), intent(in), value :: markMask
  type(c_ptr), intent(in), value :: uPtr

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr
  integer(ip) :: mark, firstMark, lastMark

  call Unwrap(meshPtr, mesh)
  call Unwrap(uPtr, uArr)

  if (markMask /= 0) then
    firstMark = markMask; lastMark = markMask
  else
    firstMark = 1; lastMark = mesh%NumBndMarks
  end if

  do mark = firstMark, lastMark
    block
      !! TODO: fix me!
      real(dp), pointer :: u(:,:)
      call uArr%Get(u)
      call FDM_ApplyBCs_InOutLet(mesh, mark, u)
    end block
  end do

end subroutine stormApplyBCs_InOutLet

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormGradient(meshPtr, vVecPtr, lambda, uPtr) &
    & bind(C, name='stormGradient')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: uPtr, vVecPtr
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr, vVecArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(uPtr, uArr); call Unwrap(vVecPtr, vVecArr)

  call FDM_Gradient(mesh, vVecArr, lambda, uArr)

end subroutine stormGradient

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormDivergence(meshPtr, vPtr, lambda, uVecPtr) &
    & bind(C, name='stormDivergence')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: uVecPtr, vPtr
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uVecArr, vArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(uVecPtr, uVecArr); call Unwrap(vPtr, vArr)

  call FDM_Divergence(mesh, vArr, lambda, uVecArr)

end subroutine stormDivergence

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormDivGrad(meshPtr, vPtr, lambda, uPtr) &
    & bind(C, name='stormDivGrad')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: uPtr, vPtr
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr, vArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(uPtr, uArr); call Unwrap(vPtr, vArr)

  call FDM_DivGrad(mesh, vArr, lambda, uArr)

end subroutine stormDivGrad

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormDivWGrad(meshPtr, vPtr, lambda, wPtr, uPtr) &
    & bind(C, name='stormDivWGrad')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: uPtr, vPtr, wPtr
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr, vArr, wArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(uPtr, uArr); call Unwrap(vPtr, vArr); call Unwrap(wPtr, wArr)

  call FDM_DivWGrad(mesh, vArr, lambda, wArr, uArr)

end subroutine stormDivWGrad

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormConvection(meshPtr, vPtr, lambda, uPtr, aPtr) &
    & bind(C, name='stormConvection')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value :: uPtr, vPtr, aPtr
  real(c_double), intent(in), value :: lambda

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: uArr, vArr, aArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(uPtr, uArr); call Unwrap(vPtr, vArr); call Unwrap(aPtr, aArr)

  call FDM_Convection_Central(mesh, vArr, lambda, uArr, aArr)

end subroutine stormConvection

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormRhieChowCorrection(meshPtr, vPtr, lambda, tau, &
    & pPtr, rhoPtr) bind(C, name='stormRhieChowCorrection')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value ::  vPtr, pPtr, rhoPtr
  real(c_double), intent(in), value :: lambda, tau

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: vArr, pArr, rhoArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(vPtr, vArr); call Unwrap(pPtr, pArr); call Unwrap(rhoPtr, rhoArr)

  call FDM_RhieChow_Correction(mesh, vArr, lambda, tau, pArr, rhoArr)

end subroutine stormRhieChowCorrection

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormLbmStream(meshPtr, gPtr, fPtr) bind(C, name='stormLbmStream')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value ::  fPtr, gPtr

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: fArr, gArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(fPtr, fArr); call Unwrap(gPtr, gArr)

  call LBM_Stream(mesh, gArr, fArr)

end subroutine stormLbmStream

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormLbmMacroscopics(meshPtr, rhoPtr, vPtr, &
    & sizeOfM, m, fPtr) bind(C, name='stormLbmMacroscopics')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value ::  fPtr, rhoPtr, vPtr
  integer(c_size_t), intent(in), value :: sizeOfM
  real(c_double), intent(in) :: m(sizeOfM)

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: fArr, rhoArr, vArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(fPtr, fArr); call Unwrap(rhoPtr, rhoArr); call Unwrap(vPtr, vArr)

  call LBM_Macroscopics(mesh, rhoArr, vArr, m, fArr)

end subroutine stormLbmMacroscopics

!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
subroutine stormLbmCollisionBGK(meshPtr, fPtr, tau, &
    & rhoPtr, vPtr) bind(C, name='stormLbmCollisionBGK')
  type(c_ptr), intent(in), value :: meshPtr
  type(c_ptr), intent(in), value ::  fPtr, rhoPtr, vPtr
  real(c_double), intent(in), value :: tau

  class(tMesh), pointer :: mesh
  class(tArray), pointer :: fArr, rhoArr, vArr

  call Unwrap(meshPtr, mesh)
  call Unwrap(fPtr, fArr); call Unwrap(rhoPtr, rhoArr); call Unwrap(vPtr, vArr)

  call LBM_Collision_BGK(mesh, fArr, tau, rhoArr, vArr)

end subroutine stormLbmCollisionBGK

end module StormRuler_API
