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
module StormRuler_Preconditioner_ILU

use StormRuler_Consts, only: ip, dp
use StormRuler_Parameters, only: gUseMKL, gDiagAdjustValueILU0, &
  & gToleranceILUT, gDiagAdjustValueILUT, gMaxFillILUT

use StormRuler_Helpers, only: ErrorStop, I2S

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_Matrix, only: tMatrix
use StormRuler_Matrix, only: tMatrix, tParallelTriangularContext, &
  & InitParallelTriangularContext, SolveTriangular
use StormRuler_Preconditioner, only: tMatrixBasedPreconditioner

#$use 'StormRuler_Macros.fi'
#$if HAS_MKL
use StormRuler_Libs_MKL, only: dcsrilu0, dcsrilut, mkl_dcsrtrsv
#$end if

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Base class for the incomplete LU preconditioners.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tMatrixBasedPreconditioner), abstract :: tBaseIluPreconditioner
  type(tMatrix), pointer, private :: Mat => null()
  type(tMatrix), private :: LuMat
  type(tParallelTriangularContext), private :: LowerLuCtx, UpperLuCtx

contains
  procedure(tMakeIluMatrixFunc), deferred :: MakeIluMatrix
  procedure, non_overridable :: SetMatrix => SetIluPreconditionerMatrix
  procedure, non_overridable :: Init => InitIluPreconditioner
  procedure, non_overridable :: Apply => ApplyIluPreconditioner

end type tBaseIluPreconditioner

abstract interface
  subroutine tMakeIluMatrixFunc(pre, mesh, luMat, mat)
    import :: tBaseIluPreconditioner, tMesh, tMatrix
    class(tBaseIluPreconditioner), intent(inout) :: pre
    class(tMesh), intent(in) :: mesh
    class(tMatrix), intent(inout) :: luMat
    class(tMatrix), intent(in) :: mat
  end subroutine tMakeIluMatrixFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! ILU(0) preconditioner.
!!
!! ILU(0) preconditioner does not preserve symmetry, so it must
!! be used with the general solvers like BiCGStab, GMRES or TFQMR. 
!!
!! Like the other triangular matrix-based preconditioners, ILU(0)
!! may suffer from poor parallel scaling.
!! Initialization of the non-block ILU(0) precondtioner may be
!! accelerated with Intel MKL.
!!
!! References:
!! [1] ???
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tBaseIluPreconditioner) :: tIlu0Preconditioner
contains
  procedure, non_overridable :: MakeIluMatrix => MakeIlu0Matrix

end type tIlu0Preconditioner

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! ILUT preconditioner.
!!
!! ILUT preconditioner does not preserve symmetry, so it must
!! be used with the general solvers like BiCGStab, GMRES or TFQMR. 
!!
!! Like the other triangular matrix-based preconditioners, ILU(0)
!! may suffer from poor parallel scaling.
!! Initialization of the ILUT precondtioner may be
!! accelerated with Intel MKL.
!!
!! References:
!! [1] ???
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tBaseIluPreconditioner) :: tIlutPreconditioner
contains
  procedure, non_overridable :: MakeIluMatrix => MakeIlutMatrix

end type tIlutPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set the ILU preconditioner matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SetIluPreconditionerMatrix(pre, mat)
  class(tBaseIluPreconditioner), intent(inout) :: pre
  class(tMatrix), intent(inout), target :: mat

  pre%Mat => mat

end subroutine SetIluPreconditionerMatrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize the ILU preconditioner: 𝓟 ← 𝘪𝘯𝘪𝘵(𝓐).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitIluPreconditioner(pre, mesh, MatVec)
  class(tBaseIluPreconditioner), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  procedure(tMatVecFunc) :: MatVec

  if (.not.associated(pre%Mat)) then
    call ErrorStop('Matrix for the ILU preconditioner is not set.')
  end if

  call pre%MakeIluMatrix(mesh, pre%LuMat, pre%Mat)

end subroutine InitIluPreconditioner

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the ILU preconditioner: 𝒚 ← 𝓟(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ApplyIluPreconditioner(pre, mesh, yArr, xArr, MatVec)
  class(tBaseIluPreconditioner), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  class(tArray), intent(inout), target :: xArr, yArr
  procedure(tMatVecFunc) :: MatVec

  type(tArray) :: tArr

  if (.not.associated(pre%Mat)) then
    call ErrorStop('Matrix for the ILU preconditioner is not set.')
  end if

  call AllocArray(tArr, mold=xArr)

  ! ----------------------
  ! Apply the ILU preconditioner:
  ! 𝒕 ← (𝓛 + 𝓘)⁻¹𝒙,
  ! 𝒚 ← (𝓓 + 𝓤)⁻¹𝒕.
  ! ----------------------
  !! TODO: reimplement with built-in triangular solvers.
#$if HAS_MKL
  block
    real(dp), pointer :: t(:), x(:)
    call tArr%Get(t); call xArr%Get(x)
    call mkl_dcsrtrsv('L', 'N', 'U', mesh%NumCells, &
      & pre%LuMat%ColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, x, t)
  end block
#$else
  error stop 'SolveTriangular with unit diagonal is not implemented yet.'
#$end if
  call SolveTriangular(mesh, 'U', pre%LuMat, pre%UpperLuCtx, yArr, tArr)

end subroutine ApplyIluPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

#$if HAS_MKL
!! ----------------------------------------------------------------- !!
!! Initizize the MKL's ILU preconditioner parameters.
!! ----------------------------------------------------------------- !!
subroutine InitMklIluParameters(iparam, dparam)
  integer(ip), intent(inout) :: iparam(128)
  real(dp), intent(inout) :: dparam(128)

  iparam(:) = 0; dparam(:) = 0.0_dp
  iparam(2) = 6 ! output messages on screen.

end subroutine InitMklIluParameters
subroutine InitMklIlu0Parameters(iparam, dparam)
  integer(ip), intent(inout) :: iparam(128)
  real(dp), intent(inout) :: dparam(128)

  call InitMklIluParameters(iparam, dparam)
  dparam(31) = gDiagAdjustValueILU0 ! diagonal threshold value.
  dparam(32) = gDiagAdjustValueILU0 ! diagonal adjustment value.

end subroutine InitMklIlu0Parameters
subroutine InitMklIlutParameters(iparam, dparam)
  integer(ip), intent(inout) :: iparam(128)
  real(dp), intent(inout) :: dparam(128)

  call InitMklIluParameters(iparam, dparam)
  iparam(7) = 1 ! output the error message and continue.
  dparam(31) = gDiagAdjustValueILUT ! diagonal threshold and adjustment value.

end subroutine InitMklIlutParameters
#$end if

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Make the ILU(0) preconditioner matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine MakeIlu0Matrix(pre, mesh, luMat, mat)
  class(tIlu0Preconditioner), intent(inout) :: pre
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(inout) :: luMat
  class(tMatrix), intent(in) :: mat

  integer(ip) :: size

  size = mat%BlockSize()
  luMat%RowAddrs => mat%RowAddrs
  luMat%ColIndices => mat%ColIndices
  allocate(luMat%ColCoeffs, mold=mat%ColCoeffs)

#$if HAS_MKL
  if (gUseMKL.and.(size == 1)) then
    block
      integer(ip) :: iparam(128), ierror
      real(dp) :: dparam(128)

      call InitMklIlu0Parameters(iparam, dparam)
      call dcsrilu0(mesh%NumCells, mat%ColCoeffs, mat%RowAddrs, &
        & mat%ColIndices, luMat%ColCoeffs, iparam, dparam, ierror)
      if (ierror /= 0) then
        call ErrorStop('MKL `dcsrilu0` has failed, ierror='//I2S(ierror))
      end if

    end block
    return
  end if
#$end if

  !! TODO: implement the ILU(0) preconditioner without MKL.
  call ErrorStop('ILU(0) preconditioner requires MKL for now.')

end subroutine MakeIlu0Matrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Make the ILUT preconditioner matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine MakeIlutMatrix(pre, mesh, luMat, mat)
  class(tIlutPreconditioner), intent(inout) :: pre
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(inout) :: luMat
  class(tMatrix), intent(in) :: mat

  integer(ip) :: size

  size = mat%BlockSize()
  allocate(luMat%RowAddrs, mold=mat%RowAddrs)
  associate(luNnz => (2*gMaxFillILUT + 1)*mesh%NumCells - &
                   & gMaxFillILUT*(gMaxFillILUT + 1) + 1)
    allocate(luMat%ColIndices(luNnz))
    allocate(luMat%ColCoeffs(size, size, luNnz))
  end associate

#$if HAS_MKL
  if (gUseMKL.and.(size == 1)) then
    block
      integer(ip) :: iparam(128), ierror
      real(dp) :: dparam(128)

      call InitMklIlutParameters(iparam, dparam)
      call dcsrilut(mesh%NumCells, &
        & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, &
        & luMat%ColCoeffs, luMat%RowAddrs, luMat%ColIndices, &
        & gToleranceILUT, gMaxFillILUT, iparam, dparam, ierror)
      if (ierror /= 0) then
        call ErrorStop('MKL `dcsrilut` has failed, ierror='//I2S(ierror))
      end if

    end block
    return
  end if
#$end if

  !! TODO: implement the ILUT preconditioner without MKL.
  call ErrorStop('ILUT preconditioner requires MKL for now.')

end subroutine MakeIlutMatrix

end module StormRuler_Preconditioner_ILU
