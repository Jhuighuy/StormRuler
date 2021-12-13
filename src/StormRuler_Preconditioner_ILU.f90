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
use StormRuler_Preconditioner, only: tMatrixPreconditioner

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
type, extends(tMatrixPreconditioner), abstract :: tBaseIluPreconditioner
  type(tMatrix), pointer, private :: Mat => null()
  type(tMatrix), private :: IluMat
  type(tParallelTriangularContext), private :: LowerIluCtx, UpperIluCtx

contains
  procedure(tMakeIluMatrixFunc), deferred :: MakeIluMatrix
  procedure, non_overridable :: SetMatrix => SetIluPreconditionerMatrix
  procedure, non_overridable :: Init => InitIluPreconditioner
  procedure, non_overridable :: Apply => ApplyIluPreconditioner

end type tBaseIluPreconditioner

abstract interface
  subroutine tMakeIluMatrixFunc(pre, mesh, iluMat, mat)
    import :: tBaseIluPreconditioner, tMesh, tMatrix
    class(tBaseIluPreconditioner), intent(inout) :: pre
    class(tMesh), intent(in) :: mesh
    class(tMatrix), intent(inout) :: iluMat
    class(tMatrix), intent(in) :: mat
  end subroutine tMakeIluMatrixFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Incomplete LU-0 (ILU0) preconditioner.
!!
!! ILU0 preconditioner does not preserve symmetry, so it must
!! be used with the general solvers like BiCGStab, GMRES or TFQMR. 
!!
!! Like the other triangular matrix-based preconditioners, ILU0
!! may suffer from poor parallel scaling.
!! Initialization of the point ILU0 precondtioner may be
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
!! Incomplete LU-threshold (ILUT) preconditioner.
!!
!! ILUT preconditioner does not preserve symmetry, so it must
!! be used with the general solvers like BiCGStab, GMRES or TFQMR. 
!!
!! Like the other triangular matrix-based preconditioners, ILU(0)
!! may suffer from poor parallel scaling.
!! Initialization of the point ILUT precondtioner may be
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

  call pre%MakeIluMatrix(mesh, pre%IluMat, pre%Mat)

  !! TODO: refactor with an initialization function.
  if (.not.allocated(pre%LowerIluCtx%LevelAddrs)) &
    & call InitParallelTriangularContext(mesh, 'LD', pre%mat, pre%LowerIluCtx)
  if (.not.allocated(pre%UpperIluCtx%LevelAddrs)) &
    & call InitParallelTriangularContext(mesh, 'IU', pre%mat, pre%UpperIluCtx)

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
  call SolveTriangular(mesh, 'IL', pre%IluMat, pre%LowerIluCtx, tArr, xArr)
  call SolveTriangular(mesh, 'DU', pre%IluMat, pre%UpperIluCtx, yArr, tArr)

end subroutine ApplyIluPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Make the ILU0 preconditioner matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine MakeIlu0Matrix(pre, mesh, iluMat, mat)
  class(tIlu0Preconditioner), intent(inout) :: pre
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(inout) :: iluMat
  class(tMatrix), intent(in) :: mat

  integer(ip) :: size

  size = mat%BlockSize()
  iluMat%RowAddrs => mat%RowAddrs
  iluMat%ColIndices => mat%ColIndices
  allocate(iluMat%ColCoeffs, mold=mat%ColCoeffs)

#$if HAS_MKL
  if (gUseMKL.and.(size == 1)) then
    block
      integer(ip) :: iparam(128), errorCode
      real(dp) :: dparam(128)

      iparam(:) = 0; dparam(:) = 0.0_dp
      iparam(2) = 6 ! output messages on screen.
      dparam(31) = gDiagAdjustValueILU0 ! diagonal threshold value.
      dparam(32) = gDiagAdjustValueILU0 ! diagonal adjustment value.

      call dcsrilu0(mesh%NumCells, &
        & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, & 
        & iluMat%ColCoeffs, iparam, dparam, errorCode)
      if (errorCode /= 0) then
        call ErrorStop('MKL `dcsrilu0` has failed, errorCode='//I2S(errorCode))
      end if

    end block
    return
  end if
#$end if

  !! TODO: implement the ILU0 preconditioner without MKL.
  call ErrorStop('ILU0 preconditioner requires MKL for now.')

end subroutine MakeIlu0Matrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Make the ILUT preconditioner matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine MakeIlutMatrix(pre, mesh, iluMat, mat)
  class(tIlutPreconditioner), intent(inout) :: pre
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(inout) :: iluMat
  class(tMatrix), intent(in) :: mat

  integer(ip) :: size

  size = mat%BlockSize()
  allocate(iluMat%RowAddrs, mold=mat%RowAddrs)
  associate(luNnz => (2*gMaxFillILUT + 1)*mesh%NumCells - &
                      & gMaxFillILUT*(gMaxFillILUT + 1) + 1)
    allocate(iluMat%ColIndices(luNnz))
    allocate(iluMat%ColCoeffs(size, size, luNnz))
  end associate

#$if HAS_MKL
  if (gUseMKL.and.(size == 1)) then
    block
      integer(ip) :: iparam(128), errorCode
      real(dp) :: dparam(128)

      iparam(:) = 0; dparam(:) = 0.0_dp
      iparam(2) = 6 ! output messages on screen.
      iparam(7) = 1 ! allow output of errors.
      iparam(31) = 1 ! adjust the small diagonal entries.
      dparam(31) = gDiagAdjustValueILUT ! diagonal threshold and adjustment value.

      call dcsrilut(mesh%NumCells, &
        & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, &
        & iluMat%ColCoeffs, iluMat%RowAddrs, iluMat%ColIndices, &
        & gToleranceILUT, gMaxFillILUT, iparam, dparam, errorCode)
      if (errorCode /= 0) then
        call ErrorStop('MKL `dcsrilut` has failed, errorCode='//I2S(errorCode))
      end if

    end block
    return
  end if
#$end if

  !! TODO: implement the ILUT preconditioner without MKL.
  call ErrorStop('ILUT preconditioner requires MKL for now.')

end subroutine MakeIlutMatrix

end module StormRuler_Preconditioner_ILU
