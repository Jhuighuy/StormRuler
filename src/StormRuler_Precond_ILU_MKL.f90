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
module StormRuler_Precond_ILU_MKL

#$use 'StormRuler_Params.fi'
#$if HAS_MKL

use StormRuler_Parameters, only: dp, ip
use StormRuler_Helpers, only: I2S

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_Precond, only: tMatrixBasedPreconditioner
use StormRuler_Matrix, only: tMatrix

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

include 'mkl_spblas.fi'
include 'mkl_rci.fi'

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Incomplete LU preconditioner, MKL based.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tMatrixBasedPreconditioner) :: tPreconditioner_ILU0_MKL
  type(tMatrix), pointer, private :: Mat => null()
  real(dp), allocatable, private :: FactorsColCoeffs(:)

contains
  procedure SetMatrix => SetPrecondMatrix_ILU0_MKL
  procedure Init => InitPrecond_ILU0_MKL
  procedure Apply => ApplyPrecond_ILU0_MKL
end type tPreconditioner_ILU0_MKL

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

subroutine SetPrecondMatrix_ILU0_MKL(pre, mat)
  class(tPreconditioner_ILU0_MKL), intent(inout) :: pre
  class(tMatrix), intent(inout), target :: mat

  pre%Mat => mat

end subroutine SetPrecondMatrix_ILU0_MKL

subroutine InitPrecond_ILU0_MKL(pre, mesh, MatVec)
  class(tPreconditioner_ILU0_MKL), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: iparam(128), ierror
  real(dp) :: dparam(128)

  iparam(:) = 0; dparam(:) = 0
  iparam(2) = 6; dparam(31) = 1.0e-16_dp 
  allocate(pre%FactorsColCoeffs(size(pre%Mat%ColCoeffs)))

  call dcsrilu0(mesh%NumCells, pre%Mat%ColCoeffs, pre%Mat%RowAddrs, &
    & pre%Mat%ColIndices, pre%FactorsColCoeffs, iparam, dparam, ierror)
  if (ierror /= 0) then
    error stop 'MKL `dcsrilu0` has failed, ierror='//I2S(ierror)
  end if

end subroutine InitPrecond_ILU0_MKL

subroutine ApplyPrecond_ILU0_MKL(pre, mesh, yArr, xArr, MatVec)
  class(tPreconditioner_ILU0_MKL), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  class(tArray), intent(inout), target :: xArr, yArr
  procedure(tMatVecFunc) :: MatVec

  type(tArray) :: tArr
  real(dp), pointer :: t(:), x(:), y(:)

  call AllocArray(tArr, mold=xArr)

  ! solve (L*U)*y = L*(U*y) = x
  ! solve L*t = x
  ! solve U*y = t

  call tArr%Get(t); call xArr%Get(x); call yArr%Get(y)

  call mkl_dcsrtrsv('L', 'N', 'U', mesh%NumCells, &
    & pre%FactorsColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, x, t)
  call mkl_dcsrtrsv('U', 'N', 'N', mesh%NumCells, &
    & pre%FactorsColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, t, y)

end subroutine ApplyPrecond_ILU0_MKL

#$end if

end module StormRuler_Precond_ILU_MKL
