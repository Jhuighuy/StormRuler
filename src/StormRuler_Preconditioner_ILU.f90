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
use StormRuler_Parameters, only: gUseMKL

use StormRuler_Helpers, only: ErrorStop, I2S

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_Matrix, only: tMatrix
use StormRuler_Preconditioner, only: tMatrixBasedPreconditioner

#$use 'StormRuler_Macros.fi'
#$if HAS_MKL
use StormRuler_Libs_MKL, only: dcsrilu0, mkl_dcsrtrsv
#$end if

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Incomplete LU preconditioner, MKL based.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tMatrixBasedPreconditioner) :: tPreconditioner_ILU0_MKL
  type(tMatrix), pointer, private :: Mat => null()
  real(dp), allocatable, private :: FactorsColCoeffs(:)

contains
  procedure, non_overridable :: SetMatrix => SetPreconditionerMatrix_ILU0_MKL
  procedure, non_overridable :: Init => InitPreconditioner_ILU0_MKL
  procedure, non_overridable :: Apply => ApplyPreconditioner_ILU0_MKL
  
end type tPreconditioner_ILU0_MKL

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

subroutine SetPreconditionerMatrix_ILU0_MKL(pre, mat)
  class(tPreconditioner_ILU0_MKL), intent(inout) :: pre
  class(tMatrix), intent(inout), target :: mat

  pre%Mat => mat

end subroutine SetPreconditionerMatrix_ILU0_MKL

subroutine InitPreconditioner_ILU0_MKL(pre, mesh, MatVec)
  class(tPreconditioner_ILU0_MKL), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: iparam(128), ierror
  real(dp) :: dparam(128)

  iparam(:) = 0; dparam(:) = 0.0_dp
  iparam(2) = 6; dparam(31) = 1.0e-16_dp 
  allocate(pre%FactorsColCoeffs(size(pre%Mat%ColCoeffs)))

#$if HAS_MKL
  call dcsrilu0(mesh%NumCells, pre%Mat%ColCoeffs, pre%Mat%RowAddrs, &
    & pre%Mat%ColIndices, pre%FactorsColCoeffs, iparam, dparam, ierror)
  if (ierror /= 0) then
    call ErrorStop('MKL `dcsrilu0` has failed, ierror='//I2S(ierror))
  end if
#$else
  error stop 'ILU0_MKL requires MKL for now..'
#$end if

end subroutine InitPreconditioner_ILU0_MKL

subroutine ApplyPreconditioner_ILU0_MKL(pre, mesh, yArr, xArr, MatVec)
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

#$if HAS_MKL
  call mkl_dcsrtrsv('L', 'N', 'U', mesh%NumCells, &
    & pre%FactorsColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, x, t)
  call mkl_dcsrtrsv('U', 'N', 'N', mesh%NumCells, &
    & pre%FactorsColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, t, y)
#$else
  error stop 'ILU0_MKL requires MKL for now..'
#$end if
    
end subroutine ApplyPreconditioner_ILU0_MKL

end module StormRuler_Preconditioner_ILU
