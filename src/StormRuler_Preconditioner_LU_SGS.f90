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
module StormRuler_Preconditioner_LU_SGS

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: tMatVecFunc, Set

use StormRuler_Preconditioner, only: tMatrixBasedPreconditioner
use StormRuler_Matrix, only: tMatrix

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Incomplete LU preconditioner, MKL based.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tMatrixBasedPreconditioner) :: tPreconditioner_LU_SGS
  type(tMatrix), pointer, private :: Mat => null()

contains
  procedure :: SetMatrix => SetPreconditionerMatrix_LU_SGS
  procedure :: Init => InitPreconditioner_LU_SGS
  procedure :: Apply => ApplyPreconditioner_LU_SGS
end type tPreconditioner_LU_SGS

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

subroutine SetPreconditionerMatrix_LU_SGS(pre, mat)
  class(tPreconditioner_LU_SGS), intent(inout) :: pre
  class(tMatrix), intent(inout), target :: mat

  pre%Mat => mat

end subroutine SetPreconditionerMatrix_LU_SGS

subroutine InitPreconditioner_LU_SGS(pre, mesh, MatVec)
  class(tPreconditioner_LU_SGS), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  procedure(tMatVecFunc) :: MatVec

end subroutine InitPreconditioner_LU_SGS

subroutine ApplyPreconditioner_LU_SGS(pre, mesh, yArr, xArr, MatVec)
  class(tPreconditioner_LU_SGS), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  class(tArray), intent(inout), target :: xArr, yArr
  procedure(tMatVecFunc) :: MatVec

  call Set(mesh, yArr, xArr)

end subroutine ApplyPreconditioner_LU_SGS

end module StormRuler_Preconditioner_LU_SGS
