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
module StormRuler_Precond

#$use 'StormRuler_Params.fi'

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray

use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_Matrix, only: tMatrix

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Preconditioner class.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, abstract :: tPreconditioner
contains
  procedure(tInitPrecondFunc), deferred :: Init
  procedure(tApplyPrecondFunc), deferred :: Apply
end type tPreconditioner

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize preconditioner: 𝓟 ← 𝘪𝘯𝘪𝘵(𝓐).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tInitPrecondFunc(pre, mesh, MatVec)
    import :: tMesh, tPreconditioner, tMatVecFunc
    class(tPreconditioner), intent(inout) :: pre
    class(tMesh), intent(in), target :: mesh
    procedure(tMatVecFunc) :: MatVec
  end subroutine tInitPrecondFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply preconditioner: 𝒚 ← 𝓟(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tApplyPrecondFunc(pre, mesh, yArr, xArr, MatVec)
    import :: tMesh, tArray, tPreconditioner, tMatVecFunc
    class(tPreconditioner), intent(inout) :: pre
    class(tMesh), intent(in), target :: mesh
    class(tArray), intent(inout), target :: xArr, yArr
    procedure(tMatVecFunc) :: MatVec
  end subroutine tApplyPrecondFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Matrix-based preconditioner class.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tPreconditioner), abstract :: tMatrixBasedPreconditioner
contains
  procedure(tPrecondSetMatrixFunc), deferred :: SetMatrix
end type tMatrixBasedPreconditioner

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set matrix for the matrix-based preconditioner.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tPrecondSetMatrixFunc(pre, mat)
    import :: tMatrixBasedPreconditioner, tMatrix
    class(tMatrixBasedPreconditioner), intent(inout) :: pre
    class(tMatrix), intent(inout), target :: mat
  end subroutine tPrecondSetMatrixFunc
end interface

end module StormRuler_Precond
