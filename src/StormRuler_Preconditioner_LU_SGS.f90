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

use StormRuler_Consts, only: ip, dp
use StormRuler_Parameters, only: gMaxIterLU_SGS

use StormRuler_Helpers, only: ErrorStop

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: tMatVecFunc, Set, Add

use StormRuler_Matrix, only: tMatrix, tParallelTriangularContext, &
  & PartialMatrixVector, SolveDiag, &
  & InitParallelTriangularContext, SolveTriangular
use StormRuler_Preconditioner, only: tMatrixBasedPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Block LU Symmetric Gauss-Seidel (LU-SGS) preconditioner.
!!
!! Consider the decomposition: 𝓐 = 𝓛 + 𝓓 + 𝓤, where 𝓓 is the 
!! (block-)diagonal of 𝓐, 𝓛 and 𝓤 are lower and upper strict
!! (block-)triangular parts of 𝓐. Then the equation 𝓐𝒚 = 𝒃 can 
!! be factored as:
!! 
!! (𝓛 + 𝓓)𝓓⁻¹(𝓓 + 𝓤)𝒚 = 𝒃 + 𝓛𝓓⁻¹𝓤𝒚.
!! 
!! The classical LU-SGS iterations for solving the equation 𝓐𝒚 = 𝒃
!! are based on following algorithm:
!! 
!! (𝓛 + 𝓓)𝒕 = 𝒃 + 𝓛𝓓⁻¹𝓤𝒚,
!! (𝓓 + 𝓤)𝒚̂ = 𝓓𝒕,
!! 
!! where 𝒚 and 𝒚̂ are the current and updated solution vectors.
!! 
!! In preconditioning, computation of the vector 𝒚 = 𝓟𝒙 ≈ 𝓐⁻¹𝒙 
!! can be organized with a few LU-SGS iterations with 𝒃 = 𝒙 
!! with zero initial guess. The number of inner iterations is 
!! controlled by the parameter `MaxIterLU_SGS`. In most practical 
!! cases, a single LU-SGS iteration is enough for the considarable 
!! preconditioning. More than 5 LU-SGS iterations is not recommended 
!! due to the numerical instability issues.
!! 
!! LU-SGS preconditioner preserves symmetry, so it may be used
!! with the self-adjoint solvers like CG or MINRES. 
!! 
!! Like the other triangular matrix-based preconditioners, LU-SGS 
!! may suffer from poor parallel scaling.
!! 
!! References:
!! [1] ???
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tMatrixBasedPreconditioner) :: tLuSgsPreconditioner
  type(tMatrix), pointer, private :: Mat => null()
  type(tParallelTriangularContext), private :: LowerCtx, UpperCtx

contains
  procedure, non_overridable :: SetMatrix => SetLuSgsPreconditionerMatrix
  procedure, non_overridable :: Init => InitLuSgsPreconditioner
  procedure, non_overridable :: Apply => ApplyLuSgsPreconditioner

end type tLuSgsPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set the LU-SGS preconditioner matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SetLuSgsPreconditionerMatrix(pre, mat)
  class(tLuSgsPreconditioner), intent(inout) :: pre
  class(tMatrix), intent(inout), target :: mat

  pre%Mat => mat

end subroutine SetLuSgsPreconditionerMatrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize the LU-SGS preconditioner: 𝓟 ← 𝘪𝘯𝘪𝘵(𝓐).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitLuSgsPreconditioner(pre, mesh, MatVec)
  class(tLuSgsPreconditioner), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  procedure(tMatVecFunc) :: MatVec

  if (.not.associated(pre%Mat)) then
    call ErrorStop('Matrix for the LU-SGS preconditioner is not set.')
  end if

  if (.not.allocated(pre%LowerCtx%LevelAddrs)) &
    & call InitParallelTriangularContext(mesh, 'L', pre%mat, pre%LowerCtx)
  if (.not.allocated(pre%UpperCtx%LevelAddrs)) &
    & call InitParallelTriangularContext(mesh, 'U', pre%mat, pre%UpperCtx)

end subroutine InitLuSgsPreconditioner

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the LU-SGS preconditioner: 𝒚 ← 𝓟(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ApplyLuSgsPreconditioner(pre, mesh, yArr, xArr, MatVec)
  class(tLuSgsPreconditioner), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  class(tArray), intent(inout), target :: xArr, yArr
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: k
  type(tArray) :: tArr

  if (.not.associated(pre%Mat)) then
    call ErrorStop('Matrix for the LU-SGS preconditioner is not set.')
  end if

  call AllocArray(tArr, mold=xArr)

  ! ----------------------
  ! First LU-SGS iteration:
  ! 𝒚 ← (𝓛 + 𝓓)⁻¹𝒙,
  ! 𝒕 ← 𝓓𝒚,
  ! 𝒚 ← (𝓓 + 𝓤)⁻¹𝒕.
  ! ----------------------
  call SolveTriangular(mesh, 'L', pre%Mat, pre%LowerCtx, yArr, xArr)
  call PartialMatrixVector(mesh, 'D', pre%Mat, tArr, yArr)
  call SolveTriangular(mesh, 'U', pre%Mat, pre%UpperCtx, yArr, tArr)

  do k = 2, gMaxIterLU_SGS
    ! ----------------------
    ! Full LU-SGS iterations:
    ! 𝒕 ← 𝓤𝒚, 
    ! 𝒚 ← 𝓓⁻¹𝒕, 
    ! 𝒕 ← 𝓛𝒚, 
    ! 𝒕 ← 𝒕 + 𝒙,
    ! ----------------------
    call PartialMatrixVector(mesh, 'U', pre%Mat, tArr, yArr)
    call SolveDiag(mesh, pre%Mat, yArr, tArr)
    call PartialMatrixVector(mesh, 'L', pre%Mat, tArr, yArr)
    call Add(mesh, tArr, tArr, xArr)

    ! ----------------------
    ! 𝒚 ← (𝓛 + 𝓓)⁻¹𝒕,
    ! 𝒕 ← 𝓓𝒚,
    ! 𝒚 ← (𝓓 + 𝓤)⁻¹𝒕.
    ! ----------------------
    call SolveTriangular(mesh, 'L', pre%Mat, pre%LowerCtx, yArr, tArr)
    call PartialMatrixVector(mesh, 'D', pre%Mat, tArr, yArr)
    call SolveTriangular(mesh, 'U', pre%Mat, pre%UpperCtx, yArr, tArr)
  end do

end subroutine ApplyLuSgsPreconditioner

end module StormRuler_Preconditioner_LU_SGS
