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
module StormRuler_Preconditioner_AINV

use StormRuler_Consts, only: ip, dp

use StormRuler_Helpers, only: ErrorStop

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: tMatVecFunc, Set, Add

use StormRuler_Matrix, only: tMatrix, &
  & PartialMatrixVector, SolveDiag
use StormRuler_Preconditioner, only: tMatrixPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Approximate-inverse (AINV) preconditioner.
!!
!! Consider the symmetric sign definite (SPD) matrix 𝓐. 
!! It can be factored as 𝓐 = 𝓩𝓓𝓩ᵀ, such that 𝓓 is the (block-)diagonal 
!! matrix, 𝓩 = {𝒛₁,…,𝒛ₙ} is a (block-)column matrix, and 𝓩ᵀ𝓐𝓩 = 𝓓. 
!! [1] proposes an algorithm for computing the {𝒛ᵢ} and {𝑑ᵢᵢ}
!! (algorithm below is slightly modified version of [1]):
!!
!! 𝗳𝗼𝗿 𝑖 = 𝟣, 𝑛 𝗱𝗼:
!!   𝒛ᵢ ← 𝒆ᵢ, // here 𝒆ᵢ is the 𝑖-th column of 𝓘.
!!   𝑑ᵢᵢ ← 𝒂ᵢᵢ, // here 𝒂ᵢᵢ is the 𝑖-th diagonal element of 𝓐.
!! 𝗲𝗻𝗱 𝗳𝗼𝗿
!! 𝗳𝗼𝗿 𝑖 = 𝟣, 𝑛 - 𝟣 𝗱𝗼:
!!   𝗳𝗼𝗿 𝑗 = 𝑖 + 𝟣, 𝑛 𝗱𝗼:
!!     𝒛ⱼ ← 𝒛ⱼ - (𝑑ⱼⱼ/𝑑ᵢᵢ)𝒛ᵢ.
!!   𝗲𝗻𝗱 𝗳𝗼𝗿
!!   𝗳𝗼𝗿 𝑗 = 𝑖, 𝑛 𝗱𝗼:
!!     𝑑ⱼⱼ ← <𝒂ᵢ⋅𝒛ⱼ>, // here 𝒂ᵢ is the 𝑖-th row of 𝓐.
!!   𝗲𝗻𝗱 𝗳𝗼𝗿
!! 𝗲𝗻𝗱 𝗳𝗼𝗿
!!
!! Then 𝓐⁻¹ = 𝓩𝓓⁻¹𝓩ᵀ, and 𝒚 = 𝓐⁻¹𝒙 = 𝓩𝓓⁻¹𝓩ᵀ𝒙.
!!
!! Due to the fill-in phenomena, matrix 𝓩 needs to be sparsified.
!! In the current implementation we limit the sparsity pattern of 𝓩 
!! to the sparsity pattern of the powers of the original matrix 𝓐.
!! Original implementation in [1] uses the a approach. 
!!
!! AINV initialization algorithm:
!!
!! • Initialize 𝓩ᵀ as the identity matrix of with a portrait of 𝓐ᵖ,
!!   initialize 𝓓 as the matrix array of the diagonal of 𝓐.
!!
!! • Compute diagonal entries {𝑑ᵢᵢ} and sparse vectors {𝒛₁,…,𝒛ₙ} as 
!!   the rows of 𝓩ᵀ with the algorithm above [1]. 
!!   Now 𝓩ᵀ is the lower triangular matrix.
!!
!! • Symmetrize 𝓩 by filling it's upper triangular part.
!!
!! Application of the AINV is requires two matrix-vector product and
!! a diagonal scaling, so a good parallel scaling is expected.
!! However the initialization phase is sequential and may be 
!! time consuming.
!!
!! References:
!! [1] Michele Benzi, Carl Dean Meyer and Miroslav Tuma. 
!!     “A Sparse Approximate Inverse Preconditioner for the 
!!      Conjugate Gradient Method.” 
!!     SIAM J. Sci. Comput. 17 (1996): 1135-1149.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tMatrixPreconditioner) :: tAinvPreconditioner
  type(tMatrix), pointer, private :: Mat => null()
  type(tArray), private :: DMat
  type(tMatrix), private :: ZMat

contains
  procedure, non_overridable :: SetMatrix => SetAinvPreconditionerMatrix
  procedure, non_overridable :: Init => InitAinvPreconditioner
  procedure, non_overridable :: Apply => ApplyAinvPreconditioner

end type tAinvPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set the AINV preconditioner matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SetAinvPreconditionerMatrix(pre, mat)
  class(tAinvPreconditioner), intent(inout) :: pre
  class(tMatrix), intent(inout), target :: mat

  pre%Mat => mat

end subroutine SetAinvPreconditionerMatrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize the AINV preconditioner: 𝓟 ← 𝘪𝘯𝘪𝘵(𝓐).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitAinvPreconditioner(pre, mesh, MatVec)
  class(tAinvPreconditioner), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  procedure(tMatVecFunc) :: MatVec

  if (.not.associated(pre%Mat)) then
    call ErrorStop('Matrix for the AINV preconditioner is not set.')
  end if

end subroutine InitAinvPreconditioner

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the AINV preconditioner: 𝒚 ← 𝓟(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ApplyAinvPreconditioner(pre, mesh, yArr, xArr, MatVec)
  class(tAinvPreconditioner), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  class(tArray), intent(inout), target :: xArr, yArr
  procedure(tMatVecFunc) :: MatVec

  if (.not.associated(pre%Mat)) then
    call ErrorStop('Matrix for the AINV preconditioner is not set.')
  end if

  call Set(mesh, yArr, xArr)

end subroutine ApplyAinvPreconditioner

end module StormRuler_Preconditioner_AINV
