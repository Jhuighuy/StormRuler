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
use StormRuler_Parameters, only: gMaxIterLU_SGS

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: tMatVecFunc, Set, Add

use StormRuler_Preconditioner, only: tMatrixBasedPreconditioner
use StormRuler_Matrix!, only: tMatrix, DiagMatrixVector

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

include 'mkl_spblas.fi'

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! LU Symmetric Gauss-Seidel (LU-SGS) preconditioner.
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
!! can be organized with a few LU-SGS iterations with 𝒃 = 𝒙.
!! For the simplicity, the term 𝓛𝓓⁻¹𝓤𝒙 is negleted 
!! on the first iteration.
!! In most practical cases, a single LU-SGS iteration
!! is enough for the considarable preconditioning.
!! 
!! Like the other triangular matrix-based preconditioners,
!! LU-SGS may suffer from poor parallel scaling.
!! 
!! References:
!! [1] ???
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

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set the LU-SGS preconditioner matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SetPreconditionerMatrix_LU_SGS(pre, mat)
  class(tPreconditioner_LU_SGS), intent(inout) :: pre
  class(tMatrix), intent(inout), target :: mat

  pre%Mat => mat

end subroutine SetPreconditionerMatrix_LU_SGS

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize the LU-SGS preconditioner: 𝓟 ← 𝘪𝘯𝘪𝘵(𝓐).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitPreconditioner_LU_SGS(pre, mesh, MatVec)
  class(tPreconditioner_LU_SGS), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  procedure(tMatVecFunc) :: MatVec

  if (.not.associated(pre%Mat)) then
    error stop 'Matrix for the LU-SGS preconditioner is not set.'
  end if

end subroutine InitPreconditioner_LU_SGS

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the LU-SGS preconditioner: 𝒚 ← 𝓟(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ApplyPreconditioner_LU_SGS(pre, mesh, yArr, xArr, MatVec)
  class(tPreconditioner_LU_SGS), intent(inout) :: pre
  class(tMesh), intent(in), target :: mesh
  class(tArray), intent(inout), target :: xArr, yArr
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: k
  type(tArray) :: tArr, sArr

  if (.not.associated(pre%Mat)) then
    error stop 'Matrix for the LU-SGS preconditioner is not set.'
  end if

  call AllocArray(tArr, mold=xArr)
  if (gMaxIterLU_SGS > 1) call AllocArray(sArr, mold=xArr)

  call Set(mesh, yArr, xArr)

  ! ----------------------
  ! First LU-SGS iteration:
  ! 𝒕 ← (𝓛 + 𝓓)⁻¹𝒙,
  ! 𝒕 ← 𝓓𝒕,
  ! 𝒚 ← (𝓓 + 𝓤)⁻¹𝒕.
  ! ----------------------
  block
    real(dp), pointer :: t(:), x(:), y(:)
    call tArr%Get(t); call xArr%Get(x); call yArr%Get(y)

    call mkl_dcsrtrsv('L', 'N', 'N', mesh%NumCells, &
      & pre%Mat%ColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, x, t)
    call DiagMatrixVector(mesh, pre%Mat, tArr, tArr)
    call mkl_dcsrtrsv('U', 'N', 'N', mesh%NumCells, &
      & pre%Mat%ColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, t, y)
  end block

  do k = 2, gMaxIterLU_SGS
    ! ----------------------
    ! Full LU-SGS iterations:
    ! 𝒔 ← 𝓤𝒚, 𝒔 ← 𝓓⁻¹𝒔, 𝒕 ← 𝓛𝒔, 𝒔 ← 𝒕 + 𝒙,
    ! 𝒕 ← (𝓛 + 𝓓)⁻¹𝒔,
    ! 𝒕 ← 𝓓𝒕,
    ! 𝒚 ← (𝓓 + 𝓤)⁻¹𝒕.
    ! ----------------------
    call UpperMatrixVector(mesh, pre%Mat, sArr, yArr)
    call InvDiagMatrixVector(mesh, pre%Mat, sArr, sArr)
    call UpperMatrixVector(mesh, pre%Mat, tArr, sArr)
    call Add(mesh, sArr, tArr, xArr)
    block
      real(dp), pointer :: t(:), s(:), y(:)
      call tArr%Get(t); call sArr%Get(s); call yArr%Get(y)
  
      call mkl_dcsrtrsv('L', 'N', 'N', mesh%NumCells, &
        & pre%Mat%ColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, s, t)
      call DiagMatrixVector(mesh, pre%Mat, tArr, tArr)
      call mkl_dcsrtrsv('U', 'N', 'N', mesh%NumCells, &
        & pre%Mat%ColCoeffs, pre%Mat%RowAddrs, pre%Mat%ColIndices, t, y)
    end block
  end do

end subroutine ApplyPreconditioner_LU_SGS

end module StormRuler_Preconditioner_LU_SGS
