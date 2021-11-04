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
module StormRuler_Solvers_Precond

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArrayR, AllocArray, ArrayAllocated

use StormRuler_BLAS, only: Fill, MatVecProd_Diagonal, Solve_Triangular
#$for T, _ in [SCALAR_TYPES[0]]
use StormRuler_BLAS, only: tMatVecFunc$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Preconditioner matrix-vector product function.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for T, typename in [SCALAR_TYPES[0]]
  subroutine tPreMatVecFunc$T(mesh, PuArr, uArr, MatVec, preEnv)
    import :: tMesh, tArray$T, tMatVecFunc$T
    class(tMesh), intent(inout) :: mesh
    class(tArray$T), intent(in), target :: uArr
    class(tArray$T), intent(inout), target :: PuArr
    procedure(tMatVecFunc$T) :: MatVec
    class(*), intent(inout), allocatable, target :: preEnv
  end subroutine tPreMatVecFunc$T
#$end for
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Jacobi preconditioner: 𝓟𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)⁻¹𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Precondition_Jacobi(mesh, PxArr, xArr, MatVec, preEnv)
  class(tMesh), intent(inout) :: mesh
  class(tArrayR), intent(in), target :: xArr
  class(tArrayR), intent(inout), target :: PxArr
  procedure(tMatVecFuncR) :: MatVec
  class(*), intent(inout), allocatable, target :: preEnv

  class(tArrayR), pointer :: diagArr
  real(dp), pointer :: x(:,:), Px(:,:), diag(:,:)

  ! ----------------------
  ! Cast preconditioner environment.
  ! ----------------------
  if (.not.allocated(preEnv)) then
    allocate(tArrayR :: preEnv)
  end if
  select type(preEnv)
    class is(tArrayR); diagArr => preEnv
    class default; error stop 'Invalid `preEnv`'
  end select

  ! ----------------------
  ! Build the preconditioner.
  ! ----------------------
  if (ArrayAllocated(diagArr)) then
    call AllocArray(diagArr, mold=xArr)

    ! TODO: this is not a correct diagonal extraction in block case!
    call Fill(mesh, PxArr, 1.0_dp)
    call MatVecProd_Diagonal(mesh, diagArr, PxArr, MatVec)
  end if

  ! ----------------------
  ! 𝓟𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)⁻¹𝒙.
  ! ----------------------
  ! TODO: this is not a correct diagonal solution in block case!
  Px(:,:) = x(:,:)/diag(:,:)

end subroutine Precondition_Jacobi

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! LU-SGS preconditioner: 𝓟𝒙 ← 𝘵𝘳𝘪𝘭(𝓐)⁻¹𝘥𝘪𝘢𝘨(𝓐)𝘵𝘳𝘪𝘶(𝓐)⁻¹𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Precondition_LU_SGS(mesh, PxArr, xArr, MatVec, preEnv)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  class(tArrayR), intent(in), target :: xArr
  class(tArrayR), intent(inout), target :: PxArr
  procedure(tMatVecFuncR) :: MatVec
  class(*), intent(inout), allocatable, target :: preEnv
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tArrayR), pointer :: diagArr
  type(tArrayR) :: yArr
  real(dp), pointer :: y(:,:), diag(:,:)

  call yArr%AllocMold(xArr)

  ! ----------------------
  ! Cast preconditioner environment.
  ! ----------------------
  if (.not.allocated(preEnv)) then
    allocate(tArrayR :: preEnv)
  end if
  select type(preEnv)
    class is(tArrayR); diagArr => preEnv
    class default; error stop 'Invalid `preEnv`'
  end select

  ! ----------------------
  ! Build the preconditioner.
  ! ----------------------
  if (ArrayAllocated(diagArr)) then
    call AllocArray(diagArr, mold=xArr)

    call Fill(mesh, PxArr, 1.0_dp)
    call MatVecProd_Diagonal(mesh, diagArr, PxArr, MatVec)
  end if

  ! ----------------------
  ! 𝒚 ← 𝘵𝘳𝘪𝘶(𝓐)⁻¹𝒙,
  ! 𝒚 ← 𝘥𝘪𝘢𝘨(𝓐)𝒚,
  ! 𝓟𝒙 ← 𝘵𝘳𝘪𝘭(𝓐)⁻¹𝒚.
  ! ----------------------
  call Solve_Triangular(mesh, 'U', yArr, xArr, diagArr, MatVec)
  ! TODO: this is not a correct diagonal multiplication in block case!
  y(:,:) = diag(:,:)*y(:,:)
  call Solve_Triangular(mesh, 'L', PxArr, yArr, diagArr, MatVec)

end subroutine Precondition_LU_SGS

end module StormRuler_Solvers_Precond
