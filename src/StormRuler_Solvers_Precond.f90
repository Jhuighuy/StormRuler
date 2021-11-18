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
use StormRuler_Array, only: tArray, AllocArray, ArrayAllocated

use StormRuler_BLAS, only: tMatVecFunc
use StormRuler_BLAS, only: Fill, MatVecProd_Diagonal, Solve_Triangular

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Precondition_Jacobi
  module procedure Precondition_Jacobi
end interface Precondition_Jacobi

interface Precondition_LU_SGS
  module procedure Precondition_LU_SGS
end interface Precondition_LU_SGS

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Preconditioner matrix-vector product function.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tPreMatVecFunc(mesh, PuArr, uArr, MatVec, preEnv)
    import :: tMesh, tArray, tMatVecFunc
    class(tMesh), intent(inout) :: mesh
    class(tArray), intent(in), target :: uArr
    class(tArray), intent(inout), target :: PuArr
    procedure(tMatVecFunc) :: MatVec
    class(*), intent(inout), allocatable, target :: preEnv
  end subroutine tPreMatVecFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Jacobi preconditioner: 𝓟𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)⁻¹𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Precondition_Jacobi(mesh, PxArr, xArr, MatVec, preEnv)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in), target :: xArr
  class(tArray), intent(inout), target :: PxArr
  procedure(tMatVecFunc) :: MatVec
  class(*), intent(inout), allocatable, target :: preEnv

  class(tArray), pointer :: diagArr
  real(dp), pointer :: x(:,:), Px(:,:), diag(:,:)

  ! ----------------------
  ! Cast preconditioner environment.
  ! ----------------------
  if (.not.allocated(preEnv)) then
    allocate(tArray :: preEnv)
  end if
  select type(preEnv)
    class is(tArray); diagArr => preEnv
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

  call xArr%Get(x); call PxArr%Get(Px); call diagArr%Get(diag)

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
  class(tArray), intent(in), target :: xArr
  class(tArray), intent(inout), target :: PxArr
  procedure(tMatVecFunc) :: MatVec
  class(*), intent(inout), allocatable, target :: preEnv
  ! >>>>>>>>>>>>>>>>>>>>>>

  type(tArray) :: yArr
  class(tArray), pointer :: diagArr
  real(dp), pointer :: y(:,:), diag(:,:)

  ! ----------------------
  ! Cast preconditioner environment.
  ! ----------------------
  if (.not.allocated(preEnv)) then
    allocate(tArray :: preEnv)
  end if
  select type(preEnv)
    class is(tArray); diagArr => preEnv
    class default; error stop 'Invalid `preEnv`'
  end select

  ! ----------------------
  ! Build the preconditioner.
  ! ----------------------
  if (.not.ArrayAllocated(diagArr)) then
    call AllocArray(diagArr, mold=xArr)

    call Fill(mesh, PxArr, 1.0_dp)
    call MatVecProd_Diagonal(mesh, diagArr, PxArr, MatVec)
  end if

  call AllocArray(yArr, mold=xArr)
  call yArr%Get(y); 
  call diagArr%Get(diag)

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
