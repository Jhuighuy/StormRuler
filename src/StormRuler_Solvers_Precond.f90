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
use StormRuler_BLAS, only: Fill, MatVecProd_Diagonal, Solve_Triangular
#$for T, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Preconditioner matrix-vector product function.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for type, typename in SCALAR_TYPES
  subroutine tPrecondFunc$type(mesh, Pu, u, MatVec, env, precond_env)
    import :: dp, tMesh, tMatVecFunc$type
    ! <<<<<<<<<<<<<<<<<<<<<<
    class(tMesh), intent(inout) :: mesh
    $typename, intent(in), target :: u(:,:)
    $typename, intent(inout), target :: Pu(:,:)
    procedure(tMatVecFunc$type) :: MatVec
    class(*), intent(inout) :: env
    class(*), intent(inout), allocatable, target :: precond_env
    ! >>>>>>>>>>>>>>>>>>>>>>
  end subroutine tPrecondFunc$type
#$end for
end interface

type :: tPrecondEnv_Diag
  real(dp), allocatable :: diag(:,:)
end type tPrecondEnv_Diag

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Jacobi preconditioner: 𝓟𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)⁻¹𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Precondition_Jacobi(mesh, Px, x, MatVec, env, precond_env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in), target :: x(:,:)
  real(dp), intent(inout), target :: Px(:,:)
  procedure(tMatVecFuncR) :: MatVec
  class(*), intent(inout) :: env
  class(*), intent(inout), allocatable, target :: precond_env
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tPrecondEnv_Diag), pointer :: diag_env

  ! ----------------------
  ! Cast preconditioner environment.
  ! ----------------------
  if (.not.allocated(precond_env)) then
    allocate(tPrecondEnv_Diag :: precond_env)
  end if
  select type(precond_env)
    class is(tPrecondEnv_Diag)
      diag_env => precond_env
  end select

  ! ----------------------
  ! Build the preconditioner.
  ! ----------------------
  if (.not.allocated(diag_env%diag)) then
    allocate(diag_env%diag, mold=x)

    ! TODO: this is not a correct diagonal extraction in block case!
    call Fill(mesh, Px, 1.0_dp)
    call MatVecProd_Diagonal(mesh, diag_env%diag, Px, MatVec, env)
  end if

  ! ----------------------
  ! 𝓟𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)⁻¹𝒙.
  ! TODO: this is not a correct diagonal solution in block case!
  ! ----------------------
  Px(:,:) = x(:,:)/diag_env%diag(:,:)

end subroutine Precondition_Jacobi

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! LU-SGS preconditioner: 𝓟𝒙 ← 𝘵𝘳𝘪𝘭(𝓐)⁻¹𝘥𝘪𝘢𝘨(𝓐)𝘵𝘳𝘪𝘶(𝓐)⁻¹𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Precondition_LU_SGS(mesh, Px, x, MatVec, env, precond_env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in), target :: x(:,:)
  real(dp), intent(inout), target :: Px(:,:)
  procedure(tMatVecFuncR) :: MatVec
  class(*), intent(inout) :: env
  class(*), intent(inout), allocatable, target :: precond_env
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tPrecondEnv_Diag), pointer :: diag_env
  real(dp), allocatable :: y(:,:)
  allocate(y, mold=x)

  ! ----------------------
  ! Cast preconditioner environment.
  ! ----------------------
  if (.not.allocated(precond_env)) then
    allocate(tPrecondEnv_Diag :: precond_env)
  end if
  select type(precond_env)
    class is(tPrecondEnv_Diag)
      diag_env => precond_env
  end select

  ! ----------------------
  ! Build the preconditioner.
  ! ----------------------
  if (.not.allocated(diag_env%diag)) then
    allocate(diag_env%diag, mold=x)

    call Fill(mesh, Px, 1.0_dp)
    call MatVecProd_Diagonal(mesh, diag_env%diag, Px, MatVec, env)
  end if

  ! ----------------------
  ! 𝒚 ← 𝘵𝘳𝘪𝘶(𝓐)⁻¹𝒙,
  ! 𝒚 ← 𝘥𝘪𝘢𝘨(𝓐)𝒚,
  ! 𝓟𝒙 ← 𝘵𝘳𝘪𝘭(𝓐)⁻¹𝒚.
  ! ----------------------
  call Solve_Triangular(mesh, y, x, diag_env%diag, 'U', MatVec, env)
  ! TODO: this is not a correct diagonal multiplication in block case!
  y(:,:) = diag_env%diag(:,:)*y(:,:)
  call Solve_Triangular(mesh, Px, y, diag_env%diag, 'L', MatVec, env)

end subroutine Precondition_LU_SGS

end module StormRuler_Solvers_Precond
