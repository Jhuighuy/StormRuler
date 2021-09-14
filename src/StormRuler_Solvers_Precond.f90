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
use StormRuler_BLAS, only: @{tMatVecFuncR$$@|@0, NUM_RANKS}@, &
  & @{tMatVecFuncC$$@|@0, NUM_RANKS}@, Fill, MatVecProd_Diagonal, Solve_Triangular

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Preconditioner matrix-vector product function.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$do rank = 0, NUM_RANKS
#$for type, typename in SCALAR_TYPES
  subroutine tPrecondFunc$type$rank(mesh, Pu, u, MatVec, env, precond_env)
    import :: dp, tMesh, tMatVecFunc$type$rank
    ! <<<<<<<<<<<<<<<<<<<<<<
    class(tMesh), intent(inout) :: mesh
    $typename, intent(in), target :: u(@:,:)
    $typename, intent(inout), target :: Pu(@:,:)
    procedure(tMatVecFunc$type$rank) :: MatVec
    class(*), intent(inout) :: env
    class(*), intent(inout), allocatable, target :: precond_env
    ! >>>>>>>>>>>>>>>>>>>>>>
  end subroutine tPrecondFunc$type$rank
#$end for
#$end do
end interface

#$do rank = 0, NUM_RANKS
type :: tPrecondEnv_Diag$rank
  real(dp), allocatable :: diag(@:,:)
end type !tPrecondEnv_Diag$rank
#$end do

interface Precondition_Jacobi
#$do rank = 0, NUM_RANKS
  module procedure Precondition_Jacobi$rank
#$end do
end interface Precondition_Jacobi

interface Precondition_LU_SGS
#$do rank = 0, NUM_RANKS
  module procedure Precondition_LU_SGS$rank
#$end do
end interface Precondition_LU_SGS

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Jacobi preconditioner: 𝓟𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)⁻¹𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Precondition_Jacobi$rank(mesh, Px, x, MatVec, env, precond_env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in), target :: x(@:,:)
  real(dp), intent(inout), target :: Px(@:,:)
  procedure(tMatVecFuncR$rank) :: MatVec
  class(*), intent(inout) :: env
  class(*), intent(inout), allocatable, target :: precond_env
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tPrecondEnv_Diag$rank), pointer :: diag_env

  ! ----------------------
  ! Cast preconditioner environment.
  ! ----------------------
  if (.not.allocated(precond_env)) then
    allocate(tPrecondEnv_Diag$rank :: precond_env)
  end if
  select type(precond_env)
    class is(tPrecondEnv_Diag$rank)
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
  Px(@:,:) = -x(@:,:)/diag_env%diag(@:,:)

end subroutine Precondition_Jacobi$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! LU-SGS preconditioner: 𝓟𝒙 ← 𝘵𝘳𝘪𝘭(𝓐)⁻¹𝘥𝘪𝘢𝘨(𝓐)𝘵𝘳𝘪𝘶(𝓐)⁻¹𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Precondition_LU_SGS$rank(mesh, Px, x, MatVec, env, precond_env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in), target :: x(@:,:)
  real(dp), intent(inout), target :: Px(@:,:)
  procedure(tMatVecFuncR$rank) :: MatVec
  class(*), intent(inout) :: env
  class(*), intent(inout), allocatable, target :: precond_env
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tPrecondEnv_Diag$rank), pointer :: diag_env
  real(dp), allocatable :: y(@:,:)
  allocate(y, mold=x)

  ! ----------------------
  ! Cast preconditioner environment.
  ! ----------------------
  if (.not.allocated(precond_env)) then
    allocate(tPrecondEnv_Diag$rank :: precond_env)
  end if
  select type(precond_env)
    class is(tPrecondEnv_Diag$rank)
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
  y(@:,:) = diag_env%diag(@:,:)*y(@:,:)
  call Solve_Triangular(mesh, Px, y, diag_env%diag, 'L', MatVec, env)

end subroutine Precondition_LU_SGS$rank
#$end do

end module StormRuler_Solvers_Precond
