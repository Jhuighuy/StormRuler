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
module StormRuler_Solvers_Base

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@, &
  & Set, Sub, Norm_2

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Preconsitioner matrix-vector product function.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$do rank = 0, NUM_RANKS
  subroutine tPrecondFunc$rank(mesh, Pu, u, MatVec, env, precond_env)
    import :: dp, tMesh, tMatVecFunc$rank
    ! <<<<<<<<<<<<<<<<<<<<<<
    class(tMesh), intent(inout) :: mesh
    real(dp), intent(in), target :: u(@:,:)
    real(dp), intent(inout), target :: Pu(@:,:)
    procedure(tMatVecFunc$rank) :: MatVec
    class(*), intent(inout) :: env
    class(*), intent(inout), allocatable, target :: precond_env
    ! >>>>>>>>>>>>>>>>>>>>>>
  end subroutine tPrecondFunc$rank
#$end do
end interface

interface ResidualNorm
#$do rank = 0, NUM_RANKS
  module procedure ResidualNorm$rank
#$end do
end interface ResidualNorm

interface ResidualNorm_Squared
#$do rank = 0, NUM_RANKS
  module procedure ResidualNorm_Squared$rank
#$end do
end interface ResidualNorm_Squared

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Compute squared residual norm: ‖b - A[P]y‖₂.
!! ----------------------------------------------------------------- !!
#$do rank = 0, NUM_RANKS
real(dp) function ResidualNorm$rank(mesh, y, b, MatVec, env, Precond, precond_env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:), y(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  procedure(tPrecondFunc$rank), optional :: Precond
  class(*), allocatable, intent(inout), optional :: precond_env
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: r(@:,:), x(@:,:)
  allocate(r, x, mold=y)

  ! ----------------------
  ! x ← Py,
  ! r ← Ax, OR: r ← Ay. 
  ! r ← b - r.
  ! OUT ← ‖r‖₂
  ! ----------------------
  if (present(Precond)) then
    call Precond(mesh, x, y, MatVec, env, precond_env)
    call MatVec(mesh, r, x, env)
  else
    call MatVec(mesh, r, y, env)
  end if
  call Sub(mesh, r, b, r)
  ResidualNorm$rank = Norm_2(mesh, r)

end function ResidualNorm$rank
#$end do

!! ----------------------------------------------------------------- !!
!! Compute squared residual norm: ‖(A[P])ᵀ(b - A[P]y)‖₂.
!! ----------------------------------------------------------------- !!
#$do rank = 0, NUM_RANKS
real(dp) function ResidualNorm_Squared$rank(mesh, y, b, &
    & MatVec, env, MatVec_T, env_T, Precond, precond_env, Precond_T, precond_env_T)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:), y(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec, MatVec_T
  class(*), intent(inout) :: env, env_T
  procedure(tPrecondFunc$rank), optional :: Precond, Precond_T
  class(*), allocatable, intent(inout), optional :: precond_env, precond_env_T
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: r(@:,:), x(@:,:)
  allocate(r, x, mold=y)

  ! ----------------------
  ! x ← Py,
  ! r ← Ax, OR: r ← Ay. 
  ! r ← b - r.
  ! ----------------------
  if (present(Precond)) then
    call Precond(mesh, x, y, MatVec, env, precond_env)
    call MatVec(mesh, r, x, env)
  else
    call MatVec(mesh, r, y, env)
  end if
  call Sub(mesh, r, b, r)

  ! ----------------------
  ! x ← Aᵀr, 
  ! r ← Pᵀx, OR: x ← r, r ← Aᵀx.
  ! OUT ← ‖r‖₂.
  ! ----------------------
  if (present(Precond)) then
    call MatVec_T(mesh, x, r, env_T)
    call Precond_T(mesh, r, x, MatVec_T, env_T, precond_env_T)
  else
    call Set(mesh, x, r)
    call MatVec_T(mesh, r, x, env_T)
  end if
  ResidualNorm_Squared$rank = Norm_2(mesh, r)

end function ResidualNorm_Squared$rank
#$end do

end module StormRuler_Solvers_Base
