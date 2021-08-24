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
module StormRuler_Solvers_LSQ

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp
use StormRuler_Helpers, only: SafeDivide
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@, &
  & Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_Base, only: @{tPrecondFunc$$@|@0, NUM_RANKS}@

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_LSQR
#$do rank = 0, NUM_RANKS
  module procedure Solve_LSQR$rank
#$end do
end interface Solve_LSQR

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear least squares problem
!! minimize ‖A[P]y - b‖₂, where x = [P]y, using the LSQR method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_LSQR$rank(mesh, x, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: x(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFunc$rank), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  real(dp) :: alpha, beta, rho, rho_bar, &
    & theta, phi, phi_bar, phi_tilde, cc, ss
  real(dp), pointer :: s(@:,:), t(@:,:), &
    & r(@:,:), u(@:,:), v(@:,:), w(@:,:)
  class(*), allocatable :: precond_env
  
  allocate(s, t, r, u, v, w, mold=x)

  ! ----------------------
  ! β ← ‖b‖, u ← b/β,
  ! +----------+----------+
  ! | t ← Aᵀu, |          |
  ! | s ← Pᵀt, | s ← Aᵀu, |
  ! +----------+----------+
  ! α ← ‖s‖, v ← s/α,
  ! w ← v.
  ! ----------------------
  beta = Norm_2(mesh, b); call Scale(mesh, u, b, 1.0_dp/beta)
  if (present(Precond)) then
    call MatVec(mesh, t, u, env)
    call Precond(mesh, s, t, MatVec, env, precond_env)
  else
    call MatVec(mesh, s, u, env)
  end if
  alpha = Norm_2(mesh, s); call Scale(mesh, v, s, 1.0_dp/alpha)
  call Set(mesh, w, v)

  ! ----------------------
  ! x ← {0},
  ! ϕ̅ ← β, ρ̅ ← α.
  ! ----------------------
  call Fill(mesh, x, 0.0_dp)
  phi_bar = beta; rho_bar = alpha

  ! ----------------------
  ! ϕ̃ ← β,
  ! check convergence for ϕ̃.
  ! ----------------------
  phi_tilde = beta
  if (params%Check(phi_tilde)) return
  
  do
    ! ----------------------
    ! +----------+----------+
    ! | s ← Pv,  |          |
    ! | t ← As,  | t ← Pv   |
    ! +----------+----------+
    ! t ← t - αu,
    ! β ← ‖t‖, u ← t/β,
    ! +----------+----------+
    ! | t ← Aᵀu, |          |
    ! | s ← Pᵀt, | s ← Aᵀu, |
    ! +----------+----------+
    ! s ← s - βv,
    ! α ← ‖s‖, v ← s/α.
    ! ----------------------
    if (present(Precond)) then
      call Precond(mesh, s, v, MatVec, env, precond_env)
      call MatVec(mesh, t, s, env)
    else
      call MatVec(mesh, t, v, env)
    end if
    call Sub(mesh, t, t, u, alpha)
    beta = Norm_2(mesh, t); call Scale(mesh, u, t, 1.0_dp/beta)
    if (present(Precond)) then
      call MatVec(mesh, t, u, env)
      call Precond(mesh, s, t, MatVec, env, precond_env)
    else
      call MatVec(mesh, s, u, env)
    end if
    call Sub(mesh, s, s, v, beta)
    alpha = Norm_2(mesh, s); call Scale(mesh, v, s, 1.0_dp/alpha)
    
    ! ----------------------
    ! ρ ← √(ρ̅² + β²),
    ! c ← ρ̅/ρ, s ← β/ρ,
    ! θ ← sα, ρ̅ ← -cα,
    ! ϕ ← cϕ̅, ϕ̅ ← sϕ̅.
    ! ----------------------
    rho = hypot(rho_bar, beta)
    cc = rho_bar/rho; ss = beta/rho
    theta = ss*alpha; rho_bar = -cc*alpha
    phi = cc*phi_bar; phi_bar = ss*phi_bar

    ! ----------------------
    ! x ← x + (ϕ/ρ)w,
    ! w ← v - (θ/ρ)w.
    ! check convergence for ϕ̅ and ϕ̅/ϕ̃.
    ! ----------------------
    call Add(mesh, x, x, w, phi/rho)
    call Sub(mesh, w, v, w, theta/rho)
    if (params%Check(phi_bar, phi_bar/phi_tilde)) exit
  end do

  ! ----------------------
  ! t ← x,
  ! x ← Pt.
  ! ----------------------
  if (present(Precond)) then
    call Set(mesh, t, x)
    call Precond(mesh, x, t, MatVec, env, precond_env)
  end if
end subroutine Solve_LSQR$rank
#$end do

end module StormRuler_Solvers_LSQ
