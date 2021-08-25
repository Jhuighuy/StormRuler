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
module StormRuler_Solvers_LSQR

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
!! Solve a right preconditioned linear least squares problem:
!! minimize ‖A[P]y - b‖₂, where x = [P]y, using the LSQR method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_LSQR$rank(mesh, x, b, &
    & MatVec, env, MatVec_T, env_T, params, Precond, Precond_T)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: x(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec, MatVec_T
  class(*), intent(inout) :: env, env_T
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFunc$rank), optional :: Precond, Precond_T
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  ! ----------------------
  ! [1] Paige, C. and M. Saunders. 
  !     “LSQR: An Algorithm for Sparse Linear Equations and Sparse Least Squares.” 
  !     ACM Trans. Math. Softw. 8 (1982): 43-71.
  ! [2] Karimi, S., D. K. Salkuyeh and F. Toutounian. 
  !     “A preconditioner for the LSQR algorithm.” 
  !     Journal of applied mathematics & informatics 26 (2008): 213-222.
  ! ----------------------

  real(dp) :: alpha, beta, rho, rho_bar, &
    & theta, phi, phi_bar, phi_tilde, cs, sn
  real(dp), pointer :: s(@:,:), t(@:,:), u(@:,:), v(@:,:), w(@:,:)
  class(*), allocatable :: precond_env, precond_env_T
  
  allocate(s, t, u, v, w, mold=x)

  ! ----------------------
  ! β ← ‖b‖, u ← b/β,
  ! t ← Aᵀu,
  ! s ← Pᵀt, OR: s ← Aᵀu,
  ! α ← ‖s‖, v ← s/α,
  ! w ← v.
  ! ----------------------
  beta = Norm_2(mesh, b); call Scale(mesh, u, b, 1.0_dp/beta)
  if (present(Precond)) then
    call MatVec_T(mesh, t, u, env_T)
    call Precond_T(mesh, s, t, MatVec_T, env_T, precond_env_T)
  else
    call MatVec_T(mesh, s, u, env_T)
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
    ! s ← Pv,
    ! t ← As, OR: t ← Pv
    ! t ← t - αu,
    ! β ← ‖t‖, u ← t/β,
    ! t ← Aᵀu,
    ! s ← Pᵀt, OR: s ← Aᵀu,
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
      call MatVec_T(mesh, t, u, env_T)
      call Precond_T(mesh, s, t, MatVec_T, env_T, precond_env_T)
    else
      call MatVec_T(mesh, s, u, env_T)
    end if
    call Sub(mesh, s, s, v, beta)
    alpha = Norm_2(mesh, s); call Scale(mesh, v, s, 1.0_dp/alpha)
    
    ! ----------------------
    ! ρ ← √(ρ̅² + β²),
    ! cs ← ρ̅/ρ, sn ← β/ρ,
    ! θ ← ssα, ρ̅ ← -ccα,
    ! ϕ ← ccϕ̅, ϕ̅ ← ssϕ̅.
    ! ----------------------
    rho = hypot(rho_bar, beta)
    cs = rho_bar/rho; sn = beta/rho
    theta = sn*alpha; rho_bar = -cs*alpha
    phi = cs*phi_bar; phi_bar = sn*phi_bar

    ! ----------------------
    ! x ← x + (ϕ/ρ)w,
    ! w ← v - (θ/ρ)w.
    ! check convergence for ϕ̅ and ϕ̅/ϕ̃.
    ! ( ϕ̅ and ϕ̃ implicitly contain residual norms. )
    ! ----------------------
    call Add(mesh, x, x, w, phi/rho)
    call Sub(mesh, w, v, w, theta/rho)
    if (params%Check(phi_bar, phi_bar/phi_tilde)) exit
  end do

  ! ----------------------
  ! t ← x,
  ! x ← Pt. OR: do nothing.
  ! ----------------------
  if (present(Precond)) then
    call Set(mesh, t, x)
    call Precond(mesh, x, t, MatVec, env, precond_env)
  end if
end subroutine Solve_LSQR$rank
#$end do

end module StormRuler_Solvers_LSQR
