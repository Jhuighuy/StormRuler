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
  & Dot, Norm_2, Fill, Set, Scale, Add, Sub
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

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve a linear symmetric indefinite least squares problem:
!! minimize ‖Ax - b‖₂ over σ(r₀,Ar₀,A²r₀,…), where r₀ = b - Ax₀, 
!! using the Minimal Residual (MINRES) method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_MINRES_DE$rank(mesh, u, b, MatVec, env, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>

  logical :: first
  real(dp) :: alpha, beta, delta, &
    & gamma, gamma_bar, gamma_bar_bar
  real(dp), pointer :: r(@:,:), t(@:,:), &
    & p(@:,:), p_bar(@:,:), p_bar_bar(@:,:), &
    & s(@:,:), s_bar(@:,:), s_bar_bar(@:,:)

  first = .true.
  allocate(r, &
    & p, p_bar, p_bar_bar, &
    & s, s_bar, s_bar_bar, mold=u)

  ! ----------------------
  ! s ← Au,
  ! r ← b - s.
  ! ----------------------
  call MatVec(mesh, s, u, env)
  call Sub(mesh, r, b, s)

  ! ----------------------
  ! δ ← <r⋅r>,
  ! check convergence for √δ.
  ! ----------------------
  delta = Dot(mesh, r, r)
  if (params%Check(sqrt(delta))) return
  
  ! ----------------------
  ! p̅ ← r,
  ! s̅ ← Ap̅,
  ! γ̅ ← <s̅⋅s̅>.
  ! ----------------------
  call Set(mesh, p_bar, r)
  call MatVec(mesh, s_bar, p_bar, env)
  gamma_bar = Dot(mesh, s_bar, s_bar)

  do
    ! ----------------------
    ! α ← <r⋅s̅>/γ̅,
    ! u ← u + αp̅,
    ! r ← r - αs̅.
    ! ----------------------
    alpha = Dot(mesh, r, s_bar)/gamma_bar
    call Add(mesh, u, u, p_bar, alpha)
    call Sub(mesh, r, r, s_bar, alpha)

    ! ----------------------
    ! α ← <r⋅r>,
    ! check convergence for √α and √α/√δ.
    ! ----------------------
    alpha = Dot(mesh, r, r)
    if (params%Check(sqrt(alpha), sqrt(alpha/delta))) exit

    ! ----------------------
    ! p ← s̅,
    ! s ← As̅.
    ! ----------------------
    call Set(mesh, p, s_bar)
    call MatVec(mesh, s, s_bar, env)
    
    ! ----------------------
    ! β ← <s⋅s̅>/γ̅,
    ! p ← p - βp̅,
    ! s ← s - βs̅.
    ! ----------------------
    beta = Dot(mesh, s, s_bar)/gamma_bar
    call Sub(mesh, p, p, p_bar, beta)
    call Sub(mesh, s, s, s_bar, beta)

    if (.not.first) then

      ! ----------------------
      ! β ← <s⋅s̿>/γ̿,
      ! p ← p - βp̿,
      ! s ← s - βs̿.
      ! ----------------------
      beta = Dot(mesh, s, s_bar_bar)/gamma_bar_bar
      call Sub(mesh, p, p, p_bar_bar, beta)
      call Sub(mesh, s, s, s_bar_bar, beta)

    end if
    first = .false.

    ! ----------------------
    ! γ ← <s⋅s>.
    ! ----------------------
    gamma = Dot(mesh, s, s)

    ! ----------------------
    ! (γ, γ̅, γ̿) ← (γ̿, γ, γ̅),
    ! (p, p̅, p̿) ← (p̿, p, p̅),
    ! (s, s̅, s̿) ← (s̿, s, s̅).
    ! ----------------------
    gamma_bar_bar = gamma_bar; gamma_bar = gamma
    t => p_bar_bar; p_bar_bar => p_bar; p_bar => p; p => t
    t => s_bar_bar; s_bar_bar => s_bar; s_bar => s; s => t
  end do
end subroutine Solve_MINRES_DE$rank
#$end do

end module StormRuler_Solvers_LSQ
