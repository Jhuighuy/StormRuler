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
module StormRuler_Solvers_MINRES

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@, &
  & Dot, Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_Base, only: @{tPrecondFunc$$@|@0, NUM_RANKS}@
use StormRuler_Solvers_Base, only: ResidualNorm

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_MINRES
#$do rank = 0, NUM_RANKS
  module procedure Solve_MINRES$rank
#$end do
end interface Solve_MINRES

#$if False
interface Solve_MINRES_QLP
#$do rank = 0, NUM_RANKS
  module procedure Solve_MINRES_QLP$rank
#$end do
end interface Solve_MINRES_QLP
#$end if

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

subroutine SymOrtho(a, b, cs, sn, rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), value :: a, b
  real(dp), intent(out) :: cs, sn, rr
  ! >>>>>>>>>>>>>>>>>>>>>>

  rr = hypot(a, b)
  if (rr > 0.0_dp) then
    cs = a/rr; sn = b/rr
  else
    cs = 1.0_dp; sn = 0.0_dp
  end if
end subroutine SymOrtho

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint indefinite operator equation: 
!! [P¹ᐟ²]A[P¹ᐟ²]y = [P¹ᐟ²]b, [P¹ᐟ²]y = x, using the MINRES method.
!! ( P = Pᵀ > 0 is required. )
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_MINRES$rank(mesh, x, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: x(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFunc$rank), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>

  ! ----------------------
  ! [1] Paige, C. and M. Saunders. 
  !     “Solution of Sparse Indefinite Systems of Linear Equations.” 
  !     SIAM Journal on Numerical Analysis 12 (1975): 617-629.
  ! [2] Choi, S.-C. T.
  !     “Iterative Methods for Singular Linear Equations and Least-Squares Problems” 
  !     PhD thesis, ICME, Stanford University.
  ! ----------------------

  real(dp) :: alpha, beta, beta_dot, gamma, &
    & delta, delta_dot, epsilon, epsilon_dot, &
    & tau, phi, phi_tilde, cs, sn
  real(dp), pointer :: tmp(@:,:), &
    & p(@:,:), q(@:,:), q_dot(@:,:), &
    & w(@:,:), w_dot(@:,:), w_ddot(@:,:), &
    & z(@:,:), z_dot(@:,:), z_ddot(@:,:)
  class(*), allocatable :: precond_env

  allocate(p, w, w_dot, w_ddot, z, z_dot, z_ddot, mold=x)
  if (present(Precond)) then
    allocate(q, q_dot, mold=x)
  end if

  ! ----------------------
  ! Initialize:
  ! ẇ ← {0}ᵀ,
  ! ẅ ← {0}ᵀ,
  ! ż ← Ax,     // Modification to
  ! ż ← b - ż,  // utilize the initial guess.
  ! z̈ ← 0,
  ! q ← Pż, OR: q ← ż,
  ! β̇ ← 1, β ← √<q⋅ż>,
  ! ϕ ← β, δ ← 0, ϵ ← 0,
  ! cs ← -1, sn ← 0.
  ! ----------------------
  call Fill(mesh, w_dot, 0.0_dp)
  call Fill(mesh, w_ddot, 0.0_dp)
  call MatVec(mesh, z_dot, x, env)
  call Sub(mesh, z_dot, b, z_dot)
  call Fill(mesh, z_ddot, 0.0_dp)
  if (present(Precond)) then
    call Precond(mesh, q, z_dot, MatVec, env, precond_env)
  else
    q => z_dot
  end if
  beta_dot = 1.0_dp; beta = sqrt(Dot(mesh, q, z_dot))
  phi = beta; delta = 0.0_dp; epsilon = 0.0_dp
  cs = -1.0_dp; sn = 0.0_dp

  ! ----------------------
  ! ϕ̃ ← ϕ,
  ! Check convergence for ϕ̃.
  ! ----------------------
  phi_tilde = phi
  if (params%Check(phi_tilde)) return

  do
    ! ----------------------
    ! Continue the Lanczos process:
    ! p ← Aq,
    ! α ← <q⋅p>/β²,
    ! z ← (1/β)p - (α/β)ż,
    ! z ← z - (β/β̇)z̈,
    ! q̇ ← q,
    ! q ← Pz, OR: q ← z,
    ! β̇ ← β, β ← √<q⋅z>,
    ! z̈ ← ż, ż ← z.
    ! ----------------------
    call MatVec(mesh, p, q, env)
    alpha = Dot(mesh, q, p)/(beta**2)
    call Sub(mesh, z, p, z_dot, alpha/beta, 1.0_dp/beta)
    call Sub(mesh, z, z, z_ddot, beta/beta_dot)
    if (present(Precond)) then
      tmp => q_dot; q_dot => q; q => tmp
      call Precond(mesh, q, z, MatVec, env, precond_env)
    else
      q_dot => q; q => z
    end if
    beta_dot = beta; beta = sqrt(Dot(mesh, q, z))
    tmp => z_ddot; z_ddot => z_dot; z_dot => z; z => tmp

    ! ----------------------
    ! Construct and apply rotations:
    ! δ̇ ← cs⋅δ + sn⋅α, γ ← sn⋅δ - cs⋅α,
    ! ϵ̇ ← ϵ, ϵ ← sn⋅β, δ ← -cs⋅β,
    ! cs, sn, γ ← SymOrtho(γ, β),
    ! τ ← cs⋅ϕ, ϕ ← sn⋅ϕ.
    ! ----------------------
    delta_dot = cs*delta + sn*alpha; gamma = sn*delta - cs*alpha
    epsilon_dot = epsilon; epsilon = sn*beta; delta = -cs*beta
    call SymOrtho(gamma*1.0_dp, beta, cs, sn, gamma)
    tau = cs*phi; phi = sn*phi
    
    ! ----------------------
    ! Update solution:
    ! w ← (1/(β̇γ))q̇ - (δ̇/γ)ẇ,
    ! w ← w - (ϵ̇/γ)ẅ,
    ! x ← x + τw,
    ! ẅ ← ẇ, ẇ ← w.
    ! ----------------------
    call Sub(mesh, w, q_dot, w_dot, delta_dot/gamma, 1.0_dp/(beta_dot*gamma))
    call Sub(mesh, w, w, w_ddot, epsilon_dot/gamma)
    call Add(mesh, x, x, w, tau)
    tmp => w_ddot; w_ddot => w_dot; w_dot => w; w => tmp

    ! ----------------------
    ! Check convergence for ϕ and ϕ/ϕ̃.
    ! ( ϕ and ϕ̃ implicitly contain residual norms. )
    ! ----------------------
    if (params%Check(phi, phi/phi_tilde)) exit
  end do
end subroutine Solve_MINRES$rank
#$end do

!!             !!
!! NOT WORKING !!
!!             !!
#$if False
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint indefinite operator equation: 
!! [P¹ᐟ²]A[P¹ᐟ²]y = [P¹ᐟ²]b, [P¹ᐟ²]y = x, using the MINRES-QLP method.
!! ( P = Pᵀ > 0 is required. )
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_MINRES_QLP$rank(mesh, x, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: x(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFunc$rank), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>

  ! ----------------------
  ! [1] Choi, S.-C. T.
  !     “Iterative Methods for Singular Linear Equations and Least-Squares Problems” 
  !     PhD thesis, ICME, Stanford University.
  ! ----------------------

  real(dp) :: alpha, &
    & beta, beta_dot, &
    & gamma(6), gamma_dot(6), gamma_ddot(6), &
    & delta(3), &
    & epsilon, epsilon_dot, &
    & tau, tau_dot, tau_ddot, &
    & eta, eta_dot, &
    & theta(2), theta_dot(2), theta_ddot(2), &
    & mu(3), mu_dot(3), mu_ddot(3), mu_dddot(3), &
    & phi, phi_tilde, &
    & cs(3), sn(3)
  real(dp), pointer :: tmp(@:,:), p(@:,:), &
    & q(@:,:), q_dot(@:,:), &
    & w(@:,:), w_dot(@:,:), w_ddot(@:,:), &
    & z(@:,:), z_dot(@:,:), z_ddot(@:,:), &
    & x_dot(@:,:), x_ddot(@:,:), x_dddot(@:,:)
  integer(ip) :: i, k
  class(*), allocatable :: precond_env

  allocate(p, z, z_dot, z_ddot, x_dot, x_ddot, x_dddot, mold=x)
  allocate(w, w_dot, w_ddot, mold=x)
  if (present(Precond)) then
    allocate(q, q_dot, mold=x)
  end if

  ! ----------------------
  ! Initialize:
  ! ẇ ← {0}ᵀ,
  ! ẅ ← {0}ᵀ,
  ! x ← {0}ᵀ,
  ! ẋ ← {0}ᵀ,
  ! ẍ ← {0}ᵀ,
  ! ż ← b,
  ! z̈ ← 0,
  ! q ← Pż, OR: q ← ż,
  ! β̇ ← 1, β ← √<q⋅ż>,
  ! ϕ ← β, δ₁ ← 0,
  ! τ̇ ← 0,
  ! γ̇ ← 0, γ̈ ← 0,
  ! η ← 0, η̇ ← 0,
  ! θ ← 0, θ̇ ← 0, θ̈ ← 0,
  ! μ̇ ← 0, μ̈ ← 0,
  ! cs₁ ← -1, cs₂ ← -1, cs₃ ← -1, 
  ! sn₁ ←  0, sn₂ ←  0, sn₃ ←  0,
  ! k ← 1
  ! ----------------------
  call Fill(mesh, w_dot, 0.0_dp)
  call Fill(mesh, w_ddot, 0.0_dp)
  call Fill(mesh, x, 0.0_dp)
  call Fill(mesh, x_dot, 0.0_dp)
  call Fill(mesh, x_ddot, 0.0_dp)
  call Set(mesh, z_dot, b)
  call Fill(mesh, z_ddot, 0.0_dp)
  if (present(Precond)) then
    call Precond(mesh, q, z_dot, MatVec, env, precond_env)
  else
    q => z_dot
  end if
  beta_dot = 1.0_dp; beta = sqrt(Dot(mesh, q, z_dot))
  phi = beta; delta(1) = 0.0_dp
  tau_dot = 0.0_dp
  gamma_dot = 0.0_dp; gamma_ddot = 0.0_dp
  eta = 0.0_dp; eta_dot = 0.0_dp
  theta = 0.0_dp; theta_dot = 0.0_dp; theta_ddot = 0.0_dp
  mu_dot = 0.0_dp; mu_ddot = 0.0_dp
  cs(:) = -1.0_dp; sn(:) =  0.0_dp
  k = 1

  ! ----------------------
  ! ϕ̃ ← ϕ,
  ! Check convergence for ϕ̃.
  ! ----------------------
  phi_tilde = phi
  print *, '000', ResidualNorm(mesh, x, b, MatVec, env)
  if (params%Check(phi_tilde)) return

  do
    ! ----------------------
    ! Continue the Lanczos process:
    ! p ← Aq,
    ! α ← <q⋅p>/β²,
    ! z ← (1/β)p - (α/β)ż,
    ! z ← z - (β/β̇)z̈,
    ! q̇ ← q,
    ! q ← Pz, OR: q ← z,
    ! β̇ ← β, β ← √<q⋅z>,
    ! z̈ ← ż, ż ← z.
    ! ----------------------
    call MatVec(mesh, p, q, env)
    alpha = Dot(mesh, q, p)/(beta**2)
    call Sub(mesh, z, p, z_dot, alpha/beta, 1.0_dp/beta)
    call Sub(mesh, z, z, z_ddot, beta/beta_dot)
    if (present(Precond)) then
      tmp => q_dot; q_dot => q; q => tmp
      call Precond(mesh, q, z, MatVec, env, precond_env)
    else
      q_dot => q; q => z
    end if
    beta_dot = beta; beta = sqrt(Dot(mesh, q, z))
    tmp => z_ddot; z_ddot => z_dot; z_dot => z; z => tmp

    ! ----------------------
    ! Construct and apply reflections:
    ! δ₂ ← cs₁⋅δ₁ + sn₁⋅α, 
    ! γ₁ ← sn₁⋅δ₁ - cs₁⋅α,
    ! ϵ̇ ← ϵ, ϵ ← sn₁⋅β, δ₁ ← -cs₁⋅β,
    !
    ! cs₁, sn₁, γ₂ ← SymOrtho(γ₁, β),
    ! τ ← cs₁⋅ϕ, ϕ ← sn₁⋅ϕ,
    !
    ! cs₂, sn₂, γ̈₆ ← SymOrtho(γ̈₅, ϵ̇),
    ! δ₃ ← sn₂⋅θ̇₁ - cs₂⋅δ₂, γ₃ ← -cs₂⋅γ₂,
    ! η₁ ← sn₂⋅γ₂, θ̇₂ ← cs₂⋅θ̇₁ + sn₂⋅δ₂,
    !
    ! cs₃, sn₃, γ̇₅ ← SymOrtho(γ̇₄, δ₃),
    ! θ₁ ← sn₃⋅γ₃, γ₄ ← -cs₃⋅γ₃.
    ! ----------------------
    delta(2) = cs(1)*delta(1) + sn(1)*alpha
    gamma(1) = sn(1)*delta(1) - cs(1)*alpha
    epsilon_dot = epsilon; epsilon = sn(1)*beta; delta(1) = -cs(1)*beta

    call SymOrtho(gamma(1), beta, cs(1), sn(1), gamma(2))
    tau = cs(1)*phi; phi = sn(1)*phi

    call SymOrtho(gamma_ddot(5), epsilon_dot, cs(2), sn(2), gamma_ddot(6))

    delta(3) = sn(2)*theta_dot(1) - cs(2)*delta(2); gamma(3) = -cs(2)*gamma(2)
    eta = sn(2)*gamma(2); theta_dot(2) = cs(2)*theta_dot(1) + sn(2)*delta(2)

    call SymOrtho(gamma_dot(4), delta(3), cs(3), sn(3), gamma_dot(5))
    theta(1) = sn(3)*gamma(3); gamma(4) = -cs(3)*gamma(3)
    
    ! ----------------------
    ! Update solution:
    ! w₁ ← sn₂ẅ₃ - (cs₂/β̇)q̇,
    ! ẅ₄ ← cs₂ẅ₃ + (sn₂/β̇)q̇,
    !
    ! IF k > 2:
    !  w₂ ← sn₃ẇ₂ - cs₃w₁,
    !  ẇ₃ ← cs₃ẇ₂ + sn₃w₁,
    !  μ̈₃ ← (τ̈ - θ̈₁⋅μ̈̇₃)/γ̈₆,
    ! END IF
    ! IF k > 1:
    !  μ̇₂ ← (τ̇ - θ̇₂⋅μ̈₃ - η̇₁⋅μ̈̇₃)/γ̇₅,
    ! END IF
    ! μ₁ ← γ₂ ≠ 0 ? (τ - θ₁⋅μ̇₂ - η₁⋅μ̈₃)/γ₄ : 0,
    !
    ! ẍ ← ẍ̇ + μ̈₃ẅ₃,
    ! x ← ẍ + μ̇₂ẇ₃,
    ! x ← x + μ₁w₂.
    ! ----------------------
    call Sub(mesh, w, w_ddot, q_dot, cs(2)/beta_dot, sn(2))
    call Add(mesh, w_ddot, w_ddot, q_dot, sn(2)/beta_dot, cs(2))
    if (k > 2) then
      call Set(mesh, p, w)
      call Sub(mesh, w, w_dot, w, cs(3), sn(3))
      call Add(mesh, w_dot, w_dot, p, sn(3), cs(3))
      mu_ddot(3) = ( tau_ddot - theta_ddot(1)*mu_dddot(3) )/gamma_ddot(6)
    end if
    if (k > 1) then
      mu_dot(2) = ( tau_dot - theta_dot(2)*mu_ddot(3) - &
                                & eta_dot*mu_dddot(3) )/gamma_dot(5)
    end if
    !mu(1) = merge(0.0_dp, &
    !  & ( tau - theta(1)*mu_dot(2) - &
    !              & eta*mu_ddot(3) )/gamma(4), gamma(2) == 0.0_dp)
    mu(1) = ( tau - theta(1)*mu_dot(2) - eta*mu_ddot(3) )/gamma(4)
    
    call Add(mesh, x_ddot, x_dddot, w_ddot, mu_ddot(3))
    call Add(mesh, x, x_ddot, w_dot, mu_dot(2))
    call Add(mesh, x, x, w, mu(1))

    ! ----------------------
    ! Shift variables:
    ! k ← k + 1,
    ! γ̈ ← γ̇, γ̇ ← γ,
    ! τ̈ ← τ̇, τ̇ ← τ,
    ! η̇ ← η,
    ! θ̈ ← θ̇, θ̇ ← θ,
    ! μ̈̇ ← μ̈, μ̈ ← μ̇, μ̇ ← μ,
    ! ẅ ← ẇ, ẇ ← w,
    ! ẍ̇ ← ẍ, ẍ ← ẋ, ẋ ← x.
    ! ----------------------
    k = k + 1
    gamma_ddot = gamma_dot; gamma_dot = gamma
    tau_ddot = tau_dot; tau_dot = tau
    eta_dot = eta
    theta_ddot = theta_dot; theta_dot = theta
    mu_dddot = mu_ddot; mu_ddot = mu_dot; mu_dot = mu
    tmp => w_ddot; w_ddot => w_dot; w_dot => w; w => tmp
    tmp => x_dddot; x_dddot => x_ddot; x_ddot => x_dot; x_dot => tmp
    call Set(mesh, x_dot, x)

    ! ----------------------
    ! Check convergence for ϕ and ϕ/ϕ̃.
    ! ( ϕ and ϕ̃ implicitly contain residual norms. )
    ! ----------------------
    if (params%Check(phi, phi/phi_tilde)) exit
  end do
end subroutine Solve_MINRES_QLP$rank
#$end do
#$end if

end module StormRuler_Solvers_MINRES