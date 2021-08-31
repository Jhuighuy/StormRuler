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

use StormRuler_Parameters, only: dp
use StormRuler_Helpers, only: SafeDivide
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@, &
  & Dot, Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_Base, only: @{tPrecondFunc$$@|@0, NUM_RANKS}@
use StormRuler_Solvers_Base, only: ResidualNorm

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_SYMMLQ
#$do rank = 0, NUM_RANKS
  module procedure Solve_SYMMLQ$rank
#$end do
end interface Solve_SYMMLQ

interface Solve_MINRES
#$do rank = 0, NUM_RANKS
  module procedure Solve_MINRES$rank
#$end do
end interface Solve_MINRES

interface Solve_GMRES
#$do rank = 0, NUM_RANKS
  module procedure Solve_GMRES$rank
#$end do
end interface Solve_GMRES

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!!
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_SYMMLQ$rank(mesh, u, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFunc$rank), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>
end subroutine Solve_SYMMLQ$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint indefinite operator equation: 
!! [P¹ᐟ²]A[P¹ᐟ²]y = [P¹ᐟ²]b, [P¹ᐟ²]y = x, using the MINRES method.
!! ( P = Pᵀ > 0 is required.  )
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

  real(dp) :: alpha, beta, beta_bar, gamma, gamma_hat, &
    & delta, delta_hat, epsilon, epsilon_bar, tau, phi, phi_tilde, cs, sn
  real(dp), pointer :: t(@:,:), p(@:,:), q(@:,:), q_bar(@:,:), &
    & d(@:,:), d_bar(@:,:), d_bar_bar(@:,:), z(@:,:), z_bar(@:,:), z_bar_bar(@:,:)
  class(*), allocatable :: precond_env

  allocate(p, d, d_bar, d_bar_bar, z, z_bar, z_bar_bar, mold=x)
  if (present(Precond)) then
    allocate(q, q_bar, mold=x)
  end if

  ! ----------------------
  ! Initialize:
  ! d̿ ← {0}ᵀ,
  ! d̅ ← {0}ᵀ,
  ! z̿ ← 0,
  ! z̅ ← Ax,
  ! z̅ ← b - z̅,
  ! q ← Pz̅, OR: q ← z̅,
  ! β̅ ← 1, β ← √<q⋅z̅>,
  ! ϕ ← β, ψ̃ ← 0, δ ← 0, ϵ ← 0,
  ! cs ← -1, sn ← 0.
  ! ----------------------
  call Fill(mesh, d_bar_bar, 0.0_dp)
  call Fill(mesh, d_bar, 0.0_dp)
  call Fill(mesh, z_bar_bar, 0.0_dp)
  call MatVec(mesh, z_bar, x, env)
  call Sub(mesh, z_bar, b, z_bar)
  if (present(Precond)) then
    call Precond(mesh, q, z_bar, MatVec, env, precond_env)
  else
    q => z_bar
  end if
  beta_bar = 1.0_dp; beta = sqrt(Dot(mesh, q, z_bar))
  phi = beta; delta = 0.0_dp; epsilon = 0.0_dp
  cs = -1.0_dp; sn = 0.0_dp

  ! ----------------------
  ! ϕ̃ ← ϕ,
  ! Check convergence for ϕ.
  ! ----------------------
  phi_tilde = phi
  if (params%Check(phi_tilde)) return

  do
    ! ----------------------
    ! Continue the Lanczos process:
    ! p ← Aq,
    ! α ← <q⋅p>/β²,
    ! z ← (1/β)p - (α/β)z̅,
    ! z ← z - (β/β̅)z̿,
    ! q̅ ← q,
    ! q ← Pz, OR: q ← z,
    ! β̅ ← β, β ← √<q⋅z>,
    ! z̿ ← z̅, z̅ ← z.
    ! ----------------------
    call MatVec(mesh, p, q, env)
    alpha = Dot(mesh, q, p)/(beta**2)
    call Sub(mesh, z, p, z_bar, alpha/beta, 1.0_dp/beta)
    call Sub(mesh, z, z, z_bar_bar, beta/beta_bar)
    if (present(Precond)) then
      t => q_bar; q_bar => q; q => t
      call Precond(mesh, q, z, MatVec, env, precond_env)
    else
      q_bar => q; q => z
    end if
    beta_bar = beta; beta = sqrt(Dot(mesh, q, z))
    t => z_bar_bar; z_bar_bar => z_bar; z_bar => z; z => t

    ! ----------------------
    ! Construct and apply rotations:
    ! δ̂ ← cs⋅δ + sn⋅α, γ ← sn⋅δ - cs⋅α,
    ! ϵ̅ ← ϵ, ϵ ← sn⋅β, δ ← -cs⋅β,
    ! γ̂ ← √(γ² + β²),
    ! cs ← γ/γ̂, sn ← β/γ̂,
    ! τ ← cs⋅ϕ, ϕ ← sn⋅ϕ.
    ! ----------------------
    delta_hat = cs*delta + sn*alpha; gamma = sn*delta - cs*alpha
    epsilon_bar = epsilon; epsilon = sn*beta; delta = -cs*beta
    gamma_hat = hypot(gamma, beta)
    cs = gamma/gamma_hat; sn = beta/gamma_hat
    tau = cs*phi; phi = sn*phi
    
    ! ----------------------
    ! Update solution:
    ! d ← (1/(β̅γ̂))q̅ - (δ̂/γ̂)d̅,
    ! d ← d - (ϵ̅/γ̂)d̿,
    ! x ← x + τd,
    ! d̿ ← d̅, d̅ ← d.
    ! ----------------------
    call Sub(mesh, d, q_bar, d_bar, delta_hat/gamma_hat, 1.0_dp/(beta_bar*gamma_hat))
    call Sub(mesh, d, d, d_bar_bar, epsilon_bar/gamma_hat)
    call Add(mesh, x, x, d, tau)
    t => d_bar_bar; d_bar_bar => d_bar; d_bar => d; d => t

    ! ----------------------
    ! Check convergence for ϕ and ϕ/ϕ̃.
    ! ( ϕ and ϕ̃ implicitly contain residual norms. )
    ! ----------------------
    if (params%Check(phi, phi/phi_tilde)) exit
  end do
end subroutine Solve_MINRES$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!!
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_GMRES$rank(mesh, u, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFunc$rank), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>
end subroutine Solve_GMRES$rank
#$end do

end module StormRuler_Solvers_MINRES
