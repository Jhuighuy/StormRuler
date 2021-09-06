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
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@, &
  & Dot, Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_Precond, only: @{tPrecondFunc$$@|@0, NUM_RANKS}@

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_MINRES
#$do rank = 0, NUM_RANKS
  module procedure Solve_MINRES$rank
#$end do
end interface Solve_MINRES

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

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
  ! ż ← Ax,     // Modification in order to
  ! ż ← b - ż,  // utilize the initial guess.
  ! z̈ ← 0,
  ! IF P THEN:
  !   q ← Pż, ELSE: q ← ż, END IF
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
    ! IF (P ≠ NONE) THEN: 
    !   q ← Pz, ELSE: q ← z, END IF
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
end subroutine Solve_MINRES$rank
#$end do

end module StormRuler_Solvers_MINRES
