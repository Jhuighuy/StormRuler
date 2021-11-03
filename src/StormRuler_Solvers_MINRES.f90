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
use StormRuler_Array, only: tArrayR, AllocArrayMold

use StormRuler_BLAS, only: Dot, Norm_2, Fill, Set, Scale, Add, Sub
#$for T, _ in [SCALAR_TYPES[0]]
use StormRuler_BLAS, only: tMatVecFunc$T
use StormRuler_Solvers_Precond, only: tPrecondFunc$T
#$end for

use StormRuler_ConvParams, only: tConvParams

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint indefinite operator equation: 
!! [𝓜]𝓐[𝓜ᵀ]𝒚 = [𝓜]𝒃, [𝓜ᵀ]𝒚 = 𝒙, [𝓜𝓜ᵀ = 𝓟], using the MINRES method.
!!
!! Despite 𝓐 may be indefinite, a positive-definite preconditioner 𝓟 
!! is explicitly required.
!!
!! MINRES may be applied to the singular problems, and the self-adjoint
!! least squares problems: ‖[𝓜](𝓐[𝓜ᵀ]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓜ᵀ]𝒚, 
!! although convergeance to minimum norm solution is not guaranteed.
!! 
!! References:
!! [1] Paige, C. and M. Saunders. 
!!     “Solution of Sparse Indefinite Systems of Linear Equations.” 
!!     SIAM Journal on Numerical Analysis 12 (1975): 617-629.
!! [2] Choi, S.-C. T.
!!     “Iterative Methods for Singular Linear Equations and 
!!     Least-Squares Problems” PhD thesis, ICME, Stanford University.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_MINRES(mesh, x, b, MatVec, params, Precond)
  class(tMesh), intent(inout) :: mesh
  class(tArrayR), intent(in) :: b
  class(tArrayR), intent(inout) :: x
  procedure(tMatVecFuncR) :: MatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond

  real(dp) :: alpha, beta, beta_bar, gamma, &
    & delta, delta_bar, epsilon, epsilon_bar, &
    & tau, phi, phi_tilde, cs, sn
  type(tArrayR) :: tmp, p, q, q_bar, w, w_bar, w_bbar, z, z_bar, z_bbar
  class(*), allocatable :: precond_env

  call AllocArrayMold(p, w, w_bar, w_bbar, z, z_bar, z_bbar, mold=x)
  if (present(Precond)) call AllocArrayMold(q, q_bar, mold=x)

  ! ----------------------
  ! Initialize:
  ! 𝒘̅ ← {0}ᵀ,
  ! 𝒘̿ ← {0}ᵀ,
  ! 𝒛̅ ← 𝓐𝒙,     // Modification in order to
  ! 𝒛̅ ← 𝒃 - 𝒛̅,  // utilize the initial guess.
  ! 𝒛̿ ← {0}ᵀ,
  ! 𝒒 ← [𝓟]𝒛̅,
  ! 𝛽̅ ← 1, 𝛽 ← √<𝒒⋅𝒛̅>,
  ! 𝜑 ← 𝛽, 𝛿 ← 0, 𝜀 ← 0,
  ! 𝑐𝑠 ← -1, 𝑠𝑛 ← 0.
  ! ----------------------
  call Fill(mesh, w_bar, 0.0_dp)
  call Fill(mesh, w_bbar, 0.0_dp)
  call MatVec(mesh, z_bar, x)
  call Sub(mesh, z_bar, b, z_bar)
  call Fill(mesh, z_bbar, 0.0_dp)
  if (present(Precond)) then
    call Precond(mesh, q, z_bar, MatVec, precond_env)
  else
    q = z_bar
  end if
  beta_bar = 1.0_dp; beta = sqrt(Dot(mesh, q, z_bar))
  phi = beta; delta = 0.0_dp; epsilon = 0.0_dp
  cs = -1.0_dp; sn = 0.0_dp

  ! ----------------------
  ! 𝜑̃ ← 𝜑,
  ! Check convergence for 𝜑̃.
  ! ----------------------
  phi_tilde = phi
  if (params%Check(phi_tilde)) return

  do
    ! ----------------------
    ! Continue the Lanczos process:
    ! 𝒑 ← 𝓐𝒒,
    ! 𝛼 ← <𝒒⋅𝒑>/𝛽²,
    ! 𝒛 ← (1/𝛽)𝒑 - (𝛼/𝛽)𝒛̅,
    ! 𝒛 ← 𝒛 - (𝛽/𝛽̅)𝒛̿,
    ! 𝒒̅ ← 𝒒, 𝒒 ← [𝓟]𝒛,
    ! 𝛽̅ ← 𝛽, 𝛽 ← √<𝒒⋅𝒛>,
    ! 𝒛̿ ← 𝒛̅, 𝒛̅ ← 𝒛.
    ! ----------------------
    call MatVec(mesh, p, q)
    alpha = Dot(mesh, q, p)/(beta**2)
    call Sub(mesh, z, p, z_bar, alpha/beta, 1.0_dp/beta)
    call Sub(mesh, z, z, z_bbar, beta/beta_bar)
    if (present(Precond)) then
      tmp = q_bar; q_bar = q; q = tmp
      call Precond(mesh, q, z, MatVec, precond_env)
    else
      q_bar = q; q = z
    end if
    beta_bar = beta; beta = sqrt(Dot(mesh, q, z))
    tmp = z_bbar; z_bbar = z_bar; z_bar = z; z = tmp

    ! ----------------------
    ! Construct and apply rotations:
    ! 𝛿̅ ← 𝑐𝑠⋅𝛿 + 𝑠𝑛⋅𝛼, 𝛾 ← 𝑠𝑛⋅𝛿 - 𝑐𝑠⋅𝛼,
    ! 𝜀̅ ← 𝜀, 𝜀 ← 𝑠𝑛⋅𝛽, 𝛿 ← -𝑐𝑠⋅𝛽,
    ! 𝑐𝑠, 𝑠𝑛, 𝛾 ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝛾, 𝛽),
    ! 𝜏 ← 𝑐𝑠⋅𝜑, 𝜑 ← 𝑠𝑛⋅𝜑.
    ! ----------------------
    delta_bar = cs*delta + sn*alpha; gamma = sn*delta - cs*alpha
    epsilon_bar = epsilon; epsilon = sn*beta; delta = -cs*beta
    call SymOrtho(gamma, beta, cs, sn, gamma)
    tau = cs*phi; phi = sn*phi
    
    ! ----------------------
    ! Update solution:
    ! 𝒘 ← (1/(𝛽̅𝛾))𝒒̅ - (𝛿̅/𝛾)𝒘̅,
    ! 𝒘 ← 𝒘 - (𝜀̅/𝛾)𝒘̿,
    ! 𝒙 ← 𝒙 + 𝜏𝒘,
    ! 𝒘̿ ← 𝒘̅, 𝒘̅ ← 𝒘.
    ! ----------------------
    call Sub(mesh, w, q_bar, w_bar, delta_bar/gamma, 1.0_dp/(beta_bar*gamma))
    call Sub(mesh, w, w, w_bbar, epsilon_bar/gamma)
    call Add(mesh, x, x, w, tau)
    tmp = w_bbar; w_bbar = w_bar; w_bar = w; w = tmp

    ! ----------------------
    ! Check convergence for 𝜑 and 𝜑/𝜑̃.
    ! ( 𝜑 and 𝜑̃ implicitly contain residual norms. )
    ! ----------------------
    if (params%Check(phi, phi/phi_tilde)) exit
  end do
  
contains
  subroutine SymOrtho(a, b, cs, sn, rr)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: cs, sn, rr

    rr = hypot(a, b)
    if (rr > 0.0_dp) then
      cs = a/rr; sn = b/rr
    else
      cs = 1.0_dp; sn = 0.0_dp
    end if
  end subroutine SymOrtho
end subroutine Solve_MINRES

end module StormRuler_Solvers_MINRES
