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
use StormRuler_Array, only: tArrayR, AllocArray

use StormRuler_BLAS, only: Dot, Norm_2, Fill, Set, Scale, Add, Sub
#$for T, _ in [SCALAR_TYPES[0]]
use StormRuler_BLAS, only: tMatVecFunc$T
use StormRuler_Solvers_Precond, only: tPreMatVecFunc$T
#$end for

use StormRuler_ConvParams, only: tConvParams

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_MINRES
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Solve_MINRES$T
#$end for
end interface Solve_MINRES

interface Solve_GMRES
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Solve_GMRES$T
#$end for
end interface Solve_GMRES

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Generate Givens rotation.
!! ----------------------------------------------------------------- !!
subroutine SymOrtho(a, b, cs, sn, rr)
  real(dp), intent(in) :: a, b
  real(dp), intent(out) :: cs, sn, rr

  ! ----------------------
  ! 𝑟𝑟 ← (𝑎² + 𝑏²)¹ᐟ²,
  ! 𝑐𝑠 ← 𝑎/𝑟𝑟, 𝑠𝑛 ← 𝑏/𝑟𝑟. 
  ! ----------------------
  rr = hypot(a, b)
  if (rr > 0.0_dp) then
    cs = a/rr; sn = b/rr
  else
    cs = 1.0_dp; sn = 0.0_dp
  end if

end subroutine SymOrtho

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
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Solve_MINRES$T(mesh, x, b, MatVec, params, PreMatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: b
  class(tArray$T), intent(inout) :: x
  procedure(tMatVecFunc$T) :: MatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPreMatVecFunc$T), optional :: PreMatVec

  real(dp) :: alpha, beta, betaBar, gamma, delta, deltaBar, &
    & epsilon, epsilonBar, tau, phi, phiTilde, cs, sn
  type(tArray$T) :: tmp, p, q, qBar, w, wBar, wBarBar, z, zBar, zBarBar
  class(*), allocatable :: preEnv

  call AllocArray(p, w, wBar, wBarBar, z, zBar, zBarBar, mold=x)
  if (present(PreMatVec)) call AllocArray(q, qBar, mold=x)

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
  call Fill(mesh, wBar, 0.0_dp)
  call Fill(mesh, wBarBar, 0.0_dp)
  call MatVec(mesh, zBar, x)
  call Sub(mesh, zBar, b, zBar)
  call Fill(mesh, zBarBar, 0.0_dp)
  if (present(PreMatVec)) then
    call PreMatVec(mesh, q, zBar, MatVec, preEnv)
  else
    q = zBar
  end if
  betaBar = 1.0_dp; beta = sqrt(Dot(mesh, q, zBar))
  phi = beta; delta = 0.0_dp; epsilon = 0.0_dp
  cs = -1.0_dp; sn = 0.0_dp

  ! ----------------------
  ! 𝜑̃ ← 𝜑,
  ! Check convergence for 𝜑̃.
  ! ----------------------
  phiTilde = phi
  if (params%Check(phiTilde)) return

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
    call Sub(mesh, z, p, zBar, alpha/beta, 1.0_dp/beta)
    call Sub(mesh, z, z, zBarBar, beta/betaBar)
    if (present(PreMatVec)) then
      tmp = qBar; qBar = q; q = tmp
      call PreMatVec(mesh, q, z, MatVec, preEnv)
    else
      qBar = q; q = z
    end if
    betaBar = beta; beta = sqrt(Dot(mesh, q, z))
    tmp = zBarBar; zBarBar = zBar; zBar = z; z = tmp

    ! ----------------------
    ! Construct and apply rotations:
    ! 𝛿̅ ← 𝑐𝑠⋅𝛿 + 𝑠𝑛⋅𝛼, 𝛾 ← 𝑠𝑛⋅𝛿 - 𝑐𝑠⋅𝛼,
    ! 𝜀̅ ← 𝜀, 𝜀 ← 𝑠𝑛⋅𝛽, 𝛿 ← -𝑐𝑠⋅𝛽,
    ! 𝑐𝑠, 𝑠𝑛, 𝛾 ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝛾, 𝛽),
    ! 𝜏 ← 𝑐𝑠⋅𝜑, 𝜑 ← 𝑠𝑛⋅𝜑.
    ! ----------------------
    deltaBar = cs*delta + sn*alpha; gamma = sn*delta - cs*alpha
    epsilonBar = epsilon; epsilon = sn*beta; delta = -cs*beta
    call SymOrtho(+gamma, beta, cs, sn, gamma)
    tau = cs*phi; phi = sn*phi
    
    ! ----------------------
    ! Update solution:
    ! 𝒘 ← (1/(𝛽̅𝛾))𝒒̅ - (𝛿̅/𝛾)𝒘̅,
    ! 𝒘 ← 𝒘 - (𝜀̅/𝛾)𝒘̿,
    ! 𝒙 ← 𝒙 + 𝜏𝒘,
    ! 𝒘̿ ← 𝒘̅, 𝒘̅ ← 𝒘.
    ! ----------------------
    call Sub(mesh, w, qBar, wBar, deltaBar/gamma, 1.0_dp/(betaBar*gamma))
    call Sub(mesh, w, w, wBarBar, epsilonBar/gamma)
    call Add(mesh, x, x, w, tau)
    tmp = wBarBar; wBarBar = wBar; wBar = w; w = tmp

    ! ----------------------
    ! Check convergence for 𝜑 and 𝜑/𝜑̃.
    ! ( 𝜑 and 𝜑̃ implicitly contain residual norms. )
    ! ----------------------
    if (params%Check(phi, phi/phiTilde)) exit
  end do

end subroutine Solve_MINRES$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using 
!! the monstrous Generalized minimal residual method (GMRES).
!! 
!! GMRES may be applied to the singular problems, and the square
!! least squares problems: ‖(𝓐[𝓜]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓜ᵀ]𝒚, 
!! although convergeance to minimum norm solution is not guaranteed 
!! (is this true?).
!!
!! References:
!! [1] Saad and M.H. Schultz, 
!!     "GMRES: A generalized minimal residual algorithm for solving 
!!      nonsymmetric linear systems", 
!!     SIAM J. Sci. Stat. Comput., 7:856–869, 1986.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Solve_GMRES$T(mesh, x, b, MatVec, params, PreMatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: b
  class(tArray$T), intent(inout) :: x
  procedure(tMatVecFunc$T) :: MatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPreMatVecFunc$T), optional :: PreMatVec

  integer(ip), parameter :: MaxIter = 500

  $typename :: chi, phi, phiTilde
  $typename, pointer :: beta(:), cs(:), sn(:), y(:), H(:,:)
  type(tArray$T) :: Q, Qq, Qi, r
  integer(ip) :: i, k

  call AllocArray(r, mold=x)
  associate(m => MaxIter)
    call AllocArray(Q, shape=[x%mShape, m])
    allocate(beta(m+1), cs(m), sn(m), y(m), H(m+1,m))
  end associate

  ! ----------------------
  ! Pre-initialize:
  ! 𝒓 ← 𝓐𝒙,
  ! 𝒓 ← 𝒃 - 𝒓,
  ! 𝜑̃ ← ‖𝒓‖,
  ! Check convergence for 𝜑̃.
  ! ----------------------

  do
    ! ----------------------
    ! Initialize:
    ! 𝒓 ← 𝓐𝒙,
    ! 𝒓 ← 𝒃 - 𝒓,
    ! 𝜑̃ ← ‖𝒓‖,
    ! Check convergence for 𝜑̃.
    ! ----------------------
    call MatVec(mesh, r, x)
    call Sub(mesh, r, b, r)
    phiTilde = Norm_2(mesh, r)
    if (params%Check(phiTilde)) return

    ! ----------------------
    ! 𝒄𝒔 ← {0}ᵀ, 𝒔𝒏 ← {0}ᵀ,
    ! 𝜷 ← {𝜑̃,0,…,0}ᵀ,
    ! 𝓠₁ ← 𝒓/𝜑̃. 
    ! ----------------------
    cs(:) = 0.0_dp; sn(:) = 0.0_dp
    beta(1) = phiTilde; beta(2:) = 0.0_dp
    Qq = Q%At(1); call Scale(mesh, Qq, r, 1.0_dp/phiTilde)

    do k = 1, MaxIter
      ! ----------------------
      ! Arnoldi iteration:
      ! 𝓠ₖ₊₁ ← 𝓐𝓠ₖ,
      ! 𝗳𝗼𝗿 𝑖 = 1, 𝑘 𝗱𝗼:
      !   𝓗ᵢₖ ← <𝓠ₖ₊₁⋅𝓠ᵢ>,
      !   𝓠ₖ₊₁ ← 𝓠ₖ₊₁ - 𝓗ᵢₖ𝓠ᵢ,
      ! 𝗲𝗻𝗱 𝗳𝗼𝗿
      ! 𝓗ₖ₊₁,ₖ ← ‖𝓠ₖ₊₁‖, 𝓠ₖ₊₁ ← 𝓠ₖ₊₁/𝓗ₖ₊₁,ₖ.  
      ! ----------------------
      Qi = Q%At(k); Qq = Q%At(k+1)
      call MatVec(mesh, Qq, Qi)
      do i = 1, k
        Qi = Q%At(i); H(i,k) = Dot(mesh, Qq, Qi)
        call Sub(mesh, Qq, Qq, Qi, H(i,k))
      end do
      H(k+1,k) = Norm_2(mesh, Qq); call Scale(mesh, Qq, Qq, 1.0_dp/H(k+1,k))

      ! ----------------------
      ! Eliminate the last element in 𝓗
      ! and and update the rotation matrix:
      ! 𝗳𝗼𝗿 𝑖 = 1, 𝑘 - 1 𝗱𝗼:
      !   𝜒 ← 𝒄𝒔ᵢ⋅𝓗ᵢₖ + 𝒔𝒏ᵢ⋅𝓗ᵢ₊₁,ₖ,
      !   𝓗ᵢ₊₁,ₖ ← -𝒔𝒏ᵢ⋅𝓗ᵢₖ + 𝒄𝒔ᵢ⋅𝓗ᵢ₊₁,ₖ 
      !   𝓗ᵢₖ ← 𝜒,
      ! 𝗲𝗻𝗱 𝗳𝗼𝗿
      ! 𝒄𝒔ₖ, 𝒔𝒏ₖ ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝓗ₖₖ, 𝓗ₖ₊₁,ₖ),
      ! 𝓗ₖₖ ← 𝒄𝒔ₖ⋅𝓗ₖₖ + 𝒔𝒏ₖ⋅𝓗ₖ₊₁,ₖ,
      ! 𝓗ₖ₊₁,ₖ ← 0.
      ! ----------------------
      do i = 1, k - 1
        chi = cs(i)*H(i,k) + sn(i)*H(i+1,k)
        H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
        H(i,k) = chi
      end do
      call SymOrtho(H(k,k), H(k+1,k), cs(k), sn(k), chi)
      H(k,k) = cs(k)*H(k,k) + sn(k)*H(k+1,k)
      H(k+1,k) = 0.0_dp

      ! ----------------------
      ! Update the residual vector:
      ! 𝜷ₖ₊₁ ← -𝒔𝒏ₖ𝜷ₖ, 𝜷ₖ ← 𝒄𝒔ₖ⋅𝜷ₖ,
      ! 𝜑 ← |𝜷ₖ₊₁|,
      ! Check convergence for 𝜑 and 𝜑/𝜑̃.
      ! TODO: is 𝜑 = ‖𝒓‖ here?
      ! ----------------------
      beta(k+1) = -sn(k)*beta(k); beta(k) = cs(k)*beta(k)
      phi = abs(beta(k+1))
      if (params%Check(phi, phi/phiTilde)) exit

    end do

    ! ----------------------
    ! Compute 𝒙-solution:
    ! 𝒚 ← (𝓗₁:ₖ,₁:ₖ)⁻¹𝜷₁:ₖ, 
    ! // TODO: here should be ‖𝓗₁:ₖ,₁:ₖ𝒚 - 𝜷₁:ₖ‖₂ → 𝘮𝘪𝘯
    ! 𝗳𝗼𝗿 𝑖 = 1, 𝑘 𝗱𝗼:
    !   𝒙 ← 𝒙 + 𝒚ᵢ𝓠ᵢ.
    ! 𝗲𝗻𝗱 𝗳𝗼𝗿
    ! ----------------------
    do i = 1, k
      Qi = Q%At(i); call Add(mesh, x, x, Qi, y(i))
    end do
    error stop 229

  end do

  error stop 'not implemented'
end subroutine Solve_GMRES$T
#$end for

end module StormRuler_Solvers_MINRES
