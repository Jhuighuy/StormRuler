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
use StormRuler_Parameters, only: gMaxIterGMRES

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: Dot, Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_BLAS, only: tMatVecFunc
use StormRuler_Solvers_Precond, only: tPreMatVecFunc

use StormRuler_ConvParams, only: tConvParams

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_MINRES
  module procedure Solve_MINRES
end interface Solve_MINRES

interface Solve_GMRES
  module procedure Solve_GMRES
end interface Solve_GMRES

interface Solve_QMR
  module procedure Solve_QMR
  module procedure Solve_SymmQMR
end interface Solve_QMR

interface Solve_TFQMR
  module procedure Solve_TFQMR
end interface Solve_TFQMR

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
subroutine Solve_MINRES(mesh, x, b, MatVec, params, PreMatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  procedure(tMatVecFunc) :: MatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPreMatVecFunc), optional :: PreMatVec

  real(dp) :: alpha, beta, betaBar, gamma, delta, deltaBar, &
    & epsilon, epsilonBar, tau, phi, phiTilde, cs, sn
  type(tArray) :: tmp, p, q, qBar, w, wBar, wBarBar, z, zBar, zBarBar
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

end subroutine Solve_MINRES

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using 
!! the monstrous Generalized Minimal Residual method (GMRES).
!! 
!! The classical GMRES(𝑚) implementation with restarts
!! after 𝑚 iterations is used.
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
subroutine Solve_GMRES(mesh, x, b, MatVec, params, PreMatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  procedure(tMatVecFunc) :: MatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPreMatVecFunc), optional :: PreMatVec

  logical :: converged
  real(dp) :: chi, phi, phiTilde
  real(dp), pointer :: beta(:), cs(:), sn(:), H(:,:)
  type(tArray) :: Q, s, t, r
  integer(ip) :: i, k

  call AllocArray(r, mold=x)
  associate(m => gMaxIterGMRES)
    call AllocArray(Q, shape=[x%mShape, m + 1])
    allocate(beta(m + 1), cs(m), sn(m), H(m + 1,m))
  end associate

  if (present(PreMatVec)) then
    error stop 'preconditioned GMRES solver is not implemented.'
  end if

  ! ----------------------
  ! Pre-initialize:
  ! 𝜑̃ ← 0,
  ! 𝘤𝘰𝘯𝘷𝘦𝘳𝘨𝘦𝘥 ← 𝙛𝙖𝙡𝙨𝙚.
  ! ----------------------
  phiTilde = 0.0_dp
  converged = .false.

  do while(.not.converged)
    ! ----------------------
    ! Initialize:
    ! 𝒓 ← 𝓐𝒙,
    ! 𝒓 ← 𝒃 - 𝒓,
    ! 𝜑 ← ‖𝒓‖,
    ! 𝗶𝗳 𝜑̃ = 0: // first non-restarted pass.
    !   𝜑̃ ← 𝜑,
    !   Check convergence for 𝜑̃.
    ! 𝗲𝗻𝗱 𝗶𝗳
    ! ----------------------
    call MatVec(mesh, r, x)
    call Sub(mesh, r, b, r)
    phi = Norm_2(mesh, r)
    if (phiTilde <= 0.0) then
      phiTilde = phi
      if (params%Check(phiTilde)) return
    end if

    ! ----------------------
    ! 𝒄𝒔 ← {0}ᵀ, 𝒔𝒏 ← {0}ᵀ,
    ! 𝜷 ← {𝜑,0,…,0}ᵀ,
    ! 𝓠₁ ← 𝒓/𝜑. 
    ! ----------------------
    cs(:) = 0.0_dp; sn(:) = 0.0_dp
    beta(1) = phi; beta(2:) = 0.0_dp
    t = Q%Slice(1); call Scale(mesh, t, r, 1.0_dp/phi)

    do k = 1, gMaxIterGMRES
      ! ----------------------
      ! Arnoldi iteration:
      ! 𝓠ₖ₊₁ ← 𝓐𝓠ₖ,
      ! 𝗳𝗼𝗿 𝑖 = 1, 𝑘 𝗱𝗼:
      !   𝓗ᵢₖ ← <𝓠ₖ₊₁⋅𝓠ᵢ>,
      !   𝓠ₖ₊₁ ← 𝓠ₖ₊₁ - 𝓗ᵢₖ𝓠ᵢ,
      ! 𝗲𝗻𝗱 𝗳𝗼𝗿
      ! 𝓗ₖ₊₁,ₖ ← ‖𝓠ₖ₊₁‖, 𝓠ₖ₊₁ ← 𝓠ₖ₊₁/𝓗ₖ₊₁,ₖ.  
      ! ----------------------
      s = Q%Slice(k); t = Q%Slice(k+1)
      call MatVec(mesh, t, s)
      do i = 1, k
        s = Q%Slice(i); H(i,k) = Dot(mesh, t, s)
        call Sub(mesh, t, t, s, H(i,k))
      end do
      H(k+1,k) = Norm_2(mesh, t); call Scale(mesh, t, t, 1.0_dp/H(k+1,k))

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
      ! ----------------------
      beta(k+1) = -sn(k)*beta(k); beta(k) = cs(k)*beta(k)
      phi = abs(beta(k+1))
      if (params%Check(phi, phi/phiTilde)) exit

    end do

    ! ----------------------
    ! Check if restart is required.
    ! ----------------------
    converged = k <= gMaxIterGMRES
    if (.not.converged) k = gMaxIterGMRES

    ! ----------------------
    ! Compute 𝒙-solution:
    ! 𝜷₁:ₖ ← (𝓗₁:ₖ,₁:ₖ)⁻¹𝜷₁:ₖ, 
    ! 𝗳𝗼𝗿 𝑖 = 1, 𝑘 𝗱𝗼:
    !   𝒙 ← 𝒙 + 𝜷ᵢ𝓠ᵢ.
    ! 𝗲𝗻𝗱 𝗳𝗼𝗿
    ! // Since 𝓗₁:ₖ is upper triangular, 
    ! // operations can be combined:
    ! 𝗳𝗼𝗿 𝑖 = 𝑘, 1, -1 𝗱𝗼:
    !   𝜷ᵢ ← (𝜷ᵢ - <𝓗ᵢ,ᵢ₊₁:ₖ⋅𝜷ᵢ₊₁:ₖ>)/𝓗ᵢᵢ,
    !   𝒙 ← 𝒙 + 𝜷ᵢ𝓠ᵢ.
    ! 𝗲𝗻𝗱 𝗳𝗼𝗿
    ! ----------------------
    do i = k, 1, -1
      beta(i) = beta(i) - dot_product(H(i,(i + 1):k), beta((i + 1):k))
      beta(i) = beta(i)/H(i,i)
      t = Q%Slice(i); call Add(mesh, x, x, t, beta(i))
    end do
  end do

end subroutine Solve_GMRES

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using 
!! the Quasi-Minimal Residual method (QMR).
!! 
!! References:
!! [1] ???
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_QMR(mesh, x, b, &
    & MatVec, ConjMatVec, params, PreMatVec, ConjPreMatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  procedure(tMatVecFunc) :: MatVec, ConjMatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPreMatVecFunc), optional :: PreMatVec, ConjPreMatVec

  error stop 'QMR solver is not implemented.'

end subroutine Solve_QMR

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using 
!! the Quasi-Minimal Residual method (QMR).
!! 
!! References:
!! [1] ???
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_SymmQMR(mesh, x, b, MatVec, params, PreMatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  procedure(tMatVecFunc) :: MatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPreMatVecFunc), optional :: PreMatVec

  if (present(PreMatVec)) then
    call Solve_QMR(mesh, x, b, MatVec, MatVec, params, PreMatVec, PreMatVec)
  else
    call Solve_QMR(mesh, x, b, MatVec, MatVec, params)
  end if

end subroutine Solve_SymmQMR

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using 
!! the Transpose-Free Quasi-Minimal Residual method (TFQMR).
!! 
!! References:
!! [1] Freund, Roland W. 
!!     “Transpose-Free Quasi-Minimal Residual Methods 
!!      for Non-Hermitian Linear Systems.” (1994).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_TFQMR(mesh, x, b, MatVec, params, PreMatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  procedure(tMatVecFunc) :: MatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPreMatVecFunc), optional :: PreMatVec

  error stop 'TFQMR solver is not implemented.'

end subroutine Solve_TFQMR

end module StormRuler_Solvers_MINRES
