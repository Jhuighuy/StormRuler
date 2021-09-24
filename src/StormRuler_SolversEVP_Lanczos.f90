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
module StormRuler_SolversEVP_Lanczos
#$if False

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: &
  & Dot, Norm_2, Fill, Fill_Random, Set, Scale, Add, Sub
#$for type_, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$type_
use StormRuler_Solvers_Precond, only: tPrecondFunc$type_
#$end for
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Tridiag, only: tTridiagMatrix
use StormRuler_Tridiag_LAPACK, only: ComputeEigenpairs_Symm_LAPACK

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute the 𝑚 eigenvalues 𝜆ₖ (and corresponding eigenvectors 𝒙ₖ) 
!! of the self-adjoint operator: 𝓐𝒙ₖ = 𝜆ₖ𝒙ₖ, 𝓐 ∊ ℝⁿ*ⁿ, 
!! using the Lanczos algorithm. 
!! Which eigenvalues are computed:
!! • `SM` (smallest magnitude) case: 
!!   ( 𝑘 ∊ [1,𝑚], 𝜆₁ ≤ … 𝜆ₘ ≤ … ),
!! • `LM` (largest magnitude) case: 
!!   ( 𝑘 ∊ [𝑛-𝑚,𝑛], … ≤ 𝜆ₙ₋ₘ ≤ … 𝜆ₙ ),
!! • `BE` (both ends, 𝑚 = 2𝑝): 
!!   ( 𝑘 ∊ [1,𝑝]∪[𝑛-𝑝,𝑛], 𝜆₁ ≤ … ≤ 𝜆ₚ ≤ … ≤ 𝜆ₙ₋ₚ ≤ … ≤ 𝜆ₙ ).
!!
!! Optionally, an initial guess for the operator equation 𝓐𝒙 = 𝒃 
!! may be refined, using the Lanczos-CG algorithm.
!! In order to enable the refiner, provide the b argument, otherwise
!! x would be used for as a mold for allocations.  
!! Beware that the error for the refinement process is not monitored.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine EigenPairs_Lanczos(mesh, x, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:)
  real(dp), intent(inout) :: x(:,:)
  procedure(tMatVecFuncR) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond !! TODO:
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: k
  real(dp) :: alpha, beta, lambda(2)
  real(dp), pointer :: tmp(:,:), v(:,:), v_dot(:,:), w(:,:)
  real(dp), allocatable :: diag_T(:), subdiag_T(:)

  allocate(v, v_dot, w, mold=x)

  ! ----------------------
  ! Initialize the Lanczos process:
  ! k ← 1,
  ! v̇ ← b,
  ! β ← ‖v̇‖, v̇ ← v̇/β,
  ! w ← Av̇,
  ! α ← <v̇⋅w>,
  ! w ← w - αv̇,
  ! diag(T)ₖ ← α,
  ! subdiag(T)ₖ ← 0,
  ! ----------------------
  k = 1
  call Set(mesh, v_dot, b)
  beta = Norm_2(mesh, v_dot); call Scale(mesh, v_dot, v_dot, 1.0_dp/beta)
  call MatVec(mesh, w, v_dot, env)
  alpha = Dot(mesh, v_dot, w)
  call Sub(mesh, w, w, v_dot, alpha)
  diag_T = [alpha]
  subdiag_T = [0.0_dp]

  do
    ! ----------------------
    ! Continue the Lanczos process:
    ! k ← k + 1,
    ! β ← ‖w‖, v ← w/β,
    ! w ← Av,
    ! α ← <v⋅w>,
    ! w ← w - αv,
    ! w ← w - βv̇,
    ! diag(T)ₖ ← α, subdiag(T)ₖ ← β. 
    ! v̇ ← v,
    ! ----------------------
    k = k + 1
    beta = Norm_2(mesh, w); call Scale(mesh, v, w, 1.0_dp/beta)
    call MatVec(mesh, w, v, env)
    alpha = Dot(mesh, v, w)
    call Sub(mesh, w, w, v, alpha)
    call Sub(mesh, w, w, v_dot, beta)
    diag_T = [diag_T, alpha]
    subdiag_T = [subdiag_T, beta]
    tmp => v_dot; v_dot => v; v => tmp

    ! ----------------------
    ! Compute eigenvalues of T.
    ! ----------------------

    !print *, 'LANCZOS', alpha, beta
    if (k > 5) then
    call TridiagEigenvalues_ContinuantRecursion(k, diag_T, subdiag_T, 7, lambda, 'LE')
    end if
    if (k == 200) then
      !call TridiagEigenvalues_ContinuantRecursion(k, diag_T, subdiag_T, 7, lambda, 'BE')
      error stop 227
    end if
  end do

end subroutine EigenPairs_Lanczos

!! ----------------------------------------------------------------- !!
!! Compute the m eigenvalues of the symmetric tridiagonal matrix: 
!! 𝓣𝒚ⱼ = 𝜃ⱼ𝒚ⱼ, 𝓣 ∊ ℝᵏ*ᵏ, using the continuant recursion algorithm. 
!! ----------------------------------------------------------------- !!
subroutine TridiagEigenvalues_ContinuantRecursion(k, diag_T, subdiag_T, &
    & m, theta, which)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: k, m
  real(dp), intent(in), target :: diag_T(:), subdiag_T(:)
  real(dp), intent(out) :: theta(:)
  character(len=*), intent(in) :: which
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: i, j, l
  real(dp) :: theta_min, theta_max, delta_theta, p(2), d

  block
    !integer(ip) :: info
    !real(dp), allocatable :: lambda(:), zeta(:)
    !lambda = [diag_T]; zeta = [subdiag_T(2:)]
    !call dsterf(k, lambda, zeta, info)
    !print *, 'DSTERF', info, lambda(1), lambda(k)
    type(tTridiagMatrix) :: T
    real(dp), allocatable :: theta(:), y(:,:)
    T%Dim = k; T%Symmetric = .true.
    T%Diag => diag_T; T%Subdiag => subdiag_T
    allocate(theta(2), y(k,2)); theta(:) = 0.0_dp
    call ComputeEigenpairs_Symm_LAPACK(T, 2, 'BE', theta, y)
    print *, 'DSTERF', theta(1), theta(2)
    return
  end block

  ! ----------------------
  ! Let 𝜒(𝜃) = 𝘥𝘦𝘵(𝓣 - 𝜃𝓘) be a characteristic polynomial of 𝓣.
  ! When 𝑘 ≤ 4, 𝜒(𝜃) = 0 may be solved with the direct algorithm. 
  ! ----------------------

  ! ----------------------
  ! 𝑘 = 1, solve directly:
  ! 𝜒(𝜃) = 𝜃 - 𝘥𝘪𝘢𝘨(𝓣)₁ = 0.
  ! ----------------------
  if (k == 1) then
    theta(1) = diag_T(1)
    return
  end if

  ! ----------------------
  ! 𝑘 = 2, solve directly:
  ! 𝜒(𝜃) = 𝜃² - 𝘵𝘳(𝓣)𝜃 + 𝘥𝘦𝘵(𝓣) = 0.
  ! ----------------------
  if (k == 2) then
    associate(tr_T => sum(diag_T), &
      &      det_T => product(diag_T) - subdiag_T(2)**2)
      associate(D => tr_T**2 - 4.0_dp*det_T)
        theta(1) = 0.5_dp*( tr_T - sqrt(D) ) 
        theta(2) = 0.5_dp*( tr_T + sqrt(D) ) 
      end associate
    end associate
    return
  end if

  ! ----------------------
  ! 𝑘 = 3, solve directly (Cardano's formula):
  ! 𝜒(𝜃) = 𝜃³ - 𝘵𝘳(𝓣)𝜃² + ½[𝘵𝘳(𝓣²) - 𝘵𝘳(𝓣)²]𝜃 - 𝘥𝘦𝘵(𝓣) = 0.
  ! ----------------------
  if (k == 3) then
    associate(tr_T => sum(diag_T), &
      &   tr_T_sqr => (sum(diag_T**2) + &
      &                sum(subdiag_T**2) + sum(subdiag_T(2:)**2)), &
      &      det_T => product(diag_T) - subdiag_T(2)**2)
    end associate
    return
  end if

  ! ----------------------
  ! 𝑘 = 4, solve directly (Ferrari's formula):
  ! 𝜒(𝜃) = 𝜃⁴ + ?𝜃³ + ?𝜃² + ?𝜃 + 𝘥𝘦𝘵(𝓣).
  ! ----------------------
  if (k == 4) then
    return
  end if

  ! ----------------------
  ! Find the initial estimates for the extremal
  ! eigenvalues using the Gershgorin circle theorem.
  ! 
  ! Theorem states that (tridiagonal case):
  ! let 𝜃 be any eigenvalue of 𝓣, then
  ! at least one of the following inequalities holds:
  ! |𝜃 - 𝘥𝘪𝘢𝘨(𝓣)₁| ≤ |𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)₂|, or:
  ! |𝜃 - 𝘥𝘪𝘢𝘨(𝓣)ₖ| ≤ |𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)ₖ|, or:
  ! ∃𝑖 ∊ [2,𝑘-1]: |𝜃 - 𝘥𝘪𝘢𝘨(𝓣)ᵢ| ≤ 2|𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)ᵢ|. 
  !
  ! Appying the theorem, we get
  ! the following general eigenvalue bounds:
  ! 𝜃ₘᵢₙ ≥ 𝘮𝘪𝘯{ 𝘥𝘪𝘢𝘨(𝓣)₁ - |𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)₂|, 
  !   𝘥𝘪𝘢𝘨(𝓣)ᵢ - 2|𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)ᵢ|, 𝘥𝘪𝘢𝘨(𝓣)ₖ - |𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)ₖ| },
  ! 𝜃ₘₐₓ ≤ 𝘮𝘢𝘹{ 𝘥𝘪𝘢𝘨(𝓣)₁ + |𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)₂|, 
  !   𝘥𝘪𝘢𝘨(𝓣)ᵢ + 2|𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)ᵢ|, 𝘥𝘪𝘢𝘨(𝓣)ₖ + |𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)ₖ| }, 𝑖 ∊ [2,𝑘-1].
  !
  ! In practice the 𝜃ₘₐₓ is rather accurate, but 𝜃ₘᵢₙ may be 
  ! really underestimated (even may be negative for the case 𝓣 > 0).
  ! ----------------------
  theta_min = diag_T(1) - abs(subdiag_T(2))
  theta_max = diag_T(1) + abs(subdiag_T(2))
  do i = 2, k-1
    theta_min = min(theta_min, diag_T(k) - 2.0_dp*abs(subdiag_T(k)))
    theta_max = max(theta_max, diag_T(k) + 2.0_dp*abs(subdiag_T(k)))
  end do
  theta_min = min(theta_min, diag_T(k) - abs(subdiag_T(k)))
  theta_max = max(theta_max, diag_T(k) + abs(subdiag_T(k)))
  print *, 'THETA_MIN=', theta_min, 'THETA_MAX=', theta_max
  delta_theta = theta_max - theta_min
  error stop 2222

  ! ----------------------
  ! ----------------------
#$if DEBUG
  block
    integer :: unit
    real(dp) :: theta, dtheta
    
    open(unit=unit, file='test/charpr.txt', status='replace')
    theta = theta_min; dtheta = 1e-5_dp*delta_theta
    do while(theta <= theta_max)
      write(unit,*) theta, T_CharPoly(theta)
      theta = theta + dtheta
    end do
    close(unit)
  end block
#$end if

  j = merge(m, 1, which == 'LM')
  do
    print *, j
    ! ----------------------
    ! Proceed to the next eigenvalue.
    ! ----------------------
    select case(which)
      case('SM')
        if (j == m) exit
        theta(j + 1) = theta(j); j = j + 1
      case('LM')
        if (j == 1) exit
        theta(j - 1) = theta(j); j = j - 1
      case('BE')
        if (j == m/2+1) exit
        if (j == m/2) then
          j = m; cycle
        end if
        if (j < m/2) then
          theta(j + 1) = theta(j); j = j + 1
        else
          theta(j - 1) = theta(j); j = j - 1
        end if
    end select
  end do
  error stop 2211

  ! ----------------------
  ! At this step we assume that some eigenvalues are already
  ! computed, and the initial estimates are present for the
  ! ongoing ones.
  !
  ! 
  ! ----------------------

  ! ----------------------
  ! Refine the estimates by applying the Newton solver
  ! to the characteristic polynomial 𝜒(𝜃) = 𝘥𝘦𝘵(𝓣 - 𝜃𝓘),
  ! computing the polynomial using the continuant recursion algorithm.
  !
  ! Assuming θⱼ is  
  ! ----------------------
  do j = 1, m
    do l = 1, 20000
      p = T_CharPoly(theta(j))
      d = p(1)/p(2)
      print *, '***', l, theta(j)
      if (abs(d) <= 1.0e-8_dp) exit
      theta(j) = theta(j) - d
    end do
    error stop 2211
  end do

contains
  function T_CharPoly(theta) result(p)
    ! <<<<<<<<<<<<<<<<<<<<<<
    real(dp), intent(in) :: theta
    real(dp) :: p(2)
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer(ip) :: i
    real(dp) :: mu, p_dot(2), p_ddot(2)

    ! ----------------------
    ! Compute
    ! p₁ = p(θ), p₂ = dp(θ)/dθ, where: p(θ) = ω(θ)⋅det(T - θI), ω(θ) > 0,
    ! using the continuant recursion algorithm.
    !
    ! One should notice that outside of values of the extremal roots,
    ! det(T - θI) explodes to infinity.
    !
    ! ω(θ) is a positive weight function which provides reasonable scaling
    ! for the characteristic polynomial such that |p(θ)| ≲ 1.
    ! This is required to maintain numerical accuracy
    ! ----------------------

    ! ----------------------
    ! With θₘᵢₙ and θₘₐₓ estimated, we can proceed to
    ! the eigenvalues computation of the shifted and rescaled
    ! matrix T̃ = (T - θₘᵢₙI)/Δθ, where Δθ = θₘₐₓ - θₘᵢₙ.
    !
    ! So:
    ! diag(T̃)ᵢ = (diag(T)ᵢ - θₘᵢₙ)/Δθ,
    ! subdiag(T̃)ᵢ = subdiag(T̃)ᵢ/Δθ,
    ! θⱼ = θₘᵢₙ + Δθ⋅eigenvalue(T̃)ⱼ.
    ! ----------------------

    ! ----------------------
    ! p ← (1,0), ṗ ← 0,
    ! ----------------------
    mu = (theta - theta_min)/delta_theta
    p = [1.0_dp, 0.0_dp]; p_dot = 0.0_dp
    do i = 1, k
      ! ----------------------
      ! p̈ ← ṗ, ṗ ← p, 
      ! p ← (diag(T)ᵢ - θ)⋅ṗ - (subdiag(T)ᵢ)²⋅p̈,
      ! p₂ ← p₂ - diag(T)ᵢ⋅ṗ₁.
      ! ----------------------
      p_ddot = p_dot; p_dot = p
      associate(T_diag_tilde => (diag_T(i) - theta_min)/delta_theta, &
        &    T_subdiag_tilde => subdiag_T(i)/delta_theta)
        p = delta_theta*( (T_diag_tilde - mu)*p_dot - (T_subdiag_tilde**2)*p_ddot )
        p(2) = p(2) - T_diag_tilde*p_dot(1)
      end associate
    end do
  end function T_CharPoly
end subroutine TridiagEigenvalues_ContinuantRecursion

#$end if
end module StormRuler_SolversEVP_Lanczos
