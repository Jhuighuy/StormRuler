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
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: Norm_2, Fill, Set, Scale, Add, Sub
#$for type_, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$type_$1
use StormRuler_Solvers_Precond, only: tPrecondFunc$type_
#$end for
use StormRuler_ConvParams, only: tConvParams

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_LSQR
  module procedure Solve_LSQR
  module procedure Solve_LSQR_Symmetric
end interface Solve_LSQR

interface Solve_LSMR
  module procedure Solve_LSMR
  module procedure Solve_LSMR_Symmetric
end interface Solve_LSMR

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a right preconditioned linear least squares problem:
!! 𝘮𝘪𝘯𝘪𝘮𝘪𝘻𝘦{‖𝓐[𝓟]𝒚 - 𝒃‖₂}, 𝒙 = [𝓟]𝒚, using the LSQR method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_LSQR(mesh, x, b, &
    & MatVec, env, MatVec_T, env_T, params, Precond, Precond_T)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:)
  real(dp), intent(inout) :: x(:,:)
  procedure(tMatVecFuncR$1) :: MatVec, MatVec_T
  class(*), intent(inout) :: env, env_T
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond, Precond_T
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
  real(dp), pointer :: s(:,:), t(:,:), &
    & r(:,:), u(:,:), v(:,:), w(:,:), z(:,:)
  class(*), allocatable :: precond_env, precond_env_T
  
  allocate(t, r, u, v, w, z, mold=x)
  if (present(Precond)) then
    allocate(s, mold=x)
  end if

  ! ----------------------
  ! Utilize the initial guess.
  ! Consider the decomposition:
  ! 𝒙 = 𝒙₀ + 𝒛. (*)
  ! Substituting (*) into the equation, we get:
  ! 𝓐[𝓟]𝒚 = 𝒓, where: 𝒛 = [𝓟]𝒚, 𝒓 = 𝒃 - 𝓐𝒙₀.
  ! The last equations can be solved with 𝒚₀ = {0}ᵀ.
  ! ----------------------

  ! ----------------------
  ! Initialize:
  ! 𝒓 ← 𝓐𝒙,
  ! 𝒓 ← 𝒃 - 𝒓,
  ! 𝛽 ← ‖𝒓‖, 𝒖 ← 𝒓/𝛽,
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
  !   𝒔 ← 𝓐ᵀ𝒖, 𝒕 ← 𝓟ᵀ𝒔, 
  ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐ᵀ𝒖, 𝗲𝗻𝗱 𝗶𝗳
  ! 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/α.
  ! ----------------------
  call MatVec(mesh, r, x, env)
  call Sub(mesh, r, b, r)
  beta = Norm_2(mesh, r); call Scale(mesh, u, r, 1.0_dp/beta)
  if (present(Precond)) then
    call MatVec_T(mesh, s, u, env_T)
    call Precond_T(mesh, t, s, MatVec_T, env_T, precond_env_T)
  else
    call MatVec_T(mesh, t, u, env_T)
  end if
  alpha = Norm_2(mesh, t); call Scale(mesh, v, t, 1.0_dp/alpha)

  ! ----------------------
  ! 𝜑̅ ← 𝛽, 𝜌̅ ← 𝛼.
  ! 𝒛 ← {0}ᵀ,
  ! 𝒘 ← 𝒗,
  ! ----------------------
  phi_bar = beta; rho_bar = alpha
  call Fill(mesh, z, 0.0_dp)
  call Set(mesh, w, v)

  ! ----------------------
  ! 𝜑̃ ← 𝜑̅,
  ! Check convergence for 𝜑̃.
  ! ----------------------
  phi_tilde = phi_bar
  if (params%Check(phi_tilde)) return
  
  do
    ! ----------------------
    ! Continue the bidiagonalization:
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
    !   𝒔 ← 𝓟𝒗, 𝒕 ← 𝓐𝒔,
    ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐𝒗, 𝗲𝗻𝗱 𝗶𝗳
    ! 𝒕 ← 𝒕 - 𝛼𝒖,
    ! 𝛽 ← ‖𝒕‖, 𝒖 ← 𝒕/𝛽,
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    !   𝒔 ← 𝓐ᵀ𝒖, 𝒕 ← 𝓟ᵀ𝒔, 
    ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐ᵀ𝒖, 𝗲𝗻𝗱 𝗶𝗳
    ! 𝒕 ← 𝒕 - 𝛽𝒗,
    ! 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
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
      call MatVec_T(mesh, s, u, env_T)
      call Precond_T(mesh, t, s, MatVec_T, env_T, precond_env_T)
    else
      call MatVec_T(mesh, t, u, env_T)
    end if
    call Sub(mesh, t, t, v, beta)
    alpha = Norm_2(mesh, t); call Scale(mesh, v, t, 1.0_dp/alpha)
    
    ! ----------------------
    ! Construct and apply rotation:
    ! 𝜌 ← (𝜌̅² + 𝛽²)¹ᐟ²,
    ! 𝑐𝑠 ← 𝜌̅/𝜌, 𝑠𝑛 ← 𝛽/𝜌,
    ! 𝜃 ← 𝑠𝑛⋅𝛼, 𝜌̅ ← -𝑐𝑠⋅𝛼,
    ! 𝜑 ← 𝑐𝑠⋅𝜑, 𝜑̅ ← 𝑠𝑛⋅𝜑̅.
    ! ----------------------
    rho = hypot(rho_bar, beta)
    cs = rho_bar/rho; sn = beta/rho
    theta = sn*alpha; rho_bar = -cs*alpha
    phi = cs*phi_bar; phi_bar = sn*phi_bar

    ! ----------------------
    ! Update 𝒛-solution:
    ! 𝒛 ← 𝒛 + (𝜑/𝜌)𝒘,
    ! 𝒘 ← 𝒗 - (𝜃/𝜌)𝒘.
    ! Check convergence for 𝜑̅ and 𝜑̅/𝜑̃.
    ! ( 𝜑̅ and 𝜑̃ implicitly contain residual norms;
    !   𝜑̅|𝜌̅| implicitly contain (𝓐[𝓟])ᵀ-residual norms, ‖(𝓐[𝓟])ᵀ𝒓‖. )
    ! ----------------------
    call Add(mesh, z, z, w, phi/rho)
    call Sub(mesh, w, v, w, theta/rho)
    if (params%Check(phi_bar, phi_bar/phi_tilde)) exit
  end do

  ! ----------------------
  ! Compute 𝒙-solution:
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  !   𝒕 ← 𝓟𝒛, 𝒙 ← 𝒙 + 𝒕.
  ! 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒙 + 𝒛. 𝗲𝗻𝗱 𝗶𝗳
  ! ----------------------
  if (present(Precond)) then
    call Precond(mesh, t, z, MatVec, env, precond_env)
    call Add(mesh, x, x, t)
  else
    call Add(mesh, x, x, z)
  end if

end subroutine Solve_LSQR
subroutine Solve_LSQR_Symmetric(mesh, x, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:)
  real(dp), intent(inout) :: x(:,:)
  procedure(tMatVecFuncR$1) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>

  ! ----------------------
  ! Using LSMR in symmetric case is not recommended,
  ! use MINRES instead.
  ! ----------------------

  if (present(Precond)) then
    call Solve_LSQR(mesh, x, b, &
      & MatVec, env, MatVec, env, params, Precond, Precond)
  else
    call Solve_LSQR(mesh, x, b, MatVec, env, MatVec, env, params)
  end if
  
end subroutine Solve_LSQR_Symmetric

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a right preconditioned linear least squares problem:
!! 𝘮𝘪𝘯𝘪𝘮𝘪𝘻𝘦{‖𝓐[𝓟]𝒚 - 𝒃‖₂}, 𝒙 = [𝓟]𝒚, using the LSMR method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_LSMR(mesh, x, b, &
    & MatVec, env, MatVec_T, env_T, params, Precond, Precond_T)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:)
  real(dp), intent(inout) :: x(:,:)
  procedure(tMatVecFuncR$1) :: MatVec, MatVec_T
  class(*), intent(inout) :: env, env_T
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond, Precond_T
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  ! ----------------------
  ! [1] Fong, D. C. and M. Saunders. 
  !     “LSMR: An Iterative Algorithm for Sparse Least-Squares Problems.” 
  !     SIAM J. Sci. Comput. 33 (2011): 2950-2971.
  ! [2] Karimi, S., D. K. Salkuyeh and F. Toutounian. 
  !     “A preconditioner for the LSQR algorithm.” 
  !     Journal of applied mathematics & informatics 26 (2008): 213-222.
  ! ----------------------

  real(dp) :: alpha, alpha_bar, beta, rho, rho_bar, &
    & theta, theta_bar, psi, psi_bar, psi_tilde, zeta, &
    & cs, sn, cs_bar, sn_bar
  real(dp), pointer :: r(:,:), s(:,:), t(:,:), &
    & w(:,:), h(:,:), u(:,:), v(:,:), z(:,:)
  class(*), allocatable :: precond_env, precond_env_T
  
  allocate(t, r, u, v, w, h, z, mold=x)
  if (present(Precond)) then
    allocate(s, mold=x)
  end if

  ! ----------------------
  ! Utilize the initial guess.
  ! Consider the decomposition:
  ! 𝒙 = 𝒙₀ + 𝒛. (*)
  ! Substituting (*) into the equation, we get:
  ! 𝓐[𝓟]𝒚 = 𝒓, where: 𝒛 = [𝓟]𝒚, 𝒓 = 𝒃 - 𝓐𝒙₀.
  ! The last equations can be solved with 𝒚₀ = {0}ᵀ.
  ! ----------------------

  ! ----------------------
  ! Initialize:
  ! 𝒓 ← 𝓐𝒙,
  ! 𝒓 ← 𝒃 - 𝒓,
  ! 𝛽 ← ‖𝒓‖, 𝒖 ← 𝒓/𝛽,
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
  !   𝒔 ← 𝓐ᵀ𝒖, 𝒕 ← 𝓟ᵀ𝒔, 
  ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐ᵀ𝒖, 𝗲𝗻𝗱 𝗶𝗳
  ! 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/α.
  ! ----------------------
  call MatVec(mesh, r, x, env)
  call Sub(mesh, r, b, r)
  beta = Norm_2(mesh, r); call Scale(mesh, u, r, 1.0_dp/beta)
  if (present(Precond)) then
    call MatVec_T(mesh, s, u, env_T)
    call Precond_T(mesh, t, s, MatVec_T, env_T, precond_env_T)
  else
    call MatVec_T(mesh, t, u, env_T)
  end if
  alpha = Norm_2(mesh, t); call Scale(mesh, v, t, 1.0_dp/alpha)

  ! ----------------------
  ! 𝛼̅ ← 𝛼, 𝜓̅ ← 𝛼𝛽,
  ! 𝜁 ← 1, 𝑐̅𝑠̅ ← 1, 𝑠̅𝑛̅ ← 0,
  ! 𝒛 ← {0}ᵀ,
  ! 𝒘 ← 𝒗, 𝒉 ← {0}ᵀ.
  ! ----------------------
  alpha_bar = alpha; psi_bar = alpha*beta
  zeta = 1.0_dp; cs_bar = 1.0_dp; sn_bar = 0.0_dp
  call Fill(mesh, z, 0.0_dp)
  call Set(mesh, w, v); call Fill(mesh, h, 0.0_dp)

  ! ----------------------
  ! 𝜓̃ ← 𝜓̅,
  ! Check convergence for 𝜓̃.
  ! ----------------------
  psi_tilde = psi_bar
  if (params%Check(psi_tilde)) return
  
  do
    ! ----------------------
    ! Continue the bidiagonalization:
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
    !   𝒔 ← 𝓟𝒗, 𝒕 ← 𝓐𝒔,
    ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐𝒗, 𝗲𝗻𝗱 𝗶𝗳
    ! 𝒕 ← 𝒕 - 𝛼𝒖,
    ! 𝛽 ← ‖𝒕‖, 𝒖 ← 𝒕/𝛽,
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    !   𝒔 ← 𝓐ᵀ𝒖, 𝒕 ← 𝓟ᵀ𝒔, 
    ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐ᵀ𝒖, 𝗲𝗻𝗱 𝗶𝗳
    ! 𝒕 ← 𝒕 - 𝛽𝒗,
    ! 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
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
      call MatVec_T(mesh, s, u, env_T)
      call Precond_T(mesh, t, s, MatVec_T, env_T, precond_env_T)
    else
      call MatVec_T(mesh, t, u, env_T)
    end if
    call Sub(mesh, t, t, v, beta)
    alpha = Norm_2(mesh, t); call Scale(mesh, v, t, 1.0_dp/alpha)
    
    ! ----------------------
    ! Construct and apply rotations:
    ! 𝜌 ← (𝛼̅² + 𝛽²)¹ᐟ²,
    ! 𝑐𝑠 ← 𝛼̅/𝜌, 𝑠𝑛 ← 𝛽/𝜌,
    ! 𝜃 ← 𝑠𝑛⋅𝛼, 𝛼̅ ← 𝑐𝑠⋅𝛼,
    ! 𝜃̅ ← 𝑠̅𝑛̅⋅𝜌, 𝜌̅ ← [(𝑐̅𝑠̅⋅𝜌)² + 𝜃²]¹ᐟ²,
    ! 𝑐̅𝑠̅ ← 𝑐̅𝑠̅⋅𝜌/𝜌̅, 𝑠̅𝑛̅ ← 𝜃/𝜌̅,
    ! 𝜓 ← 𝑐̅𝑠̅⋅𝜓̅, 𝜓̅ ← -𝑠̅𝑛̅⋅𝜓̅.
    ! ----------------------
    rho = hypot(alpha_bar, beta)
    cs = alpha_bar/rho; sn = beta/rho
    theta = sn*alpha; alpha_bar = cs*alpha
    theta_bar = sn_bar*rho; rho_bar = hypot(cs_bar*rho, theta)
    cs_bar = cs_bar*rho/rho_bar; sn_bar = theta/rho_bar
    psi = cs_bar*psi_bar; psi_bar = -sn_bar*psi_bar

    ! ----------------------
    ! Update 𝒛-solution:
    ! 𝒉 ← 𝒘 - (𝜃𝜌/𝜁)𝒉,
    ! 𝜁 ← 𝜌𝜌̅,
    ! 𝒛 ← 𝒛 + (𝜓/𝜁)𝒉,
    ! 𝒘 ← 𝒗 - (𝜃/𝜌)𝒘.
    ! Check convergence for |𝜓̅| and |𝜓̅/𝜓̃|.
    ! ( |𝜓̅| and |𝜓̃| implicitly contain (𝓐[𝓟])ᵀ-residual norms, ‖(𝓐[𝓟])ᵀ𝒓‖. )
    ! ----------------------
    call Sub(mesh, h, w, h, theta_bar*rho/zeta)
    zeta = rho*rho_bar
    call Add(mesh, z, z, h, psi/zeta)
    call Sub(mesh, w, v, w, theta/rho)
    if (params%Check(abs(psi_bar), abs(psi_bar/psi_tilde))) exit
  end do

  ! ----------------------
  ! Compute 𝒙-solution:
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  !   𝒕 ← 𝓟𝒛, 𝒙 ← 𝒙 + 𝒕.
  ! 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒙 + 𝒛. 𝗲𝗻𝗱 𝗶𝗳
  ! ----------------------
  if (present(Precond)) then
    call Precond(mesh, t, z, MatVec, env, precond_env)
    call Add(mesh, x, x, t)
  else
    call Add(mesh, x, x, z)
  end if

end subroutine Solve_LSMR
subroutine Solve_LSMR_Symmetric(mesh, x, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:)
  real(dp), intent(inout) :: x(:,:)
  procedure(tMatVecFuncR$1) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>

  ! ----------------------
  ! Using LSMR in symmetric case is not recommended,
  ! use MINRES instead.
  ! ----------------------

  if (present(Precond)) then
    call Solve_LSMR(mesh, x, b, &
      & MatVec, env, MatVec, env, params, Precond, Precond)
  else
    call Solve_LSMR(mesh, x, b, MatVec, env, MatVec, env, params)
  end if
  
end subroutine Solve_LSMR_Symmetric

end module StormRuler_Solvers_LSQR
