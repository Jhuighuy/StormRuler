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
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_ConvParams, only: tConvParams
use StormRuler_Precond, only: tPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_LSQR
  module procedure Solve_LSQR
  module procedure Solve_SymmLSQR
end interface Solve_LSQR

interface Solve_LSMR
  module procedure Solve_LSMR
  module procedure Solve_SymmLSMR
end interface Solve_LSMR

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a right preconditioned linear least squares problem:
!! ‖𝓐[𝓟]𝒚 - 𝒃‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚, using the LSQR method.
!!
!! LSQR is algebraically equivalent to applying CG
!! to the normal equations: (𝓐[𝓟])*𝓐[𝓟]𝒚 = (𝓐[𝓟])*𝒃, 𝒙 = [𝓟]𝒚,
!! (or, equivalently, [𝓟*]𝓐*𝓐[𝓟]𝒚 = [𝓟*]𝓐*𝒃, 𝒙 = [𝓟]𝒚),
!! but but has better numerical properties.
!!
!! The residual norm ‖𝓐[𝓟]𝒚 - 𝒃‖₂ decreases monotonically, 
!! while the normal equation's residual norm ‖(𝓐[𝓟])*(𝓐[𝓟]𝒚 - 𝒃)‖ 
!! is not guaranteed to decrease.
!!
!! References:
!! [1] Paige, C. and M. Saunders. 
!!     “LSQR: An Algorithm for Sparse Linear Equations and 
!!     Sparse Least Squares.” ACM Trans. Math. Softw. 8 (1982): 43-71.
!! [2] Karimi, S., D. K. Salkuyeh and F. Toutounian. 
!!     “A preconditioner for the LSQR algorithm.” 
!!     Journal of applied mathematics & informatics 26 (2008): 213-222.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_LSQR(mesh, x, b, MatVec, &
    & ConjMatVec, params, precond, conjPrecond)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  class(tConvParams), intent(inout) :: params
  class(tPreconditioner), intent(inout), optional :: precond, conjPrecond
  procedure(tMatVecFunc) :: MatVec, ConjMatVec
  
  real(dp) :: alpha, beta, rho, rhoBar, theta, phi, phiBar, phiTilde, cs, sn
  type(tArray) :: s, t, r, u, v, w, z
  
  call AllocArray(t, r, u, v, w, z, mold=x)
  if (present(precond)) then
    call AllocArray(s, mold=x)
    call precond%Init(mesh, MatVec)
    call conjPrecond%Init(mesh, MatVec)
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
  !   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔, 
  ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
  ! 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
  ! ----------------------
  call MatVec(mesh, r, x)
  call Sub(mesh, r, b, r)
  beta = Norm_2(mesh, r); call Scale(mesh, u, r, 1.0_dp/beta)
  if (present(precond)) then
    call ConjMatVec(mesh, s, u)
    call conjPrecond%Apply(mesh, t, s, ConjMatVec)
  else
    call ConjMatVec(mesh, t, u)
  end if
  alpha = Norm_2(mesh, t); call Scale(mesh, v, t, 1.0_dp/alpha)

  ! ----------------------
  ! 𝜑̅ ← 𝛽, 𝜌̅ ← 𝛼.
  ! 𝒛 ← {0}ᵀ,
  ! 𝒘 ← 𝒗,
  ! ----------------------
  phiBar = beta; rhoBar = alpha
  call Fill(mesh, z, 0.0_dp)
  call Set(mesh, w, v)

  ! ----------------------
  ! 𝜑̃ ← 𝜑̅,
  ! Check convergence for 𝜑̃.
  ! ----------------------
  phiTilde = phiBar
  if (params%Check(phiTilde)) return
  
  do
    ! ----------------------
    ! Continue the bidiagonalization:
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
    !   𝒔 ← 𝓟𝒗, 𝒕 ← 𝓐𝒔,
    ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐𝒗, 𝗲𝗻𝗱 𝗶𝗳
    ! 𝒕 ← 𝒕 - 𝛼𝒖,
    ! 𝛽 ← ‖𝒕‖, 𝒖 ← 𝒕/𝛽,
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    !   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔, 
    ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
    ! 𝒕 ← 𝒕 - 𝛽𝒗,
    ! 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
    ! ----------------------
    if (present(precond)) then
      call precond%Apply(mesh, s, v, MatVec)
      call MatVec(mesh, t, s)
    else
      call MatVec(mesh, t, v)
    end if
    call Sub(mesh, t, t, u, alpha)
    beta = Norm_2(mesh, t); call Scale(mesh, u, t, 1.0_dp/beta)
    if (present(precond)) then
      call ConjMatVec(mesh, s, u)
      call conjPrecond%Apply(mesh, t, s, ConjMatVec)
    else
      call ConjMatVec(mesh, t, u)
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
    rho = hypot(rhoBar, beta)
    cs = rhoBar/rho; sn = beta/rho
    theta = sn*alpha; rhoBar = -cs*alpha
    phi = cs*phiBar; phiBar = sn*phiBar

    ! ----------------------
    ! Update 𝒛-solution:
    ! 𝒛 ← 𝒛 + (𝜑/𝜌)𝒘,
    ! 𝒘 ← 𝒗 - (𝜃/𝜌)𝒘.
    ! Check convergence for 𝜑̅ and 𝜑̅/𝜑̃.
    ! ( 𝜑̅ and 𝜑̃ implicitly contain residual norms;
    !   𝜑̅|𝜌̅| implicitly contain (𝓐[𝓟])*-residual norms, ‖(𝓐[𝓟])*𝒓‖. )
    ! ----------------------
    call Add(mesh, z, z, w, phi/rho)
    call Sub(mesh, w, v, w, theta/rho)
    if (params%Check(phiBar, phiBar/phiTilde)) exit
  end do

  ! ----------------------
  ! Compute 𝒙-solution:
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  !   𝒕 ← 𝓟𝒛, 𝒙 ← 𝒙 + 𝒕.
  ! 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒙 + 𝒛. 𝗲𝗻𝗱 𝗶𝗳
  ! ----------------------
  if (present(precond)) then
    call precond%Apply(mesh, t, z, MatVec)
    call Add(mesh, x, x, t)
  else
    call Add(mesh, x, x, z)
  end if

end subroutine Solve_LSQR

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a right preconditioned linear self-adjoint: least squares
!! problem: ‖𝓐[𝓟]𝒚 - 𝒃‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚, using the LSQR method.
!!
!! LSQR is not recommended in the self-adjoint case,
!! please consider MINRES instead.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_SymmLSQR(mesh, x, b, MatVec, params, precond)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  class(tConvParams), intent(inout) :: params
  class(tPreconditioner), intent(inout), optional :: precond
  procedure(tMatVecFunc) :: MatVec

  if (present(precond)) then
    call Solve_LSQR(mesh, x, b, MatVec, MatVec, params, precond, precond)
  else
    call Solve_LSQR(mesh, x, b, MatVec, MatVec, params)
  end if
  
end subroutine Solve_SymmLSQR

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a right preconditioned linear least squares problem:
!! ‖𝓐[𝓟]𝒚 - 𝒃‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚, using the LSMR method.
!!
!! LSMR is algebraically equivalent to applying MINRES 
!! to the normal equations: (𝓐[𝓟])*𝓐[𝓟]𝒚 = (𝓐[𝓟])*𝒃, 𝒙 = [𝓟]𝒚, 
!! (or, equivalently, [𝓟*]𝓐*𝓐[𝓟]𝒚 = [𝓟*]𝓐*𝒃, 𝒙 = [𝓟]𝒚),
!! but but has better numerical properties.
!! 
!! The normal equation's residual norm ‖(𝓐[𝓟])*(𝓐[𝓟]𝒚 - 𝒃)‖ 
!! decreases monotonically, while the residual norm ‖𝓐[𝓟]𝒚 - 𝒃‖₂   
!! is not guaranteed to decrease (but decreases on practice).
!!
!! References:
!! [1] Fong, D. C. and M. Saunders. 
!!     “LSMR: An Iterative Algorithm for Sparse Least-Squares Problems.” 
!!     SIAM J. Sci. Comput. 33 (2011): 2950-2971.
!! [2] Karimi, S., D. K. Salkuyeh and F. Toutounian. 
!!     “A preconditioner for the LSQR algorithm.” 
!!     Journal of applied mathematics & informatics 26 (2008): 213-222.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_LSMR(mesh, x, b, MatVec, &
    & ConjMatVec, params, precond, conjPrecond)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  class(tConvParams), intent(inout) :: params
  class(tPreconditioner), intent(inout), optional :: precond, conjPrecond
  procedure(tMatVecFunc) :: MatVec, ConjMatVec

  real(dp) :: alpha, alphaBar, beta, rho, rhoBar, cs, sn, &
    & theta, thetaBar, psi, psiBar, psiTilde, zeta, csBar, snBar
  type(tArray) :: r, s, t, w, h, u, v, z
  
  call AllocArray(t, r, u, v, w, h, z, mold=x)
  if (present(precond)) then
    call AllocArray(s, mold=x)
    call precond%Init(mesh, MatVec)
    call conjPrecond%Init(mesh, MatVec)
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
  !   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔, 
  ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
  ! 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
  ! ----------------------
  call MatVec(mesh, r, x)
  call Sub(mesh, r, b, r)
  beta = Norm_2(mesh, r); call Scale(mesh, u, r, 1.0_dp/beta)
  if (present(precond)) then
    call ConjMatVec(mesh, s, u)
    call conjPrecond%Apply(mesh, t, s, ConjMatVec)
  else
    call ConjMatVec(mesh, t, u)
  end if
  alpha = Norm_2(mesh, t); call Scale(mesh, v, t, 1.0_dp/alpha)

  ! ----------------------
  ! 𝛼̅ ← 𝛼, 𝜓̅ ← 𝛼𝛽,
  ! 𝜁 ← 1, 𝑐̅𝑠̅ ← 1, 𝑠̅𝑛̅ ← 0,
  ! 𝒛 ← {0}ᵀ,
  ! 𝒘 ← 𝒗, 𝒉 ← {0}ᵀ.
  ! ----------------------
  alphaBar = alpha; psiBar = alpha*beta
  zeta = 1.0_dp; csBar = 1.0_dp; snBar = 0.0_dp
  call Fill(mesh, z, 0.0_dp)
  call Set(mesh, w, v); call Fill(mesh, h, 0.0_dp)

  ! ----------------------
  ! 𝜓̃ ← 𝜓̅,
  ! Check convergence for 𝜓̃.
  ! ----------------------
  psiTilde = psiBar
  if (params%Check(psiTilde)) return
  
  do
    ! ----------------------
    ! Continue the bidiagonalization:
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲: 
    !   𝒔 ← 𝓟𝒗, 𝒕 ← 𝓐𝒔,
    ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐𝒗, 𝗲𝗻𝗱 𝗶𝗳
    ! 𝒕 ← 𝒕 - 𝛼𝒖,
    ! 𝛽 ← ‖𝒕‖, 𝒖 ← 𝒕/𝛽,
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    !   𝒔 ← 𝓐*𝒖, 𝒕 ← 𝓟*𝒔, 
    ! 𝗲𝗹𝘀𝗲: 𝒕 ← 𝓐*𝒖, 𝗲𝗻𝗱 𝗶𝗳
    ! 𝒕 ← 𝒕 - 𝛽𝒗,
    ! 𝛼 ← ‖𝒕‖, 𝒗 ← 𝒕/𝛼.
    ! ----------------------
    if (present(precond)) then
      call precond%Apply(mesh, s, v, MatVec)
      call MatVec(mesh, t, s)
    else
      call MatVec(mesh, t, v)
    end if
    call Sub(mesh, t, t, u, alpha)
    beta = Norm_2(mesh, t); call Scale(mesh, u, t, 1.0_dp/beta)
    if (present(precond)) then
      call ConjMatVec(mesh, s, u)
      call conjPrecond%Apply(mesh, t, s, ConjMatVec)
    else
      call ConjMatVec(mesh, t, u)
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
    rho = hypot(alphaBar, beta)
    cs = alphaBar/rho; sn = beta/rho
    theta = sn*alpha; alphaBar = cs*alpha
    thetaBar = snBar*rho; rhoBar = hypot(csBar*rho, theta)
    csBar = csBar*rho/rhoBar; snBar = theta/rhoBar
    psi = csBar*psiBar; psiBar = -snBar*psiBar

    ! ----------------------
    ! Update 𝒛-solution:
    ! 𝒉 ← 𝒘 - (𝜃𝜌/𝜁)𝒉, 𝜁 ← 𝜌𝜌̅,
    ! 𝒛 ← 𝒛 + (𝜓/𝜁)𝒉,
    ! 𝒘 ← 𝒗 - (𝜃/𝜌)𝒘.
    ! Check convergence for |𝜓̅| and |𝜓̅/𝜓̃|.
    ! ( |𝜓̅| and |𝜓̃| implicitly contain (𝓐[𝓟])*-residual norms, ‖(𝓐[𝓟])*𝒓‖. )
    ! ----------------------
    call Sub(mesh, h, w, h, thetaBar*rho/zeta); zeta = rho*rhoBar
    call Add(mesh, z, z, h, psi/zeta)
    call Sub(mesh, w, v, w, theta/rho)
    if (params%Check(abs(psiBar), abs(psiBar/psiTilde))) exit
  end do

  ! ----------------------
  ! Compute 𝒙-solution:
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  !   𝒕 ← 𝓟𝒛, 𝒙 ← 𝒙 + 𝒕.
  ! 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒙 + 𝒛. 𝗲𝗻𝗱 𝗶𝗳
  ! ----------------------
  if (present(precond)) then
    call precond%Apply(mesh, t, z, MatVec)
    call Add(mesh, x, x, t)
  else
    call Add(mesh, x, x, z)
  end if

end subroutine Solve_LSMR

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a right preconditioned linear self-adjoint: least squares
!! problem: ‖𝓐[𝓟]𝒚 - 𝒃‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚, using the LSMR method.
!!
!! Using LSMR is not recommended in the self-adjoint case,
!! please consider MINRES instead.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_SymmLSMR(mesh, x, b, MatVec, params, precond)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  class(tConvParams), intent(inout) :: params
  class(tPreconditioner), intent(inout), optional :: precond
  procedure(tMatVecFunc) :: MatVec

  if (present(precond)) then
    call Solve_LSMR(mesh, x, b, MatVec, MatVec, params, precond, precond)
  else
    call Solve_LSMR(mesh, x, b, MatVec, MatVec, params)
  end if
  
end subroutine Solve_SymmLSMR

end module StormRuler_Solvers_LSQR
