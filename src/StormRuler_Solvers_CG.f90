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
module StormRuler_Solvers_CG

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp
use StormRuler_Helpers, only: SafeDivide

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: Fill, Set, Dot, Add, Sub
use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_ConvParams, only: tConvParams
use StormRuler_Precond, only: tPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_CG
  module procedure Solve_CG
end interface Solve_CG

interface Solve_BiCGStab
  module procedure Solve_BiCGStab
end interface Solve_BiCGStab

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve a linear self-adjoint definite operator equation: 
!! [𝓜]𝓐[𝓜ᵀ]𝒚 = [𝓜]𝒃, 𝒙 = [𝓜ᵀ]𝒚, [𝓜𝓜ᵀ = 𝓟], using the 
!! Conjugate Gradients (CG) method.
!!
!! CG may be applied to the consistent singular problems, 
!! it converges towards..
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_CG(mesh, x, b, MatVec, params, pre)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  class(tConvParams), intent(inout) :: params
  class(tPreconditioner), intent(inout), optional :: pre
  procedure(tMatVecFunc) :: MatVec
  
  real(dp) :: alpha, beta, gamma, delta
  type(tArray) :: p, r, t, z
  
  call AllocArray(p, r, t, mold=x)
  if (present(pre)) then
    call AllocArray(z, mold=x)
    call pre%Init(mesh, MatVec)
  else
    z = r
  end if

  ! ----------------------
  ! 𝒓 ← 𝓐𝒙,
  ! 𝒓 ← 𝒃 - 𝒕.
  ! ----------------------
  call MatVec(mesh, r, x)
  call Sub(mesh, r, b, r)

  ! ----------------------
  ! 𝛿 ← <𝒓⋅𝒓>,
  ! Check convergence for √𝛿.
  ! ----------------------
  delta = Dot(mesh, r, r)
  if (params%Check(sqrt(delta))) return
  
  ! ----------------------
  ! 𝒛 ← 𝓟𝒓,
  ! 𝒑 ← 𝒛,
  ! 𝛾 ← <𝒓⋅𝒛>,
  ! ----------------------
  if (present(pre)) then
    call pre%Apply(mesh, z, r, MatVec)
  end if
  call Set(mesh, p, z)
  gamma = Dot(mesh, r, z)

  do
    ! ----------------------
    ! 𝒕 ← 𝓐𝒑,
    ! 𝛼 ← 𝛾/<𝒑⋅𝒕>,
    ! 𝒙 ← 𝒙 + 𝛼𝒑,
    ! 𝒓 ← 𝒓 - 𝛼𝒕,
    ! ----------------------
    call MatVec(mesh, t, p)
    alpha = SafeDivide(gamma, Dot(mesh, p, t))
    call Add(mesh, x, x, p, alpha)
    call Sub(mesh, r, r, t, alpha)

    ! ----------------------
    ! 𝛼 ← <𝒓⋅𝒓>,
    ! Check convergence for √𝛼 and √𝛼/√𝛿.
    ! ----------------------
    alpha = Dot(mesh, r, r)
    if (params%Check(sqrt(alpha), sqrt(alpha/delta))) exit

    ! ----------------------
    ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
    !   𝒛 ← 𝓟𝒓,
    !   𝛼 ← <𝒓⋅𝒛>,
    ! 𝗲𝗻𝗱 𝗶𝗳 // otherwise 𝒛 ≡ 𝒓, 𝛼 unchanged.  
    ! ----------------------
    if (present(pre)) then
      call pre%Apply(mesh, z, r, MatVec)
      alpha = Dot(mesh, r, z)
    end if

    ! ----------------------
    ! 𝛽 ← 𝛼/𝛾,
    ! 𝒑 ← 𝒛 + 𝛽𝒑,
    ! 𝛾 ← 𝛼.
    ! ----------------------
    beta = SafeDivide(alpha, gamma)
    call Add(mesh, p, z, p, beta)
    gamma = alpha

  end do

end subroutine Solve_CG

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, using 
!! the good old Biconjugate Gradients (stabilized) method (BiCGStab).
!!
!! BiCGStab may be applied to the consistent singular problems,
!! it converges towards..
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_BiCGStab(mesh, x, b, MatVec, params, pre)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  class(tConvParams), intent(inout) :: params
  class(tPreconditioner), intent(inout), optional :: pre
  procedure(tMatVecFunc) :: MatVec

  real(dp) :: alpha, beta, gamma, delta, mu, rho, omega
  type(tArray) :: p, r, rTilde, s, t, v, w, y, z

  call AllocArray(p, r, rTilde, s, t, v, mold=x)
  if (present(pre)) then
    call AllocArray(w, y, z, mold=x)
    call pre%Init(mesh, MatVec)
  else
    w = t; y = p; z = s
  end if

  ! ----------------------
  ! 𝒓 ← 𝓐𝒙,
  ! 𝒓 ← 𝒃 - 𝒓.
  ! ----------------------
  call MatVec(mesh, r, x)
  call Sub(mesh, r, b, r)

  ! ----------------------
  ! 𝛿 ← <𝒓⋅𝒓>,
  ! Check convergence for √𝛿.
  ! ----------------------
  delta = Dot(mesh, r, r)
  if (params%Check(sqrt(delta))) return
  
  ! ----------------------
  ! 𝒓̃ ← 𝒓,
  ! 𝒑 ← {0}ᵀ, 𝒗 ← {0}ᵀ,
  ! 𝜌 ← 1, 𝛼 ← 1, 𝜔 ← 1. 
  ! ----------------------
  call Set(mesh, rTilde, r)
  call Fill(mesh, p, 0.0_dp)
  call Fill(mesh, v, 0.0_dp)
  rho = 1.0_dp; alpha = 1.0_dp; omega = 1.0_dp

  do
    ! ----------------------
    ! 𝜇 ← <𝒓̃⋅𝒓>,
    ! 𝛽 ← (𝜇/𝜌)⋅(𝛼/𝜔),
    ! 𝜌 ← 𝜇.
    ! ----------------------
    mu = Dot(mesh, rTilde, r)
    beta = SafeDivide(mu, rho)*SafeDivide(alpha, omega)
    rho = mu
    
    ! ----------------------
    ! 𝒑 ← 𝒑 - 𝜔𝒗,
    ! 𝒑 ← 𝒓 + 𝛽𝒑,
    ! 𝒚 ← 𝓟𝒑,
    ! 𝒗 ← 𝓐𝒚.
    ! ----------------------
    call Sub(mesh, p, p, v, omega)
    call Add(mesh, p, r, p, beta)
    if (present(pre)) then
      call pre%Apply(mesh, y, p, MatVec)
    end if
    call MatVec(mesh, v, y)
    
    ! ----------------------
    ! 𝛼 ← 𝜌/<𝒓̃⋅𝒗>,
    ! 𝒔 ← 𝒓 - 𝛼𝒗,
    ! 𝒛 ← 𝓟𝒔,
    ! 𝒕 ← 𝓐𝒛.
    ! ----------------------
    alpha = SafeDivide(rho, Dot(mesh, rTilde, v))
    call Sub(mesh, s, r, v, alpha)
    if (present(pre)) then
      call pre%Apply(mesh, z, s, MatVec)
    end if
    call MatVec(mesh, t, z)
    
    ! ----------------------
    ! 𝒘 ← 𝓟𝒕,
    ! 𝜔 ← <𝒘⋅𝒛>/<𝒘⋅𝒘>,
    ! 𝒓 ← 𝒔 - 𝜔𝒕,
    ! 𝒙 ← 𝒙 + 𝜔𝒛,
    ! 𝒙 ← 𝒙 + 𝛼𝒚,
    ! ----------------------
    if (present(pre)) then
      call pre%Apply(mesh, w, t, MatVec)
    end if
    omega = SafeDivide(Dot(mesh, w, z), Dot(mesh, w, w))
    call Sub(mesh, r, s, t, omega)
    call Add(mesh, x, x, z, omega)
    call Add(mesh, x, x, y, alpha)
    
    ! ----------------------
    ! 𝛾 ← <𝒓⋅𝒓>,
    ! Check convergence for √𝛾 and √𝛾/√𝛿.
    ! ----------------------
    gamma = Dot(mesh, r, r)
    if (params%Check(sqrt(gamma), sqrt(gamma/delta))) exit

  end do

end subroutine Solve_BiCGStab

end module StormRuler_Solvers_CG
