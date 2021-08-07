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
submodule (StormRuler_MatrixFreeSolvers) StormRuler_MatrixFreeSolvers_Builtin

#$use 'StormRuler_Parameters.f90'

use StormRuler_Helpers, only: SafeDivide
use StormRuler_BLAS, only: Fill, Set, Dot, Add, Sub

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint definite 
!! operator equation: Au = b, using the Conjugate Gradients method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
module subroutine Solve_CG$rank(mesh, u, b, MatMul, env, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatMulFunc$rank) :: MatMul
  class(*), intent(in) :: env
  type(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  real(dp) :: alpha, beta, gamma, delta
  real(dp), allocatable :: p(@:,:), r(@:,:), t(@:,:)
  allocate(p, r, t, mold=u)

  ! ----------------------
  ! t ← Au,
  ! r ← b - t.
  ! ----------------------
  call MatMul(mesh, t, u, env)
  call Sub(mesh, r, b, t)
  ! ----------------------
  ! δ ← <r⋅r>,
  ! check convergence for √δ.
  ! ----------------------
  delta = Dot(mesh, r, r)
  if (params%Check(sqrt(delta))) return
  ! ----------------------
  ! p ← r.
  ! γ ← δ.
  ! ----------------------
  call Set(mesh, p, r)
  gamma = delta

  do
    ! ----------------------
    ! t ← Ap,
    ! α ← γ/<p⋅t>,
    ! u ← u + α⋅z,
    ! r ← r - α⋅g,
    ! ----------------------
    call MatMul(mesh, t, p, env)
    alpha = SafeDivide(gamma, Dot(mesh, p, t))
    call Add(mesh, u, u, p, alpha)
    call Sub(mesh, r, r, t, alpha)

    ! ----------------------
    ! α ← <r, r>,
    ! check convergence for √α and √α/√δ.
    ! ----------------------
    alpha = Dot(mesh, r, r)
    if (params%Check(sqrt(alpha), sqrt(alpha/delta))) return

    ! ----------------------
    ! β ← α/γ,
    ! p ← β⋅z + r.
    ! γ ← α.
    ! ----------------------
    beta = SafeDivide(alpha, gamma)
    call Add(mesh, p, r, p, beta)
    gamma = alpha
  end do
end subroutine Solve_CG$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: Au = b, using 
!! the good old Biconjugate Gradients (stabilized) method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
module subroutine Solve_BiCGStab$rank(mesh, u, b, MatMul, env, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatMulFunc$rank) :: MatMul
  class(*), intent(in) :: env
  type(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  real(dp) :: alpha, beta, gamma, delta, mu, rho, omega
  real(dp), allocatable :: &
    & h(@:,:), p(@:,:), r(@:,:), s(@:,:), t(@:,:), v(@:,:)
  allocate(h, p, r, s, t, v, mold=u)

  ! ----------------------
  ! t ← Au,
  ! r ← b - t.
  ! ----------------------
  call MatMul(mesh, t, u, env)
  call Sub(mesh, r, b, t)
  ! ----------------------
  ! δ ← <r⋅r>,
  ! check convergence for √δ.
  ! ----------------------
  delta = Dot(mesh, r, r)
  if (params%Check(sqrt(delta))) return
  ! ----------------------
  ! h ← r,
  ! p ← 0, v ← 0,
  ! ρ ← 1, α ← 1, ω ← 1. 
  ! ----------------------
  call Set(mesh, h, r)
  call Fill(mesh, p, 0.0_dp)
  call Fill(mesh, v, 0.0_dp)
  rho = 1.0_dp; alpha = 1.0_dp; omega = 1.0_dp

  do
    ! ----------------------
    ! μ ← <h⋅r>
    ! β ← (μ/ρ)⋅(α/ω),
    ! ρ ← μ.
    ! ----------------------
    mu = Dot(mesh, h, r)
    beta = SafeDivide(mu, rho)*SafeDivide(alpha, omega)
    rho = mu
    
    ! ----------------------
    ! p ← p - ω⋅v,
    ! p ← r + β⋅p,
    ! v ← Ap.
    ! ----------------------
    call Sub(mesh, p, p, v, omega)
    call Add(mesh, p, r, p, beta)
    call MatMul(mesh, v, p, env)
    
    ! ----------------------
    ! α ← ρ/<h⋅v>,
    ! s ← r - α⋅v,
    ! t ← As.
    ! ----------------------
    alpha = SafeDivide(rho, Dot(mesh, h, v))
    call Sub(mesh, s, r, v, alpha)
    call MatMul(mesh, t, s, env)
    
    ! ----------------------
    ! ω ← <t⋅s>/<t⋅t>,
    ! r ← s - ω⋅t,
    ! u ← u - ω⋅s,
    ! u ← u + α⋅p,
    ! ----------------------
    omega = SafeDivide(Dot(mesh, t, s), Dot(mesh, t, t))
    call Sub(mesh, r, s, t, omega)
    call Sub(mesh, u, u, s, omega)
    call Add(mesh, u, u, p, alpha)
    
    ! ----------------------
    ! γ ← <r⋅r>,
    ! check convergence for √γ and √γ/√δ.
    ! ----------------------
    gamma = Dot(mesh, r, r)
    if (params%Check(sqrt(gamma), sqrt(gamma/delta))) return
  end do
end subroutine Solve_BiCGStab$rank
#$end do

end submodule StormRuler_MatrixFreeSolvers_Builtin
