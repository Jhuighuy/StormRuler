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
module StormRuler_Solvers_Krylov

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp
use StormRuler_Helpers, only: SafeDivide
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: Fill, Set, Dot, Add, Sub
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_Base, only: &
  & @{tMatVecFunc$$, tPreconditionerFunc$$@|@0, NUM_RANKS}@

#$if HAS_OpenMP
use :: omp_lib
#$endif

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_CG
#$do rank = 0, NUM_RANKS
  module procedure Solve_CG$rank
  module procedure Solve_PCG$rank
#$end do
end interface Solve_CG

interface Solve_BiCGStab
#$do rank = 0, NUM_RANKS
  module procedure Solve_BiCGStab$rank
#$end do
end interface Solve_BiCGStab

#$do rank = 0, NUM_RANKS
type :: tPrecondEnv_Jacobi$rank
  real(dp), allocatable :: diag_inv(@:,:)
end type !tPrecondEnv_Jacobi$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!!
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Precondition_Jacobi$rank(mesh, Pu, u, MatVec, env, precond_env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in), target :: u(@:,:)
  real(dp), intent(inout), target :: Pu(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  class(*), intent(inout), allocatable, target :: precond_env
  ! >>>>>>>>>>>>>>>>>>>>>>

  class(tPrecondEnv_Jacobi$rank), pointer :: jacobi_env

  ! ----------------------
  ! Cast Jacobi preconditioner environment.
  ! ----------------------
  if (.not.allocated(precond_env)) then
    allocate(tPrecondEnv_Jacobi$rank :: precond_env)
  end if
  select type(precond_env)
    class is(tPrecondEnv_Jacobi$rank)
      jacobi_env => precond_env
  end select

  ! ----------------------
  ! Build the Jacobi preconditioner.
  ! ----------------------
  if (.not.allocated(jacobi_env%diag_inv)) then; block
    ! ----------------------
    ! TODO: we definitely need some more fashionable API to do this:
    ! ----------------------
    integer :: iCell
    
    real(dp), allocatable :: e(@:,:)

    allocate(e, jacobi_env%diag_inv, mold=u)
    call Fill(mesh, e, 0.0_dp)
    call Fill(mesh, jacobi_env%diag_inv, 0.0_dp)

    !$omp parallel do firstprivate(e)
    do iCell = 1, mesh%NumCells 
      
      e(@:,iCell) = 1.0_dp
      call mesh%SetRange(iCell)
      call MatVec(mesh, jacobi_env%diag_inv, e, env)
      e(@:,iCell) = 0.0_dp
      jacobi_env%diag_inv(@:,iCell) = 1.0_dp/jacobi_env%diag_inv(@:,iCell)
    end do
    !$omp end parallel do
    call mesh%SetRange()

  end block; end if

  ! ----------------------
  ! Apply the Jacobi preconditioner.
  ! ----------------------
  Pu(@:,:) = jacobi_env%diag_inv(@:,:)*u(@:,:)

end subroutine Precondition_Jacobi$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint definite 
!! operator equation: Au = b, using the Conjugate Gradients method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_CG$rank(mesh, u, b, MatVec, env, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  type(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  real(dp) :: alpha, beta, gamma, delta
  real(dp), allocatable :: p(@:,:), r(@:,:), t(@:,:)
  allocate(p, r, t, mold=u)

  ! ----------------------
  ! t ← Au,
  ! r ← b - t.
  ! ----------------------
  call MatVec(mesh, t, u, env)
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
    ! u ← u + α⋅p,
    ! r ← r - α⋅t,
    ! ----------------------
    call MatVec(mesh, t, p, env)
    alpha = SafeDivide(gamma, Dot(mesh, p, t))
    call Add(mesh, u, u, p, alpha)
    call Sub(mesh, r, r, t, alpha)

    ! ----------------------
    ! α ← <r⋅r>,
    ! check convergence for √α and √α/√δ.
    ! ----------------------
    alpha = Dot(mesh, r, r)
    if (params%Check(sqrt(alpha), sqrt(alpha/delta))) return

    ! ----------------------
    ! β ← α/γ,
    ! p ← r + β⋅p.
    ! γ ← α.
    ! ----------------------
    beta = SafeDivide(alpha, gamma)
    call Add(mesh, p, r, p, beta)
    gamma = alpha
  end do
end subroutine Solve_CG$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint definite operator equation: PAu = Pb, 
!! using the Preconditioned Conjugate Gradients method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_PCG$rank(mesh, u, b, MatVec, Precond, env, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  procedure(tPreconditionerFunc$rank) :: Precond
  class(*), intent(inout) :: env
  type(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  real(dp) :: alpha, beta, gamma, delta, theta
  real(dp), allocatable :: p(@:,:), r(@:,:), z(@:,:), t(@:,:)
  class(*), allocatable :: precond_env
  allocate(p, r, z, t, mold=u)

  ! ----------------------
  ! t ← Au,
  ! r ← b - t.
  ! ----------------------
  call MatVec(mesh, t, u, env)
  call Sub(mesh, r, b, t)

  ! ----------------------
  ! δ ← <r⋅r>,
  ! check convergence for √δ.
  ! ----------------------
  delta = Dot(mesh, r, r)
  if (params%Check(sqrt(delta))) return
  
  ! ----------------------
  ! z ← Pr,
  ! p ← z,
  ! γ ← <r⋅z>,
  ! ----------------------
  call Precond(mesh, z, r, MatVec, env, precond_env)
  call Set(mesh, p, z)
  gamma = Dot(mesh, r, z)

  do
    ! ----------------------
    ! t ← Ap,
    ! α ← γ/<p⋅t>,
    ! u ← u + α⋅p,
    ! r ← r - α⋅t,
    ! ----------------------
    call MatVec(mesh, t, p, env)
    alpha = SafeDivide(gamma, Dot(mesh, p, t))
    call Add(mesh, u, u, p, alpha)
    call Sub(mesh, r, r, t, alpha)

    ! ----------------------
    ! α ← <r⋅r>,
    ! check convergence for √α and √α/√δ.
    ! ----------------------
    alpha = Dot(mesh, r, r)
    if (params%Check(sqrt(alpha), sqrt(alpha/delta))) return

    ! ----------------------
    ! z ← Pr
    ! α ← <r⋅z>,
    ! ----------------------
    call Precond(mesh, z, r, MatVec, env, precond_env)
    alpha = Dot(mesh, r, z)

    ! ----------------------
    ! β ← α/γ,
    ! p ← z + β⋅p.
    ! γ ← α.
    ! ----------------------
    beta = SafeDivide(alpha, gamma)
    call Add(mesh, p, z, p, beta)
    gamma = alpha
  end do
end subroutine Solve_PCG$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: Au = b, using 
!! the good old Biconjugate Gradients (stabilized) method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_BiCGStab$rank(mesh, u, b, MatVec, env, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  type(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  real(dp) :: alpha, beta, gamma, delta, mu, rho, omega
  real(dp), allocatable :: h(@:,:), p(@:,:), r(@:,:), s(@:,:), t(@:,:), v(@:,:)
  allocate(h, p, r, s, t, v, mold=u)

  ! ----------------------
  ! t ← Au,
  ! r ← b - t.
  ! ----------------------
  call MatVec(mesh, t, u, env)
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
    call MatVec(mesh, v, p, env)
    
    ! ----------------------
    ! α ← ρ/<h⋅v>,
    ! s ← r - α⋅v,
    ! t ← As.
    ! ----------------------
    alpha = SafeDivide(rho, Dot(mesh, h, v))
    call Sub(mesh, s, r, v, alpha)
    call MatVec(mesh, t, s, env)
    
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

end module StormRuler_Solvers_Krylov
