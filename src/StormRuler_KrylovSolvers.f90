!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module StormRuler_KrylovSolvers

use StormRuler_Helpers
use StormRuler_Arithmetics
use StormRuler_Mesh
#$use 'StormRuler_Parameters.f90'

implicit none

!! -----------------------------------------------------------------  
type :: ConvParameters
  integer :: Iteration
  integer :: NumIterations
  real(dp) :: AbsoluteTolerance
  real(dp) :: RelativeTolerance
contains
  procedure :: Init => ConvParameters_Init
  procedure :: Check => ConvParameters_Check
end type
private ConvParameters_Init, ConvParameters_Check

abstract interface
#$do rank = 0, NUM_RANKS
  subroutine MeshOperator$rank(mesh,u,c,opParams)
    import Mesh2D, dp
    class(Mesh2D), intent(in) :: mesh
    real(dp), intent(inout) :: u(@:,:),c(@:,:)
    class(*), intent(in) :: opParams
  end subroutine MeshOperator$rank
#$end do
end interface

interface Solve_CG
#$do rank = 0, NUM_RANKS
  module procedure Solve_CG$rank
#$end do
end interface Solve_CG

interface Solve_BiCGStab
#$do rank = 0, NUM_RANKS
  module procedure Solve_BiCGStab$rank
#$end do
end interface Solve_BiCGStab

contains

!! -----------------------------------------------------------------  
!! Initialize iteration parameters.
subroutine ConvParameters_Init(params &
    ,absoluteTolerance,relativeTolerance,numIterations)
  class(ConvParameters), intent(out) :: params
  real(dp), intent(in) :: absoluteTolerance
  real(dp), intent(in), optional :: relativeTolerance
  integer, intent(in), optional :: numIterations
  ! ----------------------
  ! Initialize errors.
  call EnsurePositive(absoluteTolerance)
  params%AbsoluteTolerance = absoluteTolerance
  if (present(relativeTolerance)) then
    call EnsurePositive(relativeTolerance)
    params%RelativeTolerance = relativeTolerance
  end if
  ! ----------------------
  ! Initialize iterations.
  params%Iteration = 0
  if (present(numIterations)) then
    call EnsurePositive(relativeTolerance)
    params%NumIterations = numIterations
  end if
end subroutine ConvParameters_Init

!! -----------------------------------------------------------------  
!! Check convergence of the iterations process.  
function ConvParameters_Check(params &
    ,absoluteError,relativeError) result(converged)
  class(ConvParameters), intent(inout) :: params
  real(dp), intent(in) :: absoluteError
  real(dp), intent(in), optional :: relativeError
  logical :: converged
  ! ----------------------
  ! Check whether number of iterations has exceeded.
  associate(iteration=>params%Iteration &
      , numIterations=>params%NumIterations)
    iteration = iteration + 1
    converged = converged.or.(iteration>=numIterations)
  end associate
  ! ----------------------
  ! Check convergence by absolute and relative errors.
  associate(absoluteTolerance=>params%AbsoluteTolerance &
          , relativeTolerance=>params%RelativeTolerance)
    call EnsureNonNegative(absoluteError)
    !print *, absoluteError
    converged = converged.or.(absoluteError<=absoluteTolerance)
    if (present(relativeError).and.(relativeTolerance > 0)) then
      call EnsureNonNegative(relativeError)
      converged = converged.or.(relativeError<=relativeTolerance)
    end if
  end associate
end function ConvParameters_Check

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
!! Solve a linear self-adjoint definite 
!! operator equation using the Conjugate Gradients method.
#$do rank = 0, NUM_RANKS
subroutine Solve_CG$rank(mesh &
                        ,u,b,LOp,opParams,convParams)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(inout) :: u(@:,:),b(@:,:)
  procedure(MeshOperator$rank) :: LOp
  class(*), intent(in) :: opParams
  type(ConvParameters), intent(inout) :: convParams
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp) :: alpha,beta,gamma,delta
  real(dp), allocatable :: p(@:,:),r(@:,:),t(@:,:)
  allocate(p,r,t, mold=u)
  ! ----------------------
  ! t ← Au,
  ! r ← b - t.
  call LOp(mesh,t,u,opParams)
  call Sub(mesh,r,b,t)
  ! δ ← <r⋅r>,
  ! check convergence for √δ.
  delta = Dot(mesh,r,r)
  if (convParams%Check(sqrt(delta))) return
  ! p ← r.
  ! γ ← δ.
  call Set(mesh,p,r)
  gamma = delta
  ! ----------------------
  do
    ! t ← Ap,
    ! α ← γ/<p⋅t>,
    ! u ← u + α⋅z,
    ! r ← r - α⋅g,
    call LOp(mesh,t,p,opParams)
    alpha = SafeDivide(gamma,Dot(mesh,p,t))
    call Add(mesh,u,u,p,alpha)
    call Sub(mesh,r,r,t,alpha)
    ! α ← <r,r>,
    ! check convergence for √α and √α/√δ.
    alpha = Dot(mesh,r,r)
    if (convParams%Check(sqrt(alpha),sqrt(alpha/delta))) return
    ! β ← α/γ,
    ! p ← β⋅z + r.
    beta = SafeDivide(alpha,gamma)
    call Add(mesh,p,r,p,beta)
    ! γ ← α.
    gamma = alpha
  end do
end subroutine Solve_CG$rank
#$end do

!! -----------------------------------------------------------------  
!! Solve a linear operator equation using 
!! the good old Biconjugate Gradients (stabilized) method.
#$do rank = 0, NUM_RANKS
subroutine Solve_BiCGStab$rank(mesh &
                              ,u,b,LOp,opParams,convParams)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(inout) :: u(@:,:),b(@:,:)
  procedure(MeshOperator$rank) :: LOp
  class(*), intent(in) :: opParams
  type(ConvParameters), intent(inout) :: convParams
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp) :: alpha,beta,gamma,delta,mu,rho,omega
  real(dp), allocatable :: h(@:,:),p(@:,:),r(@:,:),s(@:,:),t(@:,:),v(@:,:)
  allocate(h,p,r,s,t,v, mold=u)
  ! ----------------------
  ! t ← Au,
  ! r ← b - t.
  call LOp(mesh,t,u,opParams)
  call Sub(mesh,r,b,t)
  ! δ ← <r⋅r>,
  ! check convergence for √δ.
  delta = Dot(mesh,r,r)
  if (convParams%Check(sqrt(delta))) return
  ! h ← r,
  ! p ← 0, v ← 0,
  ! ρ ← 1, α ← 1, ω ← 1. 
  call Set(mesh,h,r)
  call Fill(mesh,p)
  call Fill(mesh,v)
  rho = 1.0_dp; alpha = 1.0_dp; omega = 1.0_dp
  ! ----------------------
  do
    ! μ ← <h⋅r>
    ! β ← (μ/ρ)⋅(α/ω),
    ! ρ ← μ.
    mu = Dot(mesh,h,r)
    beta = SafeDivide(mu,rho)*SafeDivide(alpha,omega)
    rho = mu
    ! p ← p - ω⋅v,
    ! p ← r + β⋅p,
    ! v ← Ap.
    call Sub(mesh,p,p,v,omega)
    call Add(mesh,p,r,p,beta)
    call LOp(mesh,v,p,opParams)
    ! α ← ρ/<h⋅v>,
    ! s ← r - α⋅v,
    ! t ← As.
    alpha = SafeDivide(rho,Dot(mesh,h,v))
    call Sub(mesh,s,r,v,alpha)
    call LOp(mesh,t,s,opParams)
    ! ω ← <t⋅s>/<t⋅t>,
    ! r ← s - ω⋅t,
    ! u ← u - ω⋅s,
    ! u ← u + α⋅p,
    omega = SafeDivide(Dot(mesh,t,s),Dot(mesh,t,t))
    call Sub(mesh,r,s,t,omega)
    call Sub(mesh,u,u,s,omega)
    call Add(mesh,u,u,p,alpha)
    ! γ ← <r,r>,
    ! check convergence for √γ and √γ/√δ.
    gamma = Dot(mesh,r,r)
    if (convParams%Check(sqrt(gamma),sqrt(gamma/delta))) return
  end do
end subroutine Solve_BiCGStab$rank
#$end do

end module StormRuler_KrylovSolvers
