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
module StormRuler_Solvers_Chebyshev

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@, &
  & Dot, Norm_2, Fill, Fill_Random, Set, Scale, Add, Sub
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_Base, only: @{tPrecondFunc$$@|@0, NUM_RANKS}@
use StormRuler_Solvers_Base, only: ResidualNorm

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_Chebyshev
#$do rank = 0, NUM_RANKS
  module procedure Solve_Chebyshev$rank
#$end do
end interface Solve_Chebyshev

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint definite operator equation: 
!! [P¹ᐟ²]A[P¹ᐟ²]y = [P¹ᐟ²]b, [P¹ᐟ²]y = u, using the 
!! Chebyshev semi-iterative method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_Chebyshev$rank(mesh, x, b, &
    & lambda_min, lambda_max, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: lambda_min, lambda_max, b(@:,:)
  real(dp), intent(inout) :: x(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFunc$rank), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>

  logical :: first, second
  real(dp) :: c, d, alpha, beta, gamma, delta, mu_min, mu_max
  real(dp), pointer :: p(@:,:), r(@:,:), t(@:,:), z(@:,:)
  class(*), allocatable :: precond_env
  
  allocate(p, r, t, mold=x)
  if (present(Precond)) then
    allocate(z, mold=x)
  else
    z => r
  end if

  first = .true.; second = .true.
  c = 0.5_dp*(lambda_max - lambda_min)
  d = 0.5_dp*(lambda_max + lambda_min)

  ! ----------------------
  ! ----------------------
  ! ----------------------

  !mu_max = 163.113893827740_dp

  !call Fill_Random(mesh, t)
  !do 
  !  alpha = Norm_2(mesh, t); call Scale(mesh, t, t, 1.0_dp/alpha)
  !  call MatVec(mesh, p, t, env)
  !  call Sub(mesh, p, p, t, mu_max)
  !  mu_min = Dot(mesh, p, t)
  !  print *, '***', mu_min
  !  call Set(mesh, t, p)
  !end do

  ! ----------------------
  ! ----------------------
  ! ----------------------

  ! ----------------------
  ! t ← Ax,
  ! r ← b - t.
  ! ----------------------
  call Fill(mesh, x, 0.0_dp)
  call MatVec(mesh, t, x, env)
  call Sub(mesh, r, b, t)

  ! ----------------------
  ! δ ← <r⋅r>,
  ! Check convergence for √δ.
  ! ----------------------
  delta = Dot(mesh, r, r)
  if (params%Check(sqrt(delta))) return

  do
    ! ----------------------
    ! z ← Pr.
    ! ----------------------
    if (present(Precond)) &
      & call Precond(mesh, z, r, MatVec, env, precond_env)

    if (first) then
      first = .false.
      ! ----------------------
      ! α ← 1/d,
      ! p ← z.
      ! ----------------------
      alpha = 1.0_dp/d
      call Set(mesh, p, z)
    else
      ! ----------------------
      ! β ← (cα)²/2, OR: β ← (cα/2)²
      ! α ← 1/(d - β/α),
      ! p ← z + βp.
      ! ----------------------
      if (second) then
        second = .false.
        beta = 0.5_dp*(c*alpha)**2
      else
        beta = (0.5_dp*c*alpha)**2
      end if
      alpha = 1.0_dp/(d - beta/alpha)
      call Add(mesh, p, z, p, beta)
    end if

    ! ----------------------
    ! x ← x + αp,
    ! t ← Ax,
    ! r ← b - t.
    ! ----------------------
    call Add(mesh, x, x, p, alpha)
    call MatVec(mesh, t, x, env)
    call Sub(mesh, r, b, t)

    ! ----------------------
    ! β ← <r⋅r>,
    ! Check convergence for √β and √β/√δ.
    ! ----------------------
    beta = Dot(mesh, r, r)
    if (params%Check(sqrt(beta), sqrt(beta/delta))) exit
  end do
end subroutine Solve_Chebyshev$rank
#$end do

end module StormRuler_Solvers_Chebyshev
