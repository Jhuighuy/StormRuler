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

use StormRuler_Parameters, only: dp
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@, &
  & Norm_2, Set, Fill, Add, Sub
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
!! PAx = Pb, using the Chebyshev semi-iterative method.
!! Some accurate estimates of spectrum of PA are required. 
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
  real(dp) :: c, d, alpha, beta, delta
  real(dp), pointer :: p(@:,:), r(@:,:), z(@:,:)
  class(*), allocatable :: precond_env
  
  allocate(p, r, mold=x)
  if (present(Precond)) then
    allocate(z, mold=x)
  else
    z => r
  end if

  ! ----------------------
  ! c ← (λₘₐₓ - λₘᵢₙ)/2,
  ! d ← (λₘₐₓ + λₘᵢₙ)/2.
  ! ----------------------
  first = .true.; second = .true.
  c = 0.5_dp*(lambda_max - lambda_min)
  d = 0.5_dp*(lambda_max + lambda_min)

  ! ----------------------
  ! r ← Ax,
  ! r ← b - r.
  ! ----------------------
  call MatVec(mesh, r, x, env)
  call Sub(mesh, r, b, r)

  ! ----------------------
  ! δ ← ‖r‖₂,
  ! Check convergence for δ.
  ! ----------------------
  delta = Norm_2(mesh, r)
  if (params%Check(delta)) return

  do
    ! ----------------------
    ! z ← Pr,
    ! IF k == 1:
    !   α ← 1/d,
    !   p ← z,
    ! ELSE:
    !   IF k == 2: β ← (c⋅α)²/2,
    !   ELSE: β ← (c⋅α/2)², END IF
    !   α ← 1/(d - β/α),
    !   p ← z + βp.
    ! END IF
    ! ----------------------
    if (present(Precond)) &
      & call Precond(mesh, z, r, MatVec, env, precond_env)
    if (first) then
      first = .false.
      alpha = 1.0_dp/d
      call Set(mesh, p, z)
    else
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
    ! r ← Ax,
    ! r ← b - r.
    ! ----------------------
    call Add(mesh, x, x, p, alpha)
    call MatVec(mesh, r, x, env)
    call Sub(mesh, r, b, r)

    ! ----------------------
    ! β ← ‖r‖₂,
    ! Check convergence for β and β/δ.
    ! ----------------------
    beta = Norm_2(mesh, r)
    if (params%Check(beta, beta/delta)) exit
  end do
end subroutine Solve_Chebyshev$rank
#$end do

end module StormRuler_Solvers_Chebyshev
