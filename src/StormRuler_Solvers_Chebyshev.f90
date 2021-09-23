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
use StormRuler_BLAS, only: Norm_2, Set, Fill, Add, Sub
#$for type_, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$type_
use StormRuler_Solvers_Precond, only: tPrecondFunc$type_
#$end for
use StormRuler_ConvParams, only: tConvParams
!use StormRuler_SolversEVP_Lanczos, only: EigenPairs_Lanczos

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint definite operator equation: 
!! [𝓟]𝓐𝒙 = [𝓟]𝒃, using the Chebyshev semi-iterative method.
!! Some accurate estimates of spectrum of [𝓟]𝓐 are required. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_Chebyshev(mesh, x, b, &
    & lambda_min, lambda_max, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: lambda_min, lambda_max, b(:,:)
  real(dp), intent(inout) :: x(:,:)
  procedure(tMatVecFuncR) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>

  logical :: first, second
  real(dp) :: c, d, alpha, beta, delta
  real(dp), pointer :: p(:,:), r(:,:), z(:,:)
  class(*), allocatable :: precond_env
  
  allocate(p, r, mold=x)
  if (present(Precond)) then
    allocate(z, mold=x)
  else
    z => r
  end if

  !call EigenPairs_Lanczos(mesh, x, b, MatVec, env, params, Precond)

  ! ----------------------
  ! 𝑐 ← ½(𝜆ₘₐₓ - 𝜆ₘᵢₙ),
  ! 𝑑 ← ½(𝜆ₘₐₓ + 𝜆ₘᵢₙ).
  ! ----------------------
  first = .true.; second = .true.
  c = 0.5_dp*(lambda_max - lambda_min)
  d = 0.5_dp*(lambda_max + lambda_min)

  ! ----------------------
  ! 𝒓 ← 𝓐𝒙,
  ! 𝒓 ← 𝒃 - 𝒓.
  ! ----------------------
  call MatVec(mesh, r, x, env)
  call Sub(mesh, r, b, r)

  ! ----------------------
  ! 𝛿 ← ‖𝒓‖₂,
  ! Check convergence for 𝛿.
  ! ----------------------
  delta = Norm_2(mesh, r)
  if (params%Check(delta)) return

  do
    ! ----------------------
    ! 𝒛 ← 𝓟𝒓,
    ! 𝗶𝗳 𝑘 == 1:
    !   𝛼 ← 1/𝑑,
    !   𝒑 ← 𝒛,
    ! 𝗲𝗹𝘀𝗲:
    !   𝗶𝗳 𝑘 == 2: 𝛽 ← ½(𝑐⋅𝛼)²,
    !   𝗲𝗹𝘀𝗲: 𝛽 ← (½⋅𝑐⋅𝛼)², 𝗲𝗻𝗱 𝗶𝗳
    !   𝛼 ← 1/(𝑑 - 𝛽/𝛼),
    !   𝒑 ← 𝒛 + 𝛽𝒑.
    ! 𝗲𝗻𝗱 𝗶𝗳
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
    ! 𝒙 ← 𝒙 + 𝛼𝒑,
    ! 𝒓 ← 𝓐𝒙,
    ! 𝒓 ← 𝒃 - 𝒓.
    ! ----------------------
    call Add(mesh, x, x, p, alpha)
    call MatVec(mesh, r, x, env)
    call Sub(mesh, r, b, r)

    ! ----------------------
    ! 𝛽 ← ‖𝒓‖₂,
    ! Check convergence for 𝛽 and 𝛽/𝛿.
    ! ----------------------
    beta = Norm_2(mesh, r)
    if (params%Check(beta, beta/delta)) exit
  end do
end subroutine Solve_Chebyshev

end module StormRuler_Solvers_Chebyshev
