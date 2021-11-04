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
use StormRuler_Array, only: tArrayR, AllocArray

use StormRuler_BLAS, only: Norm_2, Set, Fill, Add, Sub
#$for T, _ in [SCALAR_TYPES[0]]
use StormRuler_BLAS, only: tMatVecFunc$T
use StormRuler_Solvers_Precond, only: tPreMatVecFunc$T
#$end for

use StormRuler_ConvParams, only: tConvParams
!use StormRuler_SolversEVP_Lanczos, only: EigenPairs_Lanczos

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_Chebyshev
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Solve_Chebyshev$T
#$end for
end interface Solve_Chebyshev

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint definite operator equation: 
!! [𝓟]𝓐𝒙 = [𝓟]𝒃, using the Chebyshev semi-iterative method.
!! Some accurate estimates of spectrum of [𝓟]𝓐 are required. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Solve_Chebyshev$T(mesh, x, b, &
    & minLambda, maxLambda, MatVec, params, PreMatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: b
  class(tArray$T), intent(inout) :: x
  procedure(tMatVecFunc$T) :: MatVec
  class(tConvParams), intent(inout) :: params
  procedure(tPreMatVecFunc$T), optional :: PreMatVec
  real(dp), intent(in) :: minLambda, maxLambda

  logical :: first, second
  real(dp) :: c, d, alpha, beta, delta
  type(tArrayR) :: p, r, z
  class(*), allocatable :: preEnv
  
  call AllocArray(p, r, mold=x)
  if (present(PreMatVec)) then
    call AllocArray(z, mold=x)
  else
    z = r
  end if

  !call EigenPairs_Lanczos(mesh, x, b, MatVec, params, PreMatVec)

  ! ----------------------
  ! 𝑐 ← ½(𝜆ₘₐₓ - 𝜆ₘᵢₙ),
  ! 𝑑 ← ½(𝜆ₘₐₓ + 𝜆ₘᵢₙ).
  ! ----------------------
  first = .true.; second = .true.
  c = 0.5_dp*(maxLambda - minLambda)
  d = 0.5_dp*(maxLambda + minLambda)

  ! ----------------------
  ! 𝒓 ← 𝓐𝒙,
  ! 𝒓 ← 𝒃 - 𝒓.
  ! ----------------------
  call MatVec(mesh, r, x)
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
    if (present(PreMatVec)) &
      & call PreMatVec(mesh, z, r, MatVec, preEnv)
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
    call MatVec(mesh, r, x)
    call Sub(mesh, r, b, r)

    ! ----------------------
    ! 𝛽 ← ‖𝒓‖₂,
    ! Check convergence for 𝛽 and 𝛽/𝛿.
    ! ----------------------
    beta = Norm_2(mesh, r)
    if (params%Check(beta, beta/delta)) exit
  end do
  
end subroutine Solve_Chebyshev$T
#$end for

end module StormRuler_Solvers_Chebyshev
