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
module StormRuler_Solvers_Newton

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: Dot, Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_BLAS, only: tMatVecFunc, tBiMatVecFunc

use StormRuler_Solvers, only: LinSolve
use StormRuler_ConvParams, only: tConvParams

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a nonlinear operator equation: 𝓐(𝒙) = 𝒃,
!! where 𝓙(𝒙) ≈ ∂𝓐(𝒙)/∂𝒙, using the Newton's method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_Newton(mesh, MatVec, JacobianMatVec, xArr, bArr, params)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: bArr
  class(tArray), intent(inout) :: xArr
  class(tConvParams), intent(inout) :: params
  procedure(tMatVecFunc) :: MatVec
  procedure(tBiMatVecFunc) :: JacobianMatVec

  type(tArray) :: tArr, rArr
  type(tConvParams) :: jacConvParams 

  call AllocArray(tArr, rArr, mold=xArr)

  ! ----------------------
  ! Newton's method:
  ! 𝓐(𝒚) ≈ 𝓐(𝒙) + 𝓙(𝒙)(𝒚 - 𝒙) = 𝒃, 
  ! 𝓙(𝒙) ≈ ∂𝓐(𝒙)/∂𝒙, 𝒚 ← 𝒙.
  ! Alternative formulation:
  ! 𝓙(𝒙)𝒕 = 𝒓, 𝒕 = 𝒚 - 𝒙, 𝒓 = 𝒃 - 𝓐(𝒙),
  ! ----------------------

  do
    ! ----------------------
    ! Compute residual:
    ! 𝒓 ← 𝓐(𝒙),
    ! 𝒓 ← 𝒃 - 𝒓,
    ! Check convergence for ‖𝒓‖.
    ! ----------------------
    call MatVec(mesh, rArr, xArr)
    call Sub(mesh, rArr, bArr, rArr)
    if (params%Check(Norm_2(mesh, rArr))) exit

    ! ----------------------
    ! Solve the Jacobian equation:
    ! 𝒕 ← 𝒓,
    ! 𝒕 ← 𝓙(𝒙)⁻¹𝒓,
    ! 𝒙 ← 𝒙 + 𝒕.
    ! ----------------------
    call Set(mesh, tArr, rArr)
    !! TODO: equation parameters!
    call jacConvParams%Init(1e-8_dp, 1e-8_dp, 2000, 'Newton')
    call LinSolve(mesh, 'BiCGStab', '', tArr, rArr, JacobianMatVecWithX, jacConvParams)
    call Add(mesh, xArr, xArr, tArr)

  end do

contains
  subroutine JacobianMatVecWithX(mesh, zArr, yArr)
    class(tMesh), intent(in), target :: mesh
    class(tArray), intent(inout), target :: yArr, zArr

    ! ----------------------
    ! 𝒛 ← 𝓙(𝒙)𝒚.
    ! ----------------------
    call JacobianMatVec(mesh, zArr, yArr, xArr)

  end subroutine JacobianMatVecWithX
end subroutine Solve_Newton

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a nonlinear operator equation: 𝓐(𝒙) = 𝒃,
!! using the Jacobian free-Newton-Krylov method.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_JFNK(mesh, MatVec, xArr, bArr, params)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: bArr
  class(tArray), intent(inout) :: xArr
  class(tConvParams), intent(inout) :: params
  procedure(tMatVecFunc) :: MatVec

  type(tArray) :: tArr, yArr, zArr

  !! TODO: reimplement me in the optimized manner!

  call AllocArray(tArr, yArr, zArr, mold=xArr)

  call Solve_Newton(mesh, MatVec, ApproxJacMatVec_1, xArr, bArr, params)

contains
  subroutine ApproxJacMatVec_1(mesh, Jx, xArr, xTildeArr)
    class(tMesh), intent(in), target :: mesh
    class(tArray), intent(inout), target :: xArr, xTildeArr, Jx

    real(dp), parameter :: epsilon = 1e-6_dp

    !! TODO: selection of 𝜀 is missing! 

    ! ----------------------
    ! Consider the first-order jacobian approximation:
    ! 𝓐(𝒙₀ + 𝜀𝒙) = 𝓐(𝒙₀) + 𝜀(∂𝓐(𝒙₀)/∂𝒙)𝒙 + 𝓞(𝜀²),
    ! 𝓙(𝒙₀)𝒙 ≈ [𝓐(𝒙₀ + 𝜀𝒙) - 𝓐(𝒙₀)]/𝜀 = (∂𝓐(𝒙₀)/∂𝒙)𝒙 + 𝓞(𝜀).
    ! ----------------------

    ! ----------------------
    ! 𝒕 ← 𝒙₀ + 𝜀𝒙.
    ! 𝒚 ← 𝓐(𝒙₀), 𝒛 ← 𝓐(𝒕),
    ! 𝓙𝒙 ← (1/𝜀)𝒛 - (1/𝜀)𝒚.
    ! ----------------------
    call Add(mesh, tArr, xTildeArr, xArr, epsilon)
    call MatVec(mesh, yArr, xTildeArr); call MatVec(mesh, zArr, tArr)
    call Sub(mesh, Jx, zArr, yArr, 1.0_dp/epsilon, 1.0_dp/epsilon)

  end subroutine ApproxJacMatVec_1
!  subroutine ApproxJacMatVec_2(mesh, Jx, xArr, xTildeArr)
!    class(tMesh), intent(inout), target :: mesh
!    class(tArray), intent(inout), target :: xArr, xTildeArr, Jx
!
!    real(dp), parameter :: epsilon = 1e-6_dp
!
!    ! ----------------------
!    ! Consider the second-order jacobian approximation:
!    ! 𝓐(𝒙₀ + 𝜀𝒙) = 𝓐(𝒙₀) + 𝜀(∂𝓐(𝒙₀)/∂𝒙)𝒙 + ½𝜀²(∂²𝓐(𝒙₀)/∂𝒙²)𝒙² + 𝓞(𝜀³),
!    ! 𝓐(𝒙₀ - 𝜀𝒙) = 𝓐(𝒙₀) - 𝜀(∂𝓐(𝒙₀)/∂𝒙)𝒙 + ½𝜀²(∂²𝓐(𝒙₀)/∂𝒙²)𝒙² + 𝓞(𝜀³),
!    ! 𝓙(𝒙₀)𝒙 ≈ [𝓐(𝒙₀ + 𝜀𝒙) - 𝓐(𝒙₀ - 𝜀𝒙)]/(2𝜀) = (∂𝓐(𝒙₀)/∂𝒙)𝒙 + 𝓞(𝜀²).
!    ! ----------------------
!
!    error stop 'unimplemented'
!
!  end subroutine ApproxJacMatVec_2
end subroutine Solve_JFNK

end module StormRuler_Solvers_Newton
