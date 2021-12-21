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

use StormRuler_Consts, only: dp
use StormRuler_Helpers, only: SafeInverse, SafeDivide

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: Dot, Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_BLAS, only: tMatVecFunc, tBiMatVecFunc

use StormRuler_Solvers, only: LinSolve
use StormRuler_ConvParams, only: tConvParams

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_Newton
  module procedure Solve_Newton
end interface Solve_Newton

interface Solve_JFNK
  module procedure Solve_JFNK
end interface Solve_JFNK

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a nonlinear operator equation: 𝓐(𝒙) = 𝒃, using the Newton's 
!! method.
!!
!! The classical Newton iterations are based on the linearization 
!! of 𝓐(𝒙) near 𝒙: 
!!
!! 𝓐(𝒙̂) ≈ 𝓐(𝒙) + [∂𝓐(𝒙)/∂𝒙](𝒙̂ - 𝒙) = 𝒃, 
!!
!! or, alternatively:
!!
!! [∂𝓐(𝒙)/∂𝒙]𝒕 = 𝒓, 𝒕 = 𝒙̂ - 𝒙, 𝒓 = 𝒃 - 𝓐(𝒙)
!!
!! where 𝒙 and 𝒙̂ are the current and updated solution vectors.
!! Therefore, a linear equation has to be solved on each iteration,
!! linear operator 𝓙(𝒙) ≈ ∂𝓐(𝒙)/∂𝒙 for computing Jacobian-vector 
!! products is required.
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
    ! Solve the Jacobian equation (using the current residual as the initial guess):
    ! 𝒕 ← 𝒓,
    ! 𝒕 ← 𝓙(𝒙)⁻¹𝒓,
    ! 𝒙 ← 𝒙 + 𝒕.
    ! ----------------------
    call Set(mesh, tArr, rArr)
    !! TODO: equation parameters!
    call jacConvParams%Init(1e-8_dp, 1e-8_dp, 2000, 'Newton')
    call LinSolve(mesh, 'GMRES', '', tArr, rArr, JacobianMatVecWithX, jacConvParams)
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
!! Solve a nonlinear operator equation: 𝓐(𝒙) = 𝒃, using the first
!! order Jacobian free-Newton-Krylov method (JFNK).
!!
!! For the Newton iterations, computing of the Jacobian-vector
!! products 𝒛 = 𝓙(𝒙)𝒚, where 𝓙(𝒙) ≈ ∂𝓐(𝒙)/∂𝒙 is required.
!! Consider the expansion:
!!
!! 𝓐(𝒙 + 𝛿⋅𝒚) = 𝓐(𝒙) + 𝛿⋅[∂𝓐(𝒙)/∂𝒙]𝒚 + 𝓞(𝛿²),
!!
!! where 𝛿 is some small number. Therefore,
!!
!! 𝓙(𝒙)𝒚 = [𝓐(𝒙 + 𝛿⋅𝒚) - 𝓐(𝒙)]/𝛿 = [∂𝓐(𝒙)/∂𝒙]𝒚 + 𝓞(𝛿).
!!
!! Expression above may be used as the formula for computing
!! the (approximate) Jacobian-vector products. Parameter 𝛿 is commonly 
!! defined as [1]:
!!
!! 𝛿 = 𝜇⋅‖𝒚‖⁺, 𝜇 = (𝜀ₘ)¹ᐟ²⋅(1 + ‖𝒙‖)¹ᐟ²,
!!
!! where 𝜀ₘ is the machine roundoff, ‖𝒚‖⁺ is the pseudo-inverse to ‖𝒚‖.
!!
!! References:
!! [1] Liu, Wei, Lilun Zhang, Ying Zhong, Yongxian Wang, 
!!     Yonggang Che, Chuanfu Xu and Xinghua Cheng. 
!!     “CFD High-order Accurate Scheme JFNK Solver.” 
!!     Procedia Engineering 61 (2013): 9-15.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_JFNK(mesh, MatVec, xArr, bArr, params)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: bArr
  class(tArray), intent(inout) :: xArr
  class(tConvParams), intent(inout) :: params
  procedure(tMatVecFunc) :: MatVec

  real(dp) :: delta, mu
  type(tArray) :: sArr, tArr, rArr, wArr
  type(tConvParams) :: jacConvParams 

  call AllocArray(sArr, tArr, rArr, wArr, mold=xArr)

  do
    ! ----------------------
    ! Compute residual:
    ! 𝒘 ← 𝓐(𝒙),
    ! 𝒓 ← 𝒃 - 𝒘,
    ! Check convergence for ‖𝒓‖.
    ! ----------------------
    call MatVec(mesh, wArr, xArr)
    call Sub(mesh, rArr, bArr, wArr)
    if (params%Check(Norm_2(mesh, rArr))) exit

    ! ----------------------
    ! Solve the Jacobian equation:
    ! 𝜇 ← (𝜀ₘ)¹ᐟ²⋅(1 + ‖𝒙‖)]¹ᐟ²,
    ! 𝒕 ← 𝒓,
    ! 𝒕 ← 𝓙(𝒙)⁻¹𝒓,
    ! 𝒙 ← 𝒙 + 𝒕.
    ! ----------------------
    mu = sqrt(epsilon(mu))*sqrt((1.0_dp + Norm_2(mesh, xArr)))
    call Set(mesh, tArr, rArr)
    !! TODO: equation parameters!
    call jacConvParams%Init(1e-8_dp, 1e-8_dp, 2000, 'Newton')
    call LinSolve(mesh, 'GMRES', '', tArr, rArr, ApproxJacobianMatVecWithX, jacConvParams)
    call Add(mesh, xArr, xArr, tArr)

  end do

contains
  subroutine ApproxJacobianMatVecWithX(mesh, zArr, yArr)
    class(tMesh), intent(in), target :: mesh
    class(tArray), intent(inout), target :: yArr, zArr

    ! ----------------------
    ! Compute the Jacobian-vector product:
    ! 𝛿 ← 𝜇⋅‖𝒚‖⁺,
    ! 𝒔 ← 𝒙 + 𝛿⋅𝒚,
    ! 𝒛 ← 𝓐(𝒔),
    ! 𝛿 ← 𝛿⁺,
    ! 𝒛 ← 𝛿⋅𝒛 - 𝛿⋅𝒘.
    ! ----------------------
    delta = SafeDivide(mu, Norm_2(mesh, yArr))
    call Add(mesh, sArr, xArr, yArr, delta)
    call MatVec(mesh, zArr, sArr)
    delta = SafeInverse(delta)
    call Sub(mesh, zArr, zArr, wArr, delta, delta)

  end subroutine ApproxJacobianMatVecWithX
end subroutine Solve_JFNK

end module StormRuler_Solvers_Newton
