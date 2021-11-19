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
module StormRuler_Solvers_Unified

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8, error_code

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: Norm_2, Fill, Set, Sub 
use StormRuler_BLAS, only: tMatVecFunc
use StormRuler_Solvers_Precond, only: tPreMatVecFunc

use StormRuler_ConvParams, only: tConvParams

use StormRuler_Solvers_CG, only: Solve_CG, Solve_BiCGStab
use StormRuler_Solvers_MINRES, only: &
  & Solve_MINRES, Solve_GMRES, Solve_QMR, Solve_TFQMR
use StormRuler_Solvers_LSQR, only: Solve_LSQR, Solve_LSMR

use StormRuler_Solvers_Precond, only: &
  & Precondition_Jacobi, Precondition_LU_SGS

use, intrinsic :: iso_fortran_env, only: error_unit

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve a linear operator equation: 
!! • self-adjoint definite operator case:
!!   [𝓜]𝓐[𝓜ᵀ]𝒚 = [𝓜]𝒃, [𝓜ᵀ]𝒚 = 𝒙, [𝓜𝓜ᵀ = 𝓟], 
!! • general nonsingular operator case:
!!   [𝓟]𝓐𝒙 = [𝓟]𝒃.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine LinSolve(mesh, method, preMethod, x, b, MatVec, params)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: b
  class(tArray), intent(inout) :: x
  procedure(tMatVecFunc) :: MatVec
  class(tConvParams), intent(inout) :: params
  character(len=*), intent(in) :: method, preMethod

  procedure(tMatVecFunc), pointer :: uMatVec
  type(tArray) :: t, f
  
  ! ----------------------
  ! Check if the operator is non-uniform, e.g. 𝓐(𝒙) = 𝓐𝒙 + 𝒕:
  ! 𝒕 ← 𝓐(0),
  ! 𝗶𝗳 𝒕 = 0:
  !   𝒇 ← 𝒕,
  ! 𝗲𝗹𝘀𝗲:
  !   𝒇 ← 𝒃 - 𝒕,
  !   𝓐(𝒙) ← 𝓐(𝒙) - 𝒕,
  ! 𝗲𝗻𝗱 𝗶𝗳,
  ! 𝘀𝗼𝗹𝘃𝗲: 𝓐(𝒙) = 𝒇.
  ! ----------------------
  call AllocArray(t, f, mold=x)
  call Fill(mesh, f, 0.0_dp)
  call MatVec(mesh, t, f)
  if (Norm_2(mesh, t) <= epsilon(1.0_dp)) then
    uMatVec => MatVec
    call FreeArray(t, f); f = b
  else
    uMatVec => MatVec_Uniformed
    call Sub(mesh, f, b, t)
  end if

  ! ----------------------
  ! Two-step call is utilized 
  ! in order to match optional preconditioning.
  ! ----------------------
  select case(preMethod)
    case('', 'none')
      params%Name = params%Name//'('
      call SelectMethod()
    case('Jacobi')
      params%Name = params%Name//'(Jacobi-'
      call SelectMethod(Precondition_Jacobi)
    case('LU_SGS')
      params%Name = params%Name//'(LU-SGS-'
      call SelectMethod(Precondition_LU_SGS)
    case default
      error stop 'invalid precond method, preMethod='//preMethod
  end select

contains
  subroutine SelectMethod(PreMatVec)
    procedure(tPreMatVecFunc), optional :: PreMatVec

    select case(method)
      case('CG')
        params%Name = params%Name//'CG)'
        call Solve_CG(mesh, x, f, uMatVec, params, PreMatVec)
      case('BiCGStab')
        params%Name = params%Name//'BiCGStab)'
        call Solve_BiCGStab(mesh, x, f, uMatVec, params, PreMatVec)
      case('MINRES')
        params%Name = params%Name//'MINRES)'
        call Solve_MINRES(mesh, x, f, uMatVec, params, PreMatVec)
      case('GMRES')
        params%Name = params%Name//'GMRES)'
        call Solve_GMRES(mesh, x, f, uMatVec, params, PreMatVec)
      case('QMR')
        params%Name = params%Name//'QMR)'
        call Solve_QMR(mesh, x, f, uMatVec, params, PreMatVec)
      case('TFQMR')
        params%Name = params%Name//'TFQMR)'
        call Solve_QMR(mesh, x, f, uMatVec, params, PreMatVec)
      case('LSQR')
        params%Name = params%Name//'LSQR)'
        call Solve_LSQR(mesh, x, f, uMatVec, params, PreMatVec)
      case('LSMR')
        params%Name = params%Name//'LSMR)'
        call Solve_LSMR(mesh, x, f, uMatVec, params, PreMatVec)
      case default
        error stop 'invalid method, method='//method
    end select

  end subroutine SelectMethod
  subroutine MatVec_Uniformed(mesh, Ax, x)
    class(tMesh), intent(inout), target :: mesh
    class(tArray), intent(inout), target :: x, Ax

    call MatVec(mesh, Ax, x)
    call Sub(mesh, Ax, Ax, t)

  end subroutine MatVec_Uniformed
end subroutine LinSolve

end module StormRuler_Solvers_Unified
