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
module StormRuler_Solvers

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray, FreeArray

use StormRuler_BLAS, only: Norm_2, Fill, Set, Sub 
use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_ConvParams, only: tConvParams
use StormRuler_Precond, only: tPreconditioner

use StormRuler_Matrix!, only: ...
use StormRuler_Matrix_Extraction!, only: ...
use StormRuler_Precond_SPAI!, only: ...

use StormRuler_Solvers_CG, only: Solve_CG, Solve_BiCGStab
use StormRuler_Solvers_MINRES, only: &
  & Solve_MINRES, Solve_GMRES, Solve_QMR, Solve_TFQMR
use StormRuler_Solvers_LSQR, only: Solve_LSQR, Solve_LSMR

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

include 'mkl.fi'

type, extends(tPreconditioner) :: tPreconditioner_ILU0_MKL
  type(tMatrix), pointer :: Mat => null()
  real(dp), allocatable :: ColCoeffsILU0(:)

contains
  procedure Init => InitPrecond_ILU0_MKL
  procedure Apply => ApplyPrecond_ILU0_MKL
end type tPreconditioner_ILU0_MKL

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

subroutine InitPrecond_ILU0_MKL(precond, mesh, MatVec)
  class(tPreconditioner_ILU0_MKL), intent(inout) :: precond
  class(tMesh), intent(inout), target :: mesh
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: iparam(128), ierror
  real(dp) :: dparam(128)

  iparam(:) = 0; dparam(:) = 0
  iparam(2) = 6; dparam(31) = 1.0e-16_dp 
  allocate(precond%ColCoeffsILU0(size(precond%Mat%ColCoeffs)))

  call dcsrilu0(mesh%NumCells, precond%Mat%ColCoeffs, precond%Mat%RowAddrs, &
    & precond%Mat%ColIndices, precond%ColCoeffsILU0, iparam, dparam, ierror)
  if (ierror /= 0) then
    error stop 'MKL `dcsrilu0` has failed, ierror='//I2S(ierror)
  end if

end subroutine InitPrecond_ILU0_MKL

subroutine ApplyPrecond_ILU0_MKL(precond, mesh, yArr, xArr, MatVec)
  class(tPreconditioner_ILU0_MKL), intent(inout) :: precond
  class(tMesh), intent(inout), target :: mesh
  class(tArray), intent(inout), target :: xArr, yArr
  procedure(tMatVecFunc) :: MatVec

  type(tArray) :: tArr
  real(dp), pointer :: t(:), x(:), y(:)

  call AllocArray(tArr, mold=xArr)

  ! solve (L*U)*y = L*(U*y) = x
  ! solve L*t = x
  ! solve U*y = t

  call tArr%Get(t); call xArr%Get(x); call yArr%Get(y)

  call mkl_dcsrtrsv('L', 'N', 'U', mesh%NumCells, &
    & precond%ColCoeffsILU0, precond%Mat%RowAddrs, precond%Mat%ColIndices, x, t)
  call mkl_dcsrtrsv('U', 'N', 'N', mesh%NumCells, &
    & precond%ColCoeffsILU0, precond%Mat%RowAddrs, precond%Mat%ColIndices, t, y)

end subroutine ApplyPrecond_ILU0_MKL

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

  type(tMatrix), target, save :: mat

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
  !select case(preMethod)
  !  case('', 'none')
      params%Name = params%Name//'('
      call SelectMethod()
  !  case('Jacobi')
  !    params%Name = params%Name//'(Jacobi-'
  !    call SelectMethod(Precondition_Jacobi)
  !  case('LU_SGS')
  !    params%Name = params%Name//'(LU-SGS-'
  !    call SelectMethod(Precondition_LU_SGS)
  !  case default
  !    error stop 'invalid precond method, preMethod='//preMethod
  !end select

contains
  subroutine SelectMethod(precond)
    class(tPreconditioner), intent(inout), optional :: precond
    
    if (preMethod == 'extr') then; block
      type(tMatrixLabeling), save :: labeling
      type(tPreconditioner_ILU0_MKL) :: ilu

      if (labeling%NumLabels == 0) then
        call InitMatrix(mesh, mat, 2)
        call LabelColumns_Patterned(mesh, mat, labeling)
      end if
      call ExtractMatrix(mesh, mat, labeling, uMatVec, mold=x)
      ilu%mat => mat

      params%Name = params%Name//'EXTR)'
      call Solve_BiCGStab(mesh, x, f, uMatVec, params, ilu)
      return

    end block; end if

    select case(method)
      case('CG')
        params%Name = params%Name//'CG)'
        call Solve_CG(mesh, x, f, uMatVec, params, precond)
      case('BiCGStab')
        params%Name = params%Name//'BiCGStab)'
        call Solve_BiCGStab(mesh, x, f, uMatVec, params, precond)
      case('MINRES')
        params%Name = params%Name//'MINRES)'
        call Solve_MINRES(mesh, x, f, uMatVec, params, precond)
      case('GMRES')
        params%Name = params%Name//'GMRES)'
        call Solve_GMRES(mesh, x, f, uMatVec, params, precond)
      case('QMR')
        params%Name = params%Name//'QMR)'
        call Solve_QMR(mesh, x, f, uMatVec, params, precond)
      case('TFQMR')
        params%Name = params%Name//'TFQMR)'
        call Solve_TFQMR(mesh, x, f, uMatVec, params, precond)
      case('LSQR')
        params%Name = params%Name//'LSQR)'
        call Solve_LSQR(mesh, x, f, uMatVec, params, precond)
      case('LSMR')
        params%Name = params%Name//'LSMR)'
        call Solve_LSMR(mesh, x, f, uMatVec, params, precond)
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

end module StormRuler_Solvers
