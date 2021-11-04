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
use StormRuler_Array, only: tArrayR, AllocArray, FreeArray

use StormRuler_BLAS, only: Norm_2, Fill, Set, Sub 
#$for T, _ in [SCALAR_TYPES[0]]
use StormRuler_BLAS, only: tMatVecFunc$T
use StormRuler_Solvers_Precond, only: tPreMatVecFunc$T
#$end for

use StormRuler_ConvParams, only: tConvParams

use StormRuler_Solvers_CG, only: Solve_CG, Solve_BiCGStab
use StormRuler_Solvers_Chebyshev, only: Solve_Chebyshev
use StormRuler_Solvers_MINRES, only: Solve_MINRES, Solve_GMRES
use StormRuler_Solvers_LSQR, only: Solve_LSQR, Solve_LSMR

use StormRuler_Solvers_Precond, only: &
  & Precondition_Jacobi, Precondition_LU_SGS

use StormRuler_Threads, only: tThread, tSem

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
  class(tArrayR), intent(in) :: b
  class(tArrayR), intent(inout) :: x
  procedure(tMatVecFuncR) :: MatVec
  class(tConvParams), intent(inout) :: params
  character(len=*), intent(in) :: method, preMethod

  procedure(tMatVecFuncR), pointer :: uMatVec
  type(tArrayR) :: t, f
  
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
  if (Norm_2(mesh, t) == 0.0_dp) then
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
      write(error_unit, *) &
        & 'invalid precond method, preMethod=', preMethod
      error stop error_code
  end select

contains
  subroutine SelectMethod(PreMatVec)
    procedure(tPreMatVecFuncR), optional :: PreMatVec

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
      case('LSQR')
        params%Name = params%Name//'LSQR)'
        call Solve_LSQR(mesh, x, f, uMatVec, params, PreMatVec)
      case('LSMR')
        params%Name = params%Name//'LSMR)'
        call Solve_LSMR(mesh, x, f, uMatVec, params, PreMatVec)
      case default
        write(error_unit, *) 'invalid method, method=', method
        error stop error_code
    end select

  end subroutine SelectMethod
  subroutine MatVec_Uniformed(mesh, Ax, x)
    class(tMesh), intent(inout), target :: mesh
    class(tArrayR), intent(inout), target :: x, Ax

    call MatVec(mesh, Ax, x)
    call Sub(mesh, Ax, Ax, t)

  end subroutine MatVec_Uniformed
end subroutine LinSolve

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve a linear operator equation, using the 
!! Reverse Communication Interface (RCI): 
!! • self-adjoint definite operator case:
!!   [𝓜]𝓐[𝓜ᵀ]𝒚 = [𝓜]𝒃, [𝓜ᵀ]𝒚 = 𝒙, [𝓜𝓜ᵀ = 𝓟], 
!! • general nonsingular operator case:
!!   [𝓟]𝓐𝒙 = [𝓟]𝒃.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
function LinSolve_RCI(mesh, method, preMethod, &
    & x, b, params, Ay, y, resetState) result(request)
  class(tMesh), intent(inout), target, optional :: mesh
  character(len=*), intent(in), target, optional :: method, preMethod
  class(tArrayR), intent(in), target, optional :: b
  class(tArrayR), intent(inout), target, optional :: x
  class(tConvParams), intent(inout), target, optional :: params
  class(tArrayR), intent(out), pointer, optional :: Ay, y
  logical, intent(in), optional :: resetState
  character(len=:), allocatable :: request

  class(tMesh), pointer, save :: sMesh
  character(len=:), allocatable, save :: sMethod, sPrecondMethod
  class(tArrayR), pointer, save :: sX, sB
  class(tArrayR), pointer, save :: sAy, sY
  type(tConvParams), save :: sParams
  character(len=:), allocatable, save :: sRequest

  type(tThread), allocatable, save :: sThread
  type(tSem), allocatable, save :: sSemYield, sSemAwait

  ! ----------------------
  ! This is the reset request.
  ! ----------------------
  if (present(resetState)) then
    if (resetState) then
      if (allocated(sThread)) deallocate(sThread)
      if (allocated(sSemYield)) deallocate(sSemYield)
      if (allocated(sSemAwait)) deallocate(sSemAwait)
      if (allocated(sMethod)) deallocate(sMethod)
      if (allocated(sPrecondMethod)) deallocate(sPrecondMethod)
      if (allocated(sRequest)) deallocate(sRequest)
      request = 'Reset'
      write(error_unit, *) '`LinSolve_RCI` reset requested.'
      return
    end if
  end if

  if (.not.allocated(sRequest)) then
    ! ----------------------
    ! Save the arguments.
    ! ----------------------
    sMesh => mesh
    sMethod = method; sPrecondMethod = preMethod
    sX => x; sB => b
    sParams = params

    ! ----------------------
    ! Create the semaphores and launch the compute thread.
    ! ----------------------
    sSemYield = tSem(value=0)
    sSemAwait = tSem(value=0)
    sThread = tThread(StartRoutine)
  end if

  ! ----------------------
  ! Await for the compute thread to yield.
  ! ----------------------
  call sSemYield%Post()
  call sSemAwait%Wait()

  ! ----------------------
  ! Read the request.
  ! ----------------------
  request = sRequest
  Ay => sAy; y => sY

  ! ----------------------
  ! Clean-up on the exit request.
  ! ----------------------
  if (request == 'Done') then
    call sThread%Join()
    deallocate(sThread)
    deallocate(sSemYield, sSemAwait)
    deallocate(sMethod, sPrecondMethod, sRequest)
  end if

contains
  subroutine StartRoutine()

    print *, 'Entering the RCI thread function..'

    ! ----------------------
    ! Wait for the main thread to resume.
    ! ----------------------
    call sSemYield%Wait()

    ! ----------------------
    ! Enter the computational routine.
    ! ----------------------
    call LinSolve(sMesh, sMethod, sPrecondMethod, sX, sB, MatVec_RCI, sParams)

    ! ----------------------
    ! Pass the exit request to the main thread and leave.
    ! ----------------------
    sRequest = 'Done'
    call sSemAwait%Post()

    print *, 'Leaving the RCI thread function..'

  end subroutine StartRoutine
  subroutine MatVec_RCI(mesh, Ay, y)
    class(tMesh), intent(inout), target :: mesh
    class(tArrayR), intent(inout), target :: y, Ay

    ! ----------------------
    ! Pass the 'MatVec' or 'ConjMatVec' request 
    ! and field pointers to the main thread and yield.
    ! ----------------------
    sRequest = 'MatVec'
    sAy => Ay; sY => y

    call sSemAwait%Post()
    call sSemYield%Wait()

  end subroutine MatVec_RCI
end function LinSolve_RCI

end module StormRuler_Solvers_Unified
