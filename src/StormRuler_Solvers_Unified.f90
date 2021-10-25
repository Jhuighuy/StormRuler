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
use StormRuler_BLAS, only: Norm_2, Fill, Set, Sub 
#$for type_, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$type_
#$end for

use StormRuler_ConvParams, only: tConvParams

use StormRuler_Solvers_CG, only: Solve_CG, Solve_BiCGStab
use StormRuler_Solvers_Chebyshev, only: Solve_Chebyshev
use StormRuler_Solvers_MINRES, only: Solve_MINRES
use StormRuler_Solvers_GMRES, only: Solve_GMRES
use StormRuler_Solvers_LSQR, only: Solve_LSQR, Solve_LSMR

#$for type_, _ in SCALAR_TYPES
use StormRuler_Solvers_Precond, only: tPrecondFunc$type_
#$end for
use StormRuler_Solvers_Precond, only: &
  & Precondition_Jacobi, Precondition_LU_SGS

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: c_int, c_long, &
  & c_size_t, c_ptr, c_funptr, c_null_ptr, c_funloc

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
subroutine LinSolve(mesh, method, precondMethod, x, b, MatVec, params)
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: method, precondMethod
  real(dp), intent(in), target :: b(:,:)
  real(dp), intent(inout) :: x(:,:)
  procedure(tMatVecFuncR) :: MatVec
  class(tConvParams), intent(inout) :: params

  procedure(tMatVecFuncR), pointer :: uMatVec
  real(dp), pointer :: t(:,:), f(:,:)
  
  ! ----------------------
  ! Check if the operator is non-uniform, e.g. 𝓐(𝒙) = 𝓐𝒙 + 𝒕:
  ! 𝒕 ← 𝓐(0),
  ! 𝗶𝗳 𝒕 = 0:
  !   𝒇 ← 𝒕,
  ! 𝗲𝗹𝘀𝗲:
  !   𝒇 ← 𝒃 - 𝒕,
  !   𝓐(𝒙) ← 𝓐(𝒙) - 𝒕,
  ! 𝗲𝗻𝗱 𝗶𝗳,
  ! 𝘀𝗼𝗹v𝗲: 𝓐(𝒙) = 𝒇.
  ! ----------------------
  allocate(t, f, mold=x)
  call Fill(mesh, f, 0.0_dp)
  call MatVec(mesh, t, f)
  if (Norm_2(mesh, t) == 0.0_dp) then
    uMatVec => MatVec
    deallocate(t, f); f => b
  else
    uMatVec => MatVec_Uniformed
    call Sub(mesh, f, b, t)
  end if

  ! ----------------------
  ! Two-step call is utilized 
  ! in order to match optional preconditioning.
  ! ----------------------
  select case(precondMethod)
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
        & 'invalid precond method, precondMethod=', precondMethod
      error stop error_code
  end select

contains
  subroutine SelectMethod(Precond)
    procedure(tPrecondFuncR), optional :: Precond

    select case(method)
      case('CG')
        params%Name = params%Name//'CG)'
        call Solve_CG(mesh, x, f, uMatVec, params, Precond)
      case('BiCGStab')
        params%Name = params%Name//'BiCGStab)'
        call Solve_BiCGStab(mesh, x, f, uMatVec, params, Precond)
      case('MINRES')
        params%Name = params%Name//'MINRES)'
        call Solve_MINRES(mesh, x, f, uMatVec, params, Precond)
      case('GMRES')
        params%Name = params%Name//'GMRES)'
        call Solve_GMRES(mesh, x, f, uMatVec, params, Precond)
      case('LSQR')
        params%Name = params%Name//'LSQR)'
        call Solve_LSQR(mesh, x, f, uMatVec, params, Precond)
      case('LSMR')
        params%Name = params%Name//'LSMR)'
        call Solve_LSMR(mesh, x, f, uMatVec, params, Precond)
      case default
        write(error_unit, *) 'invalid method, method=', method
        error stop error_code
    end select

  end subroutine SelectMethod
  subroutine MatVec_Uniformed(mesh, Ax, x)
    class(tMesh), intent(inout), target :: mesh
    real(dp), intent(inout), target :: x(:,:), Ax(:,:)

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
function LinSolve_RCI(mesh, method, precondMethod, x, b, params, Ay, y) result(request)
  class(tMesh), intent(inout), target :: mesh
  character(len=*), intent(in), target :: method, precondMethod
  real(dp), intent(in), target :: b(:,:)
  real(dp), intent(inout), target :: x(:,:)
  class(tConvParams), intent(inout), target :: params
  real(dp), intent(out), pointer :: Ay(:,:), y(:,:)
  character(len=:), allocatable :: request

  !! --------------------------------------------------------------- !!
  !! PThread API.
  !! --------------------------------------------------------------- !!
  type, bind(C) :: ctPThread
    integer(c_long) :: mPrivate
  end type ctPThread
  interface
    subroutine cPThreadCreate(pThread, pAttr, pFunc, env) bind(C, name='pthread_create')
      import :: c_ptr, c_funptr, ctPThread
      type(ctPThread), intent(in) :: pThread
      type(c_ptr), intent(in), value :: pAttr
      type(c_funptr), intent(in), value :: pFunc
      type(c_ptr), intent(in), value :: env
    end subroutine cPThreadCreate
    subroutine cPThreadJoin(pThread, pRetval) bind(C, name='pthread_join')
      import :: c_ptr, ctPThread
      type(ctPThread), intent(in), value :: pThread
      type(c_ptr), intent(in), value :: pRetval
    end subroutine cPThreadJoin
  end interface

  !! --------------------------------------------------------------- !!
  !! Unix semaphores API.
  !! --------------------------------------------------------------- !!
  type, bind(C) :: ctSem
    integer(c_size_t) :: mPrivate(4)
  end type ctSem
  interface
    subroutine cSemInit(pSem, shared, value) bind(C, name='sem_init')
      import :: c_int, ctSem
      type(ctSem), intent(in) :: pSem
      integer(c_int), intent(in), value :: shared, value
    end subroutine cSemInit
    subroutine cSemClose(pSem) bind(C, name='sem_close')
      import :: ctSem
      type(ctSem), intent(in) :: pSem
    end subroutine cSemClose
    subroutine cSemPost(pSem) bind(C, name='sem_post')
      import :: ctSem
      type(ctSem), intent(in) :: pSem
    end subroutine cSemPost
    subroutine cSemWait(pSem) bind(C, name='sem_wait')
      import :: ctSem
      type(ctSem), intent(in) :: pSem
    end subroutine cSemWait
  end interface

  class(tMesh), pointer, save :: sMesh
  character(len=:), allocatable, save :: sMethod, sPrecondMethod
  real(dp), pointer, save :: sX(:,:), sB(:,:)
  real(dp), pointer, save :: sAy(:,:), sY(:,:)
  type(tConvParams), save :: sParams
  character(len=:), allocatable, save :: sRequest

  logical, save :: sFirstCall = .true.
  type(ctPThread), save :: sThread
  type(ctSem), save :: sSemYield, sSemAwait

  if (sFirstCall) then
    ! ----------------------
    ! Save the arguments.
    ! ----------------------
    sMesh => mesh
    sMethod = method; sPrecondMethod = precondMethod
    sX => x; sB => b
    sParams = params

    ! ----------------------
    ! Create the semaphores and launch the compute thread.
    ! ----------------------
    sFirstCall = .false.
    call cSemInit(sSemYield, 0, 0)
    call cSemInit(sSemAwait, 0, 0)
    call cPThreadCreate(sThread, c_null_ptr, c_funloc(cThreadFunc), c_null_ptr)
  end if

  ! ----------------------
  ! Await for the compute thread to yield.
  ! ----------------------
  call cSemPost(sSemYield)
  call cSemWait(sSemAwait)

  ! ----------------------
  ! Read the request.
  ! ----------------------
  request = sRequest
  Ay => sAy; y => sY

  ! ----------------------
  ! Clean-up on the exit request.
  ! ----------------------
  if (request == 'Done') then
    sFirstCall = .true.
    call cPThreadJoin(sThread, c_null_ptr)
    call cSemClose(sSemYield)
    call cSemClose(sSemAwait)
    deallocate(sMethod, sPrecondMethod, sRequest)
  end if

contains
  subroutine cThreadFunc(env) bind(C)
    type(c_ptr), intent(in), value :: env

    print *, 'Entering the RCI thread function..'

    ! ----------------------
    ! Wait for the main thread to resume.
    ! ----------------------
    call cSemWait(sSemYield)

    ! ----------------------
    ! Enter the computational routine.
    ! ----------------------
    call LinSolve(sMesh, sMethod, sPrecondMethod, sX, sB, MatVec_RCI, sParams)

    ! ----------------------
    ! Pass the exit request to the main thread and leave.
    ! ----------------------
    sRequest = 'Done'
    call cSemPost(sSemAwait)

    print *, 'Leaving the RCI thread function..'

  end subroutine cThreadFunc
  subroutine MatVec_RCI(mesh, Ay, y)
    class(tMesh), intent(inout), target :: mesh
    real(dp), intent(inout), target :: y(:,:), Ay(:,:)

    ! ----------------------
    ! Pass the 'MatVec' or 'MatVec_H' request 
    ! and field pointers to the main thread and yield.
    ! ----------------------
    sRequest = 'MatVec'
    sAy => Ay; sY => y 

    call cSemPost(sSemAwait)
    call cSemWait(sSemYield)

  end subroutine MatVec_RCI
end function LinSolve_RCI

end module StormRuler_Solvers_Unified
