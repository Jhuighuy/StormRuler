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
module StormRuler_Solvers_MKL

#$use 'StormRuler_Params.fi'
#$if HAS_MKL or True

use StormRuler_Parameters, only: dp, ip
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: Set
#$for type_, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$type_
use StormRuler_Solvers_Precond, only: tPrecondFunc$type_
#$end for
use StormRuler_ConvParams, only: tConvParams

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

include 'mkl_rci.fi'

integer(ip), parameter :: &
  & gMKL_RCI_DebugLevel = 2_ip

integer(ip), parameter :: &
  & gFGMRES_MKL_MaxNumNonRestartedIterations = 150_ip

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint definite 
!! operator equation: [𝓜]𝓐[𝓜ᵀ]𝒚 = [𝓜]𝒃, [𝓜ᵀ]𝒚 = 𝒙, [𝓜𝓜ᵀ = 𝓟], 
!! using the MKL CG RCI solver.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_CG_MKL(mesh, u, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:)
  real(dp), intent(inout) :: u(:,:)
  procedure(tMatVecFuncR) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: rci_request
  integer(ip) :: n, iparams(128)
  real(dp) :: dparams(128)
  real(dp), allocatable :: tmp(:,:,:)
  class(*), allocatable :: precond_env

  ! ----------------------
  ! Preallocate storage.
  ! ----------------------
  n = size(u)
  allocate( tmp(size(u, dim=1), size(u, dim=2), merge(4, 3, present(Precond))) )

  ! ----------------------
  ! Initialize and configure MKL CG.
  ! ----------------------
  call dcg_init(n, u, b, rci_request, iparams, dparams, tmp)
  if (rci_request /= 0) then
    write(error_unit, *) 'MKL dcg_init failed, RCI_REQUEST=', rci_request
    error stop 1
  end if
  associate( &
    & outputWarningsAndErrors => iparams( 2), &
    & enableNumIterationsTest => iparams( 8), &
    &        enableCustomTest => iparams(10), &
    &   enablePreconditioning => iparams(11) )
    outputWarningsAndErrors = 6
    enableNumIterationsTest = 0
    enableCustomTest = 1
    enablePreconditioning = merge(1, 0, present(Precond))
  end associate
  call dcg_check(n, u, b, rci_request, iparams, dparams, tmp)
  if (rci_request /= 0) then
    write(error_unit, *) 'MKL dcg_check failed, RCI_REQUEST=', rci_request
    error stop 1
  end if

  ! ----------------------
  ! Iterate MKL CG.
  ! ----------------------
  do
    call dcg(n, u, b, rci_request, iparams, dparams, tmp)
    if (rci_request == 1) then
      call MatVec(mesh, tmp(:,:,2), tmp(:,:,1), env)
      cycle
    end if
    if (rci_request == 2) then
      associate( &
        & initialResidualNorm => sqrt(dparams(3)), &
        & currentResidualNorm => sqrt(dparams(5)) )
        if (params%Check(currentResidualNorm, &
          & currentResidualNorm/initialResidualNorm)) exit
      end associate
      cycle
    end if
    if (rci_request == 3) then
      call Precond(mesh, tmp(:,:,4), tmp(:,:,3), MatVec, env, precond_env)
      cycle
    end if
    write(error_unit, *) 'MKL dcg failed, RCI_REQUEST=', rci_request
    error stop 1
  end do
end subroutine Solve_CG_MKL

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: [𝓟]𝓐𝒙 = [𝓟]𝒃, 
!! using the MKL FGMRES RCI solver.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_FGMRES_MKL(mesh, u, b, MatVec, env, params, Precond)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:)
  real(dp), intent(inout) :: u(:,:)
  procedure(tMatVecFuncR) :: MatVec
  class(*), intent(inout) :: env
  class(tConvParams), intent(inout) :: params
  procedure(tPrecondFuncR), optional :: Precond
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: rci_request
  integer(ip) :: n, iparams(128), itercount
  real(dp) :: dparams(128)
  real(dp), pointer :: tmp(:), in(:,:), out(:,:)
  real(dp) :: trueInitialResidualNorm
  class(*), allocatable :: precond_env

  ! ----------------------
  ! Preallocate storage.
  ! ----------------------
  ! Formula is used:
  ! (2*ipar(15) + 1)*n + ipar(15)*(ipar(15) + 9)/2 + 1)
  n = size(u)
  associate(k => gFGMRES_MKL_MaxNumNonRestartedIterations)
    allocate(tmp((2*k + 1)*n + k*(k + 9)/2 + 1))
  end associate
  trueInitialResidualNorm = 0.0_dp

  ! ----------------------
  ! Initialize and configure MKL FGMRES.
  ! ----------------------
  call dfgmres_init(n, u, b, rci_request, iparams, dparams, tmp)
  if (rci_request /= 0) then
    write(error_unit, *) 'MKL dfgmres_init failed, RCI_REQUEST=', rci_request
    error stop 1
  end if
  associate( &
    &     outputWarningsAndErrors => iparams( 2), &
    &     enableNumIterationsTest => iparams( 8), &
    &            enableCustomTest => iparams(10), &
    &       enableCGVZeroNormTest => iparams(12), &
    &       enablePreconditioning => iparams(11), &
    & dfgmres_getRoutineBehaviour => iparams(13), &
    &   numNonRestartedIterations => iparams(15) )
    outputWarningsAndErrors = 6
    enableNumIterationsTest = 0
    enableCustomTest = 1
    ! CGV stands for 'currently generated vector'.
    enableCGVZeroNormTest = 1
    enablePreconditioning = merge(1, 0, present(Precond))
    ! value = 0: update inplace,
    ! value > 0: write to RHS vector,
    ! value < 0: only get the iteration number.
    dfgmres_getRoutineBehaviour = 0
    numNonRestartedIterations = gFGMRES_MKL_MaxNumNonRestartedIterations
  end associate
  call dfgmres_check(n, u, b, rci_request, iparams, dparams, tmp)
  if (rci_request == -1100) then
    write(error_unit, *) 'MKL dfgmres_check failed, RCI_REQUEST=', rci_request
    error stop 1
  else if ((rci_request /= 0).and.(gMKL_RCI_DebugLevel > 0_ip)) then
    write(error_unit, *) 'MKL dfgmres_check altered parameters, RCI_REQUEST=', rci_request
  end if

  ! ----------------------
  ! Iterate MKL FGMRES.
  ! ----------------------
  do
    call dfgmres(n, u, b, rci_request, iparams, dparams, tmp)
    if (rci_request == 1) then
      associate(inOffset => iparams(22), outOffset => iparams(23))
        call c_f_pointer(cptr=c_loc(tmp(inOffset)), fptr=in, shape=shape(u))
        call c_f_pointer(cptr=c_loc(tmp(outOffset)), fptr=out, shape=shape(u))
      end associate
      call MatVec(mesh, out, in, env)
      cycle
    end if
    if (rci_request == 2) then
      associate(initialResidualNorm => sqrt(dparams(3)), &
              & currentResidualNorm => sqrt(dparams(5)))
        ! Initial residual norm is reset after the restart, so we
        ! should keep the true one for the correct relative error test.
        if (trueInitialResidualNorm == 0.0_dp) &
          & trueInitialResidualNorm = initialResidualNorm
        if (params%Check(currentResidualNorm, &
          & currentResidualNorm/trueInitialResidualNorm)) exit
      end associate
      cycle
    end if
    if (rci_request == 3) then
      associate(inOffset => iparams(22), outOffset => iparams(23))
        call c_f_pointer(cptr=c_loc(tmp(inOffset)), fptr=in, shape=shape(u))
        call c_f_pointer(cptr=c_loc(tmp(outOffset)), fptr=out, shape=shape(u))
      end associate
      call Precond(mesh, out, in, MatVec, env, precond_env)
      cycle
    end if
    write(error_unit, *) 'MKL dfgmres failed, RCI_REQUEST=', rci_request
    error stop 1
  end do

  ! ----------------------
  ! Get MKL FGMRES solution.
  ! ----------------------
  call dfgmres_get(n, u, b, rci_request, iparams, dparams, tmp, itercount)
  if (rci_request /= 0) then
    write(error_unit, *) 'MKL dfgmres_get failed, RCI_REQUEST=', rci_request
    error stop 1
  end if
end subroutine Solve_FGMRES_MKL

#$end if
end module StormRuler_Solvers_MKL
