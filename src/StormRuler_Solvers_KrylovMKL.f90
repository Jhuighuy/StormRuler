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
module StormRuler_Solvers_KrylovMKL

#$use 'StormRuler_Params.fi'
#$if HAS_MKL

use StormRuler_Parameters, only: dp, ip
use StormRuler_Mesh, only: tMesh
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Solvers_Krylov, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

include 'mkl_rci.fi'

integer(ip), parameter :: &
  & gMKL_RCI_DebugLevel = 1_ip

integer(ip), parameter :: &
  & gFGMRES_MKL_MaxNumNonRestartedIterations = 150_ip

interface Solve_CG_MKL
#$do rank = 0, NUM_RANKS
  module procedure Solve_CG_MKL$rank
#$end do
end interface Solve_CG_MKL

interface Solve_FGMRES_MKL
#$do rank = 0, NUM_RANKS
  module procedure Solve_FGMRES_MKL$rank
#$end do
end interface Solve_FGMRES_MKL

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear self-adjoint definite 
!! operator equation: Au = b, using the MKL CG RCI solver.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_CG_MKL$rank(mesh, u, b, MatVec, env, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(in) :: env
  type(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: rci_request
  integer(ip) :: n, iparams(0:127)
  real(dp) :: dparams(0:127)
  real(dp), allocatable, target :: tmp(:)
  real(dp), pointer :: in(@:,:), out(@:,:)

  ! ----------------------
  ! Preallocate storage.
  ! ----------------------
  n = size(u) 
  allocate(tmp(0:(3*n - 1)))
  call c_f_pointer( &
    & cptr=c_loc(tmp(0)), fptr=in, shape=shape(u))
  call c_f_pointer( &
    & cptr=c_loc(tmp(n)), fptr=out, shape=shape(b))

  ! ----------------------
  ! Initialize and configure MKL CG.
  ! ----------------------
  call dcg_init(n, u, b, rci_request, iparams, dparams, tmp)
  if (rci_request /= 0) then
    write(error_unit, *) 'MKL DCG_INIT FAILED'
    error stop 1
  end if
  associate( &
    &       outputWarningsAndErrors => iparams(1), &
    &              maxNumIterations => iparams(4), &
    & enableStopAtMaxIterationsTest => iparams(7), &
    &            enableResidualTest => iparams(8), &
    &         relativeTolerance_sqr => dparams(0), &
    &         absoluteTolerance_sqr => dparams(1), &
    &              enableCustomTest => iparams(9), &
    &         enablePreconditioning => iparams(10) )
    outputWarningsAndErrors = 6
    if (params%MaxNumIterations /= 0) then
      enableStopAtMaxIterationsTest = 1
      maxNumIterations = params%MaxNumIterations
    end if
    if (params%RelativeTolerance * &
      & params%AbsoluteTolerance > 0.0_dp) then
      enableResidualTest = 1
      relativeTolerance_sqr = params%RelativeTolerance
      absoluteTolerance_sqr = params%AbsoluteTolerance
    end if
    enableCustomTest = 0
    enablePreconditioning = 0
  end associate
  call dcg_check(n, u, b, rci_request, iparams, dparams, tmp)
  if (rci_request /= 0) then
    write(error_unit, *) 'MKL DCG_CHECK FAILED'
    error stop 1
  end if

  ! ----------------------
  ! Iterate MKL CG.
  ! ----------------------
  associate( &
    & initialResidualNorm_sqr => dparams(2), &
    & currentResidualNorm_sqr => dparams(4) )
    do
      call dcg(n, u, b, rci_request, iparams, dparams, tmp)
      if (rci_request == 0) exit
      if (rci_request == 1) then
        if (gMKL_RCI_DebugLevel > 0_ip) then
          write(error_unit, *) &
            & 'AE=', currentResidualNorm_sqr, &
            & 'RE=', currentResidualNorm_sqr/initialResidualNorm_sqr
        end if
        call MatVec(mesh, out, in, env)
      else
        write(error_unit, *) 'MKL DCG FAILED'
        error stop 1
      end if
    end do
  end associate
end subroutine Solve_CG_MKL$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: Au = b, 
!! using the MKL FGMRES RCI solver.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine Solve_FGMRES_MKL$rank(mesh, u, b, MatVec, env, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(in) :: env
  type(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: rci_request
  integer(ip) :: n, iparams(0:127), itercount
  real(dp) :: dparams(0:127)
  real(dp), allocatable, target :: tmp(:)
  real(dp), pointer :: in(@:,:), out(@:,:)

  ! ----------------------
  ! Preallocate storage.
  ! ----------------------
  ! Formula is used:
  ! (2*ipar[14] + 1)*n + ipar[14]*(ipar[14] + 9)/2 + 1)
  n = size(u)
  associate(k => gFGMRES_MKL_MaxNumNonRestartedIterations)
    associate(m => (2_ip*k + 1_ip)*n + k*(k + 9_ip)/2_ip + 1_ip)
      allocate(tmp(0_ip:(m-1_ip)))
    end associate
  end associate

  ! ----------------------
  ! Initialize and configure MKL FGMRES.
  ! ----------------------
  call dfgmres_init(n, u, b, rci_request, iparams, dparams, tmp)
  if (rci_request /= 0) then
    write(error_unit, *) &
      & 'MKL DFGMRES_INIT FAILED, RCI_REQUEST=', rci_request
    error stop 1
  end if
  associate( &
    &       outputWarningsAndErrors => iparams( 1), &
    &              maxNumIterations => iparams( 4), &
    & enableStopAtMaxIterationsTest => iparams( 7), &
    &            enableResidualTest => iparams( 8), &
    &         relativeTolerance_sqr => dparams( 0), &
    &         absoluteTolerance_sqr => dparams( 1), &
    &              enableCustomTest => iparams( 9), &
    &         enableCGVZeroNormTest => iparams(11), &
    &      cgvZeroNormTolerance_sqr => dparams( 7), &
    &         enablePreconditioning => iparams(10), &
    &   dfgmres_getRoutineBehaviour => iparams(12), &
    &     numNonRestartedIterations => iparams(14) )
    outputWarningsAndErrors = 6
    if (params%MaxNumIterations /= 0) then
      enableStopAtMaxIterationsTest = 1
      maxNumIterations = params%MaxNumIterations
    end if
    if (params%RelativeTolerance * &
      & params%AbsoluteTolerance > 0.0_dp) then
      enableResidualTest = 1
      relativeTolerance_sqr = params%RelativeTolerance
      absoluteTolerance_sqr = params%AbsoluteTolerance
    end if
    enablePreconditioning = 0
    enableCustomTest = 0
    ! What the fuck is the 'currently generated vector'?
    enableCGVZeroNormTest = 1
    cgvZeroNormTolerance_sqr = 1.0e-12_dp
    ! value = 0: update inplace,
    ! value > 0: write to RHS vector,
    ! value < 0: only get the iteration number.
    dfgmres_getRoutineBehaviour = 0
    numNonRestartedIterations = &
      & gFGMRES_MKL_MaxNumNonRestartedIterations
  end associate
  call dfgmres_check(n, u, b, rci_request, iparams, dparams, tmp)
  if (rci_request == -1100) then
    write(error_unit, *) &
      & 'MKL DFGMRES_CHECK FAILED, RCI_REQUEST=', rci_request
    error stop 1
  else if ((rci_request /= 0).and.(gMKL_RCI_DebugLevel > 0_ip)) then
    write(error_unit, *) &
      & 'MKL DFGMRES_CHECK ALTERED PARAMETERS, RCI_REQUEST=', rci_request
  end if

  ! ----------------------
  ! Iterate MKL FGMRES.
  ! ----------------------
  associate( &
    & initialResidualNorm_sqr => dparams(2), &
    & currentResidualNorm_sqr => dparams(4) )
    do
      call dfgmres(n, u, b, rci_request, iparams, dparams, tmp)
      if (rci_request == 0) exit
      if (rci_request == 1) then
        associate(inOffset => (iparams(21) - 1), &
          &      outOffset => (iparams(22) - 1))
          if (gMKL_RCI_DebugLevel > 1_ip) then
            write(error_unit, *) &
              & 'IN/OUT=', inOffset/n, outOffset/n
          end if
          call c_f_pointer( &
            & cptr=c_loc(tmp(inOffset)), fptr=in, shape=shape(u))
          call c_f_pointer( &
            & cptr=c_loc(tmp(outOffset)), fptr=out, shape=shape(b))
        end associate
        if (gMKL_RCI_DebugLevel > 0_ip) then
          write(error_unit, *) &
            & 'AE=', currentResidualNorm_sqr, &
            & 'RE=', currentResidualNorm_sqr/initialResidualNorm_sqr
        end if
        call MatVec(mesh, out, in, env)
      else
        write(error_unit, *) &
          & 'MKL DFGMRES FAILED, RCI_REQUEST=', rci_request
        error stop 1
      end if
    end do
  end associate

  ! ----------------------
  ! Get MKL FGMRES solution.
  ! ----------------------
  call dfgmres_get(n, u, b, rci_request, iparams, dparams, tmp, itercount)
  if (rci_request /= 0) then
    write(error_unit, *) &
      & 'MKL DFGMRES_GET FAILED, RCI_REQUEST=', rci_request
    error stop 1
  end if
  if (gMKL_RCI_DebugLevel > 1) then
    write(error_unit, *) 'ITERCOUNT=', itercount
  end if
end subroutine Solve_FGMRES_MKL$rank
#$end do

#$end if
end module StormRuler_Solvers_KrylovMKL
