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
module StormRuler_ConvParams

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Helpers, only: EnsurePositive, EnsureNonNegative

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! ----------------------------------------------------------------- !!
!! A class that controls convergence for the iterative algorithms. 
!! ----------------------------------------------------------------- !!
type :: tConvParams
  character(len=:), allocatable :: Name

  ! ----------------------
  ! Number of iterations counter.
  ! ----------------------
  integer(ip) :: NumIterations
  integer(ip) :: MaxNumIterations

  ! ----------------------
  ! Convergence tolerances.
  ! ----------------------
  real(dp) :: AbsoluteTolerance
  real(dp) :: RelativeTolerance

contains
  procedure :: Init => tConvParams_Init
  generic :: Check => CheckR, CheckC
  procedure :: CheckR => tConvParams_Check
  procedure :: CheckC => tConvParams_CheckC
end type tConvParams

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Initialize iteration parameters.
!! ----------------------------------------------------------------- !!
subroutine tConvParams_Init(params, &
    & absoluteTolerance, relativeTolerance, maxNumIterations, name)
  class(tConvParams), intent(out) :: params
  real(dp), intent(in) :: absoluteTolerance
  real(dp), intent(in), optional :: relativeTolerance
  integer(ip), intent(in), optional :: maxNumIterations
  character(len=*), intent(in), optional :: name

  params%Name = ''
  if (present(name)) params%Name = name

  ! ----------------------
  ! Initialize errors.
  ! ----------------------
  call EnsurePositive(absoluteTolerance)
  params%AbsoluteTolerance = absoluteTolerance
  if (present(relativeTolerance)) then
    call EnsurePositive(relativeTolerance)
    params%RelativeTolerance = relativeTolerance
  end if

  ! ----------------------
  ! Initialize iterations.
  ! ----------------------
  params%NumIterations = 0
  if (present(maxNumIterations)) then
    call EnsurePositive(relativeTolerance)
    params%MaxNumIterations = maxNumIterations
  end if
  !print *, 'Init iterations:', & 
  !  & params%AbsoluteTolerance, &
  !  & params%RelativeTolerance, params%MaxNumIterations
end subroutine tConvParams_Init

!! ----------------------------------------------------------------- !!
!! Check convergence of the iterations process.  
!! ----------------------------------------------------------------- !!
logical function tConvParams_Check(params, absoluteError, relativeError)
  class(tConvParams), intent(inout) :: params
  real(dp), intent(in) :: absoluteError
  real(dp), intent(in), optional :: relativeError

  tConvParams_Check = .false.

  ! ----------------------
  ! Check whether number of iterations has exceeded.
  ! ----------------------
  associate(numIterations => params%NumIterations, &
      &  maxNumIterations => params%MaxNumIterations)
    numIterations = numIterations + 1
    if ((maxNumIterations /= 0).and.(numIterations >= maxNumIterations)) then
      tConvParams_Check = .true.; return
    end if
  end associate
  
  ! ----------------------
  ! Check convergence by absolute and relative errors.
  ! ----------------------
  associate(absoluteTolerance => params%AbsoluteTolerance, &
    &       relativeTolerance => params%RelativeTolerance)
    call EnsureNonNegative(absoluteError)
    if (absoluteError <= absoluteTolerance) then
      print *, '----------------', params%Name, params%NumIterations
      tConvParams_Check = .true.; return
    end if
    if (present(relativeError).and.(relativeTolerance > 0)) then
      call EnsureNonNegative(relativeError)
      print *, params%NumIterations, relativeError, absoluteError
      if (relativeError <= relativeTolerance) then
        print *, '----------------', params%Name, params%NumIterations
        tConvParams_Check = .true.; return
      end if
    end if
  end associate
end function tConvParams_Check
logical function tConvParams_CheckC(params, absoluteError, relativeError)
  class(tConvParams), intent(inout) :: params
  complex(dp), intent(in) :: absoluteError
  complex(dp), intent(in), optional :: relativeError

  if (present(relativeError)) then
    tConvParams_CheckC = params%CheckR(abs(absoluteError), abs(relativeError))
  else
    tConvParams_CheckC = params%CheckR(abs(absoluteError))
  end if
end function tConvParams_CheckC

end module StormRuler_ConvParams
