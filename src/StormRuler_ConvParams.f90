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

use StormRuler_Parameters, only: dp
use StormRuler_Helpers, only: EnsurePositive, EnsureNonNegative

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! ----------------------------------------------------------------- !!
!! A class that controls convergence for the iterative algorithms. 
!! ----------------------------------------------------------------- !!
type :: tConvParams
  integer :: NumIterations
  integer :: MaxNumIterations
  real(dp) :: AbsoluteTolerance
  real(dp) :: RelativeTolerance
contains
  procedure :: Init => tConvParams_Init
  procedure :: Check => tConvParams_Check
end type tConvParams

private :: tConvParams_Init, tConvParams_Check

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Initialize iteration parameters.
!! ----------------------------------------------------------------- !!
subroutine tConvParams_Init(params, &
    & absoluteTolerance, relativeTolerance, maxNumIterations)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tConvParams), intent(out) :: params
  real(dp), intent(in) :: absoluteTolerance
  real(dp), intent(in), optional :: relativeTolerance
  integer, intent(in), optional :: maxNumIterations
  ! >>>>>>>>>>>>>>>>>>>>>>
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
  else
    params%MaxNumIterations = huge(params%MaxNumIterations)
  end if
  !print *, 'Init iterations:', & 
  !  & params%AbsoluteTolerance, &
  !  & params%RelativeTolerance, params%MaxNumIterations
end subroutine tConvParams_Init

!! ----------------------------------------------------------------- !!
!! Check convergence of the iterations process.  
!! ----------------------------------------------------------------- !!
logical function tConvParams_Check(params, absoluteError, relativeError)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tConvParams), intent(inout) :: params
  real(dp), intent(in) :: absoluteError
  real(dp), intent(in), optional :: relativeError
  ! >>>>>>>>>>>>>>>>>>>>>>
  tConvParams_Check = .false.
  ! ----------------------
  ! Check whether number of iterations has exceeded.
  ! ----------------------
  associate(numIterations=>params%NumIterations, &
      &  maxNumIterations=>params%MaxNumIterations)
    numIterations = numIterations + 1
    if (numIterations >= params%MaxNumIterations) then
      tConvParams_Check = .true.; return
    end if
  end associate
  ! ----------------------
  ! Check convergence by absolute and relative errors.
  ! ----------------------
  associate(absoluteTolerance => params%AbsoluteTolerance, &
    &       relativeTolerance => params%RelativeTolerance)
    call EnsureNonNegative(absoluteError)
    !print *, absoluteError
    if (absoluteError <= absoluteTolerance) then
      !print *, '----------------'
      tConvParams_Check = .true.; return
    end if
    if (present(relativeError).and.(relativeTolerance > 0)) then
      call EnsureNonNegative(relativeError)
      if (relativeError <= relativeTolerance) then
        !print *, '----------------'
        tConvParams_Check = .true.; return
      end if
    end if
  end associate
end function tConvParams_Check

end module StormRuler_ConvParams
