!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module StormRuler_Helpers

use StormRuler_Parameters

implicit none  

abstract interface
  pure function MathFunc(x) result(f)
    import dp
    real(dp), intent(in) :: x
    real(dp) :: f
  end function MathFunc
end interface

contains

subroutine sort(n, a)
  integer :: n, i, j
  integer :: a(n), x
  do i = 2, n
    x = a(i)
    j = i - 1
    do while (j >= 1)
      if (a(j) <= x) exit
      a(j + 1) = a(j)
      j = j - 1
    end do
    a(j + 1) = x
  end do
end subroutine

!! -----------------------------------------------------------------
!! Integer increment.
subroutine Increment(value)
  integer, intent(inout) :: value
  value = value + 1
end subroutine Increment

!! -----------------------------------------------------------------  
!! Convert an integer to string.
function IntToString(value) result(string)
  integer, intent(in) :: value
  character(Len=20) :: string
  write(string, *) value
  string = adjustl(string)
end function IntToString

!! -----------------------------------------------------------------  
subroutine EnsurePositive(value)
  real(dp), intent(in) :: value
  if (isnan(value) .or. (value <= 0)) then
    print *, 'EnsurePositive', value
    error stop 1
  end if
end subroutine EnsurePositive

!! -----------------------------------------------------------------  
subroutine EnsureNonNegative(value)
  real(8), intent(in) :: value
  if (isnan(value) .or. (value < 0)) then
    print *, 'EnsureNonNegative', value
    error stop 1
  end if
end subroutine EnsureNonNegative

function SafeDivide(a,b) result(c)
  real(8), intent(in) :: a, b
  real(8) :: c 
  Call EnsureNonNegative(abs(a))
  Call EnsurePositive(abs(b))
  c = a/b
end function SafeDivide

elemental function SafeInverse(a) result(c)
  real(dp), intent(in) :: a
  real(dp) :: c 
  c = 0.0_dp
  if (a/=0.0_dp) c = 1.0_dp/a
end function SafeInverse

end module StormRuler_Helpers
