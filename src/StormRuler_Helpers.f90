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
#$use 'StormRuler_Parameters.f90'

implicit none  

interface operator(.inner.)
#$do rank = 0, NUM_RANKS-1
  module procedure Inner$rank
#$end do
end interface

interface operator(.outer.)
#$do rank = 0, NUM_RANKS-1
  module procedure Outer$rank
#$end do
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
!! Convert an integer to string.
function I2S(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer, intent(in) :: value
  character(len=:), allocatable :: I2S
  ! >>>>>>>>>>>>>>>>>>>>>>
  character(len=20) :: buffer
  write(buffer, *) value
  I2S = trim(adjustl(buffer))
end function I2S

!! -----------------------------------------------------------------  
!! Convert a real number to string.
function R2S(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: value
  character(len=:), allocatable :: R2S
  ! >>>>>>>>>>>>>>>>>>>>>>
  character(len=50) :: buffer
  write(buffer, *) value
  R2S = trim(adjustl(buffer))
end function R2S

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pure function Find(array,value) result(index)
  integer, intent(in) :: array(:), value
  integer :: index
  index = findloc(array,value,dim=1)
end function Find

function MergeString(trueString,falseString,condition) result(string)
  ! <<<<<<<<<<<<<<<<<<<<<<
  character(len=*), intent(in) :: trueString,falseString
  logical, intent(in) :: condition
  character(len=:), allocatable :: string
  ! >>>>>>>>>>>>>>>>>>>>>>
  if (condition)  then
    string = trueString
  else
    string = falseString
  end if
end function MergeString

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
subroutine EnsurePositive(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: value
  ! >>>>>>>>>>>>>>>>>>>>>>
  if (isnan(value).or.(value <= 0)) then
    print *, 'EnsurePositive', value
    error stop 1
  end if
end subroutine EnsurePositive

!! -----------------------------------------------------------------  
subroutine EnsureNonNegative(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(8), intent(in) :: value
  ! >>>>>>>>>>>>>>>>>>>>>>
  if (isnan(value).or.(value < 0)) then
    print *, 'EnsureNonNegative', value
    error stop 1
  end if
end subroutine EnsureNonNegative

!! -----------------------------------------------------------------  
function SafeDivide(a,b) result(c)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(8), intent(in) :: a, b
  real(8) :: c 
  ! >>>>>>>>>>>>>>>>>>>>>>
  call EnsureNonNegative(abs(a))
  call EnsurePositive(abs(b))
  c = a/b
end function SafeDivide

!! -----------------------------------------------------------------  
!! Compute pseudo inverse: i ← 1/a if a ≠ 0 else 0.
elemental function SafeInverse(a) result(i)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: a
  real(dp) :: i
  ! >>>>>>>>>>>>>>>>>>>>>>
  i = merge(0.0_dp, 1.0_dp/a, a==0.0_dp)
end function SafeInverse
!! -----------------------------------------------------------------  

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
pure function Inner$0(v,w) result(u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: v(:),w(:)
  real(dp) :: u
  ! >>>>>>>>>>>>>>>>>>>>>>
  u = dot_product(v, w)
end function Inner$0
#$do rank = 1, NUM_RANKS-1
pure function Inner$rank(v,w) result(u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: v(:), w(:,@:)
  real(dp) :: u(@{size(w, dim=$$)}@)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: i
  u(@:) = 0.0_dp
  do i = 1, size(v)
    u(@:) += v(i)*w(i,@:)
  end do
end function Inner$rank
#$end do
!! -----------------------------------------------------------------

!! -----------------------------------------------------------------  
pure function Outer$0(v,w) result(u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: v(:),w
  real(dp) :: u(size(v))
  ! >>>>>>>>>>>>>>>>>>>>>>
  u(:) = v(:)*w
end function Outer$0
#$do rank = 1, NUM_RANKS-1
pure function Outer$rank(v,w) result(u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: v(:), w(@:)
  real(dp) :: u(size(v,dim=1),@{size(w,dim=$$)}@)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: i
  do i = 1, size(v,dim=1)
    u(i,@:) = v(i)*w(@:)
  end do
end function Outer$rank
#$end do
!! -----------------------------------------------------------------  

end module StormRuler_Helpers
