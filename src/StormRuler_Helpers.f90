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
module StormRuler_Helpers

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
use, intrinsic :: iso_fortran_env, only: error_unit

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

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

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Flip integer oddity, e.g. 1→2, 2→1, 3→4, 4→3, ...
!! ----------------------------------------------------------------- !!
integer(ip) pure function Flip(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: value
  ! >>>>>>>>>>>>>>>>>>>>>>

  Flip = merge(value+1, value-1, mod(value, 2) == 1)
end function Flip

!! ----------------------------------------------------------------- !!
!! Find value index in the array.
!! ----------------------------------------------------------------- !!
integer(ip) pure function IndexOf(value, array)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: value, array(:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  do IndexOf = 1, size(array)
    if (array(IndexOf) == value) return
  end do
  IndexOf = 0
end function IndexOf

!! ----------------------------------------------------------------- !!
!! Bubble-sort the integer array.
!! ----------------------------------------------------------------- !!
subroutine BubbleSort(array)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(inout) :: array(:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: temp
  integer(ip) :: i, j
  logical :: swapped
 
  do j = size(array)-1, 1, -1
    swapped = .false.
    do i = 1, j
      if (array(i) > array(i+1)) then
        temp = array(i)
        array(i) = array(i+1)
        array(i+1) = temp
        swapped = .true.
      end if
    end do
    if (.not.swapped) exit
  end do
end subroutine BubbleSort

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Ensure the value is positive.
!! ----------------------------------------------------------------- !!
subroutine EnsurePositive(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: value
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  if (ieee_is_nan(value).or.(value <= 0)) then
    write(error_unit, *) 'NEGATIVE, ZERO OR NaN VALUE', value
    error stop 1
  end if
end subroutine EnsurePositive

!! ----------------------------------------------------------------- !!
!! Ensure the value is positive or zero.
!! ----------------------------------------------------------------- !!
subroutine EnsureNonNegative(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: value
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  if (ieee_is_nan(value).or.(value < 0)) then
    write(error_unit, *) 'NEGATIVE OR NaN VALUE', value
    error stop 1
  end if
end subroutine EnsureNonNegative

!! ----------------------------------------------------------------- !!
!! ----------------------------------------------------------------- !!
real(dp) function SafeDivide(a, b)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: a, b
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  call EnsurePositive(abs(b))
  SafeDivide = a/b
end function SafeDivide

!! ----------------------------------------------------------------- !!
!! Compute pseudo inverse: i ← 1/a if a ≠ 0 else 0.
!! ----------------------------------------------------------------- !!
real(dp) elemental function SafeInverse(a)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: a
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  SafeInverse = merge(0.0_dp, 1.0_dp/a, a==0.0_dp)
end function SafeInverse

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

impure elemental subroutine Assert(b)

  logical, intent(in) :: b

  if (.not.b) then
    error stop 'assertion failed'
  end if

end subroutine Assert

function Reshape2D(u) result(v)
  use, intrinsic :: iso_c_binding

  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in), target :: u(..)
  real(dp), pointer :: v(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  select rank(u)
    rank(1)
      call c_f_pointer(cptr=c_loc(u), fptr=v, &
        & shape=[ 1, size(u) ])
#$do rank = 1, NUM_RANKS
    rank($rank + 1)
      call c_f_pointer(cptr=c_loc(u), fptr=v, &
        & shape=[ size(u(@:,1)), size(u(@1,:)) ])
#$end do
    rank default
      error stop 'Invalid rank'
  end select

end function Reshape2D

function Reshape3D(dim, u) result(v)
  use, intrinsic :: iso_c_binding

  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: dim
  real(dp), intent(in), target :: u(..)
  real(dp), pointer :: v(:,:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  select rank(u)
#$do rank = 1, NUM_RANKS
    rank($rank + 1)
      call c_f_pointer(cptr=c_loc(u), fptr=v, &
        & shape=[ dim, size(u(@:,1))/dim, size(u(@1,:)) ])
#$end do
    rank default
      error stop 'Invalid rank'
  end select

end function Reshape3D

function Reshape4D(dim, u) result(v)
  use, intrinsic :: iso_c_binding

  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: dim
  real(dp), intent(in), target :: u(..)
  real(dp), pointer :: v(:,:,:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  select rank(u)
#$do rank = 1, NUM_RANKS
    rank($rank + 1)
      call c_f_pointer(cptr=c_loc(u), fptr=v, &
        & shape=[ dim, dim, size(u(@:,1))/dim**2, size(u(@1,:)) ])
#$end do
    rank default
      error stop 'Invalid rank'
  end select

end function Reshape4D

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Inner vector-vector product: z ← y⋅x.
!! ----------------------------------------------------------------- !!
pure function Inner$0(y, x) result(z)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: y(:), x(:)
  real(dp) :: z
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  z = dot_product(y, x)

end function Inner$0
#$do rank = 1, NUM_RANKS-1
pure function Inner$rank(y, x) result(z)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: y(:), x(:,@:)
  real(dp) :: z(@{size(x, dim=$$+1)}@)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: i

  z(@:) = 0.0_dp
  do i = 1, size(y)
    z(@:) = z(@:) + y(i)*x(i,@:)
  end do

end function Inner$rank
#$end do

!! ----------------------------------------------------------------- !!
!! Outer vector-tensor product: z ← y⊗x.
!! ----------------------------------------------------------------- !!
pure function Outer$0(y, x) result(z)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: y(:), x
  real(dp) :: z(size(y))
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  z(:) = y(:)*x

end function Outer$0
#$do rank = 1, NUM_RANKS-1
pure function Outer$rank(y, x) result(z)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: y(:), x(@:)
  real(dp) :: z(size(y, dim=1), @{size(x, dim=$$)}@)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: i
  
  do i = 1, size(y, dim=1)
    z(i,@:) = y(i)*x(@:)
  end do

end function Outer$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Get real part of a complex number.
!! ----------------------------------------------------------------- !!
real(dp) elemental function Re(z)
  ! <<<<<<<<<<<<<<<<<<<<<<
  complex(dp), intent(in) :: z
  ! >>>>>>>>>>>>>>>>>>>>>>

  Re = real(z, kind=dp)
end function Re

!! ----------------------------------------------------------------- !!
!! Get imaginary part of a complex number.
!! ----------------------------------------------------------------- !!
real(dp) elemental function Im(z)
  ! <<<<<<<<<<<<<<<<<<<<<<
  complex(dp), intent(in) :: z
  ! >>>>>>>>>>>>>>>>>>>>>>

  Im = aimag(z) !! TODO: is this operation double-precision preserving?
end function Im

!! ----------------------------------------------------------------- !!
!! Convert real number (or a pair of two) into the complex number.
!! ----------------------------------------------------------------- !!
complex(dp) elemental function R2C(x,y)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: x
  real(dp), intent(in), optional :: y
  ! >>>>>>>>>>>>>>>>>>>>>>

  if (present(y)) then
    R2C = cmplx(x, y, kind=dp)
  else
    R2C = cmplx(x, 0.0_dp, kind=dp)
  end if
end function R2C

!! ----------------------------------------------------------------- !!
!! Convert an integer to string.
!! ----------------------------------------------------------------- !!
function I2S(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: value
  character(len=:), allocatable :: I2S
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  character(len=256) :: buffer
  write(buffer, *) value
  I2S = trim(adjustl(buffer))
end function I2S

!! ----------------------------------------------------------------- !!
!! Convert a real number to string.
!! ----------------------------------------------------------------- !!
function R2S(value)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: value
  character(len=:), allocatable :: R2S
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  character(len=256) :: buffer
  write(buffer, *) value
  R2S = trim(adjustl(buffer))
end function R2S

!! ----------------------------------------------------------------- !!
!! Ternary operator for strings.
!! ----------------------------------------------------------------- !!
function MergeString(trueString, falseString, condition)
  ! <<<<<<<<<<<<<<<<<<<<<<
  character(len=*), intent(in) :: trueString, falseString
  logical, intent(in) :: condition
  character(len=:), allocatable :: MergeString
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  if (condition) then
    MergeString = trueString
  else
    MergeString = falseString
  end if
end function MergeString

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Convert RGB pixel to integer.
!! ----------------------------------------------------------------- !!
pure function PixelToInt(colorChannels) result(int)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: colorChannels(3)
  integer(ip) :: int
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  int = ior(iand(255, colorChannels(1)), &
    &       ior(ishft(iand(255, colorChannels(2)), 8), &
    &           ishft(iand(255, colorChannels(3)), 16)))
end function PixelToInt

!! ----------------------------------------------------------------- !!
!! Convert integer(ip) to RGB pixel.
!! ----------------------------------------------------------------- !!
pure function IntToPixel(int) result(colorChannels)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: int
  integer(ip) :: colorChannels(3)
  ! >>>>>>>>>>>>>>>>>>>>>>

  colorChannels(1) = iand(255, int)
  colorChannels(2) = iand(255, ishft(int, -8))
  colorChannels(3) = iand(255, ishft(int, -16))
end function IntToPixel

end module StormRuler_Helpers
