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

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Print an epic banner.
!! ----------------------------------------------------------------- !!
subroutine PrintBanner
  
  print *, ''
  print *, '//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\\'
  print *, '|      _____ __                       ____        __            |'
  print *, '|     / ___// /_____  _________ ___  / __ \__  __/ ___  _____   |'
  print *, '|     \__ \/ __/ __ \/ ___/ __ `__ \/ /_/ / / / / / _ \/ ___/   |'
  print *, '|    ___/ / /_/ /_/ / /  / / / / / / _, _/ /_/ / /  __/ /       |'
  print *, '|   /____/\__/\____/_/  /_/ /_/ /_/_/ |_|\__,_/_/\___/_/        |'
  print *, '|                                                               |'
  print *, '\\-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//'
  print *, ''

end subroutine PrintBanner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Swap function.
!! ----------------------------------------------------------------- !!
recursive subroutine Swap(i, j)
  integer(ip), intent(inout) :: i, j

  integer(ip) :: k

  k = i; i = j; j = k

end subroutine Swap

!! ----------------------------------------------------------------- !!
!! Flip integer oddity, e.g. 1→2, 2→1, 3→4, 4→3, ...
!! ----------------------------------------------------------------- !!
integer(ip) pure function Flip(value)
  integer(ip), intent(in) :: value

  Flip = merge(value+1, value-1, mod(value, 2) == 1)

end function Flip

!! ----------------------------------------------------------------- !!
!! Find value index in the array.
!! ----------------------------------------------------------------- !!
integer(ip) pure function IndexOf(value, array)
  integer(ip), intent(in) :: value, array(:)
  
  do IndexOf = 1, size(array)
    if (array(IndexOf) == value) return
  end do
  IndexOf = 0

end function IndexOf

!! ----------------------------------------------------------------- !!
!! Find value index in the sorted array.
!! ----------------------------------------------------------------- !!
integer(ip) pure function SortedIndexOf(value, array) result(pivot)
  integer(ip), intent(in) :: value, array(:)
  
  integer(ip) :: first, last

  first = 1
  if (value < array(first)) then
    pivot = 0; return  
  end if

  last = size(array)
  if (value > array(last)) then
    pivot = 0; return  
  end if
  
  do while (first + 1 < last)
    pivot = (first + last)/2
    if (array(pivot) == value) then
      return
    else if (array(pivot) > value) then
      last = pivot
    else if (array(pivot) < value) then
      first = pivot
    end if
  end do
  pivot = 0

end function SortedIndexOf

!! ----------------------------------------------------------------- !!
!! Check if two sorted arrays inetersect.
!! ----------------------------------------------------------------- !!
logical pure function ArraysIntersect(left, right) result(intersects)
  integer(ip), intent(in) :: left(:), right(:)

  integer(ip) :: i, j

  i = 1; j = 1; intersects = .false.
  do while((i <= size(left)).and.(j <= size(right)))
    if (left(i) == right(j)) then
      intersects = .true.; return
    else if (left(i) < right(j)) then
      i = i + 1
    else if (left(i) > right(j)) then
      j = j + 1
    end if
  end do

end function ArraysIntersect

!! ----------------------------------------------------------------- !!
!! Bubble-sort the integer array.
!! ----------------------------------------------------------------- !!
subroutine BubbleSort(array)
  integer(ip), intent(inout) :: array(:)
  
  logical :: swapped
  integer(ip) :: i, j
 
  do j = size(array) - 1, 1, -1
    swapped = .false.
    do i = 1, j
      if (array(i) > array(i+1)) then
        call Swap(array(i), array(i+1))
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
  real(dp), intent(in) :: value
  
  if (ieee_is_nan(value).or.(value <= 0)) then
    error stop 'NEGATIVE, ZERO OR NaN VALUE: '//R2S(value)
  end if

end subroutine EnsurePositive

!! ----------------------------------------------------------------- !!
!! Ensure the value is positive or zero.
!! ----------------------------------------------------------------- !!
subroutine EnsureNonNegative(value)
  real(dp), intent(in) :: value
  
  if (ieee_is_nan(value).or.(value < 0)) then
    error stop 'NEGATIVE OR NaN VALUE: '//R2S(value)
  end if

end subroutine EnsureNonNegative

!! ----------------------------------------------------------------- !!
!! Compute pseudo inverse: 𝑎⁺ ← 1/𝑎 𝗶𝗳 𝑎 ≠ 0 𝗲𝗹𝘀𝗲 0.
!! ----------------------------------------------------------------- !!
pure real(dp) elemental function SafeInverse(a)
  real(dp), intent(in) :: a
  
  SafeInverse = merge(0.0_dp, 1.0_dp/a, abs(a) <= tiny(1.0_dp))

end function SafeInverse

!! ----------------------------------------------------------------- !!
!! Divide with pseudo inverse: 𝑑 ← 𝑏⁺⋅𝑎.
!! ----------------------------------------------------------------- !!
real(dp) function SafeDivide(a, b)
  real(dp), intent(in) :: a, b
  
  SafeDivide = SafeInverse(b)*a

end function SafeDivide

!! ----------------------------------------------------------------- !!
!! Solve a small dense linear SLAE: 𝓐𝒙 = 𝒃.
!! ----------------------------------------------------------------- !!
pure function DenseSolve(A, b) result(x)
  real(dp), intent(in) :: A(:,:), b(:)
  real(dp) :: x(size(b))

  !! TODO: real implementation is missing!
  x(:) = b(:)/A(1,1)

end function DenseSolve

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Convert an integer to string.
!! ----------------------------------------------------------------- !!
pure function I2S(value)
  integer(ip), intent(in) :: value
  character(len=:), allocatable :: I2S
  
  character(len=256) :: buffer
  
  write(buffer, *) value
  I2S = trim(adjustl(buffer))

end function I2S

!! ----------------------------------------------------------------- !!
!! Convert a real number to string.
!! ----------------------------------------------------------------- !!
pure function R2S(value)
  real(dp), intent(in) :: value
  character(len=:), allocatable :: R2S
  
  character(len=256) :: buffer

  write(buffer, *) value
  R2S = trim(adjustl(buffer))

end function R2S

!! ----------------------------------------------------------------- !!
!! Ternary operator for strings.
!! ----------------------------------------------------------------- !!
pure function MergeString(trueString, falseString, condition)
  character(len=*), intent(in) :: trueString, falseString
  logical, intent(in) :: condition
  character(len=:), allocatable :: MergeString
  
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
  integer(ip), intent(in) :: colorChannels(3)
  integer(ip) :: int
  
  int = ior(iand(255, colorChannels(1)), &
    &       ior(ishft(iand(255, colorChannels(2)), 8), &
    &           ishft(iand(255, colorChannels(3)), 16)))

end function PixelToInt

!! ----------------------------------------------------------------- !!
!! Convert integer(ip) to RGB pixel.
!! ----------------------------------------------------------------- !!
pure function IntToPixel(int) result(colorChannels)
  integer(ip), intent(in) :: int
  integer(ip) :: colorChannels(3)

  colorChannels(1) = iand(255, int)
  colorChannels(2) = iand(255, ishft(int, -8))
  colorChannels(3) = iand(255, ishft(int, -16))

end function IntToPixel

end module StormRuler_Helpers
