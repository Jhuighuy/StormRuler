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

use StormRuler_Consts, only: ip, dp

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface L2S
  module procedure L2S
  module procedure LArr2S
end interface L2S

interface I2S
  module procedure I2S
  module procedure IArr2S
end interface I2S

interface R2S
  module procedure R2S
  module procedure RArr2S
end interface R2S

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print an epic banner.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine PrintBanner
  
  print *, ''
  print *, '//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\\'
  print *, '|      _____ __                       ____        __            |'
  print *, '|     / ___// /_____   ____ ___ ___  / __ \__  __/ ___  _____   |'
  print *, '|     \__ \/ __/`__ \/ ___/`__ `__ \/ /_/ / / / / / _ \/ ___/   |'
  print *, '|    ___/ / /_/ /_/ / /  / / / / / / _, _/ /_/ / /  __/ /       |'
  print *, '|   /____/\__/\____/_/  /_/ /_/ /_/_/ |_|\__,_/_/\___/_/        |'
  print *, '|                                                               |'
  print *, '\\-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//'
  print *, ''

end subroutine PrintBanner

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Error stop with message.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ErrorStop(message)
  character(len=*), intent(in), optional :: message

  if (present(message)) write(error_unit,*) message
  error stop 1

end subroutine ErrorStop

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print a log message.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine PrintLog(message)
  character(len=*), intent(in) :: message

  print *, '[LOG] ', message

end subroutine PrintLog

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print a warning.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine PrintWarning(warning)
  character(len=*), intent(in) :: warning

  write(error_unit,*) '[WRN] ', warning

end subroutine PrintWarning

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

!! ----------------------------------------------------------------- !!
!! Insert value into the allocatable array.
!! ----------------------------------------------------------------- !!
subroutine InsertTo(array, value)
  integer(ip), intent(inout), allocatable :: array(:)
  integer(ip), intent(in) :: value

  if (allocated(array)) then
    array = [array, value]
  else
    array = [value]
  end if

end subroutine InsertTo

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Inverse and compress a mapping: 𝜑(𝑖) = 𝑗 into 𝜑⁻¹(𝑗) = {𝑖: 𝜑(𝑖) = 𝑗}. 
!! ----------------------------------------------------------------- !!
subroutine InverseCompressMapping(addrs, indices, num, mapping, width)
  integer(ip), intent(inout), allocatable :: addrs(:), indices(:)
  integer(ip), intent(in) :: num, mapping(:)
  integer(ip), intent(out), optional :: width

  integer(ip) :: i, j

  allocate(addrs(num + 1), indices(size(mapping)))

  addrs(1) = 1; addrs(2:) = 0
  do i = 1, size(mapping)
    j = mapping(i)
    addrs(j + 1) = addrs(j + 1) + 1
  end do
  if (present(width)) width = 1
  do j = 1, num
    if (present(width)) width = max(width, addrs(j + 1))
    addrs(j + 1) = addrs(j + 1) + addrs(j)
  end do

  do i = 1, size(mapping)
    j = mapping(i)
    indices(addrs(j)) = i
    addrs(j) = addrs(j) + 1
  end do
  addrs = cshift(addrs, shift=-1); addrs(1) = 1
  
end subroutine InverseCompressMapping

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Ensure the value is positive.
!! ----------------------------------------------------------------- !!
subroutine EnsurePositive(value)
  real(dp), intent(in) :: value
  
  if (ieee_is_nan(value).or.(value <= 0)) then
    call ErrorStop('unexpected negative, zero or NaN value: '//R2S(value))
  end if

end subroutine EnsurePositive

!! ----------------------------------------------------------------- !!
!! Ensure the value is positive or zero.
!! ----------------------------------------------------------------- !!
subroutine EnsureNonNegative(value)
  real(dp), intent(in) :: value
  
  if (ieee_is_nan(value).or.(value < 0)) then
    call ErrorStop('unexpected negative or NaN value: '//R2S(value))
  end if

end subroutine EnsureNonNegative

!! ----------------------------------------------------------------- !!
!! Compute pseudo-inverse: 𝑎⁺ ← 1/𝑎 𝗶𝗳 𝑎 ≠ 0 𝗲𝗹𝘀𝗲 0.
!! ----------------------------------------------------------------- !!
pure real(dp) elemental function SafeInverse(a)
  real(dp), intent(in) :: a
  
  SafeInverse = merge(0.0_dp, 1.0_dp/a, abs(a) <= tiny(1.0_dp))

end function SafeInverse

!! ----------------------------------------------------------------- !!
!! Divide with pseudo-inverse: 𝑑 ← 𝑏⁺⋅𝑎.
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

!! ----------------------------------------------------------------- !!
!! Convert a logical value to string.
!! ----------------------------------------------------------------- !!
pure function L2S(value)
  logical, intent(in) :: value
  character(len=:), allocatable :: L2S
  
  L2S = MergeString('true', 'false', value)

end function L2S
pure function LArr2S(array)
  logical, intent(in) :: array(:)
  character(len=:), allocatable :: LArr2S

  integer(ip) :: index

  LArr2S = L2S(array(1))
  do index = 2, size(array)
    LArr2S = LArr2S//' '//L2S(array(index))
  end do

end function LArr2S

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
pure function IArr2S(array)
  integer(ip), intent(in) :: array(:)
  character(len=:), allocatable :: IArr2S

  integer(ip) :: index

  IArr2S = I2S(array(1))
  do index = 2, size(array)
    IArr2S = IArr2S//' '//I2S(array(index))
  end do

end function IArr2S

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
pure function RArr2S(array)
  real(dp), intent(in) :: array(:)
  character(len=:), allocatable :: RArr2S

  integer(ip) :: index

  RArr2S = R2S(array(1))
  do index = 2, size(array)
    RArr2S = RArr2S//' '//R2S(array(index))
  end do

end function RArr2S

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Convert RGB pixel to integer.
!! ----------------------------------------------------------------- !!
pure function RgbToInt(rgb) result(int)
  integer(ip), intent(in) :: rgb(3)
  integer(ip) :: int
  
  int = ior(iand(255, rgb(1)), &
    &   ior(ishft(iand(255, rgb(2)), 8), &
    &       ishft(iand(255, rgb(3)), 16)))

end function RgbToInt

!! ----------------------------------------------------------------- !!
!! Convert integer(ip) to RGB pixel.
!! ----------------------------------------------------------------- !!
pure function IntToRgb(int) result(rgb)
  integer(ip), intent(in) :: int
  integer(ip) :: rgb(3)

  rgb(1) = iand(255, int)
  rgb(2) = iand(255, ishft(int, -8))
  rgb(3) = iand(255, ishft(int, -16))

end function IntToRgb

end module StormRuler_Helpers
