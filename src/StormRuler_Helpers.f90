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
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module StormRuler_Helpers

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: ieee_arithmetic

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none  

abstract interface
  pure function tMapFunc$0(u) result(v)
    import dp
    real(dp), intent(in) :: u
    real(dp) :: v
  end function tMapFunc$0
#$do rank = 1, NUM_RANKS
  pure function tMapFunc$rank(u) result(v)
    import dp
    real(dp), intent(in) :: u(@:)
    real(dp) :: v(@{size(u, dim=$$)}@)
  end function tMapFunc$rank
#$end do
end interface

abstract interface
  pure function tSMapFunc$0(x, u) result(v)
    import dp
    real(dp), intent(in) :: x(:), u
    real(dp) :: v
  end function tSMapFunc$0
#$do rank = 1, NUM_RANKS
  pure function tSMapFunc$rank(x, u) result(v)
    import dp
    real(dp), intent(in) :: x(:), u(@:)
    real(dp) :: v(@{size(u, dim=$$)}@)
  end function tSMapFunc$rank
#$end do
end interface

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

interface
  logical function SameAddresses(a, b) bind(c, name='SameAddresses')
    type(*), intent(in) :: a(*), b(*)
  end function SameAddresses
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
!! ----------------------------------------------------------------- !!
subroutine BubbleSort(a)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(inout) :: a(:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: temp
  integer(ip) :: i, j
  logical :: swapped
 
  DO j = SIZE(a)-1, 1, -1
    swapped = .FALSE.
    DO i = 1, j
      IF (a(i) > a(i+1)) THEN
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END subroutine BubbleSort

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

!! ----------------------------------------------------------------- !!
!! Inner inner vector-tensor product: û ← v̅⋅ŵ.
!! ----------------------------------------------------------------- !!
pure function Inner$0(vBar, wHat) result(uHat)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: vBar(:), wHat(:)
  real(dp) :: uHat
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  uHat = dot_product(vBar, wHat)
end function Inner$0
#$do rank = 1, NUM_RANKS-1
pure function Inner$rank(vBar, wHat) result(uHat)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: vBar(:), wHat(:,@:)
  real(dp) :: uHat(@{size(wHat, dim=$$)}@)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: i

  uHat(@:) = 0.0_dp
  do i = 1, size(vBar)
    uHat(@:) = uHat(@:) + vBar(i)*wHat(i, @:)
  end do
end function Inner$rank
#$end do

!! ----------------------------------------------------------------- !!
!! Outer vector-tensor product: û ← v̅⊗ŵ.
!! ----------------------------------------------------------------- !!
pure function Outer$0(vBar, wHat) result(uHat)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: vBar(:), wHat
  real(dp) :: uHat(size(vBar))
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  uHat(:) = vBar(:)*wHat
end function Outer$0
#$do rank = 1, NUM_RANKS-1
pure function Outer$rank(vBar, wHat) result(uHat)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: vBar(:), wHat(@:)
  real(dp) :: uHat(size(vBar, dim=1), @{size(wHat, dim=$$)}@)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: i
  
  do i = 1, size(vBar, dim=1)
    uHat(i, @:) = vBar(i)*wHat(@:)
  end do
end function Outer$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Convert an integer(ip) to string.
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
!! Convert RGB pixel to integer(ip).
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

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

elemental function INTR_2(u_l, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l, u_r
  real(dp) :: INTR_2
  ! >>>>>>>>>>>>>>>>>>>>>>
  INTR_2 = 0.5_dp*(u_l + u_r)
end function INTR_2

elemental function EXTR_2(u_l, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l, u_r
  real(dp) :: EXTR_2
  ! >>>>>>>>>>>>>>>>>>>>>>
  EXTR_2 = 2.0_dp*u_l - u_r
end function EXTR_2

!! ----------------------------------------------------------------- !!
!! Interpolate cell-centered 2D values to node-centered values.
!! ----------------------------------------------------------------- !!
#$do rank = 0, NUM_RANKS
subroutine Interpolate2D_CellToNode$rank(u_2D_nc, u_2D_cc)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_2D_cc(@:,:,:)
  real(dp), intent(out) :: u_2D_nc(@:,:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: i, j, n, m

  n = size(u_2D_cc, dim=$rank+1)
  m = size(u_2D_cc, dim=$rank+2)

  ! ----------------------
  ! Corner interpolation.
  ! ----------------------
  !
  !       |       |       |
  !       o-------o-------o--
  !       |       |  2,2  |
  !       o-------o-------o--
  !       |  1,1  |       |
  ! (1,1) o-------o-------o--
  !
  u_2D_nc(@:,1,1) = EXTR_2(u_2D_cc(@:,1,1), u_2D_cc(@:,2,2))
  ! 
  !   |       |       |
  ! --o-------o-------o
  !   | n-1,2 |       |
  ! --o-------o-------o
  !   |       |  n,1  |
  ! --o-------o-------o (n+1,1)
  ! 
  u_2D_nc(@:,n+1,1) = EXTR_2(u_2D_cc(@:,n,1), u_2D_cc(@:,n-1,2))
  !
  ! --o-------o-------o (n+1,m+1)
  !   |       |  n,m  |
  ! --o-------o-------o
  !   |n-1,m-1|       |
  ! --o-------o-------o
  !   |       |       |
  !
  u_2D_nc(@:,n+1,m+1) = EXTR_2(u_2D_cc(@:,n,m), u_2D_cc(@:,n-1,m-1))
  !
  ! (1,m+1) o-------o-------o--
  !         |  1,m  |       |
  !         o-------o-------o--
  !         |       | 2,m-1 |
  !         o-------o-------o--
  !         |       |       |
  !
  u_2D_nc(@:,1,m+1) = EXTR_2(u_2D_cc(@:,1,m), u_2D_cc(@:,2,m-1))

  ! ----------------------
  ! Border interpolation.
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(n, m, u_2D_nc, u_2D_cc)
  do j = 2, m
    !
    !   m   |       |       |
    !   ^   o-------o-------o--
    !   |   |  1,j  |  2,j  |
    ! (1,j) o-------o-------o--
    !   |   | 1,j-1 | 2,j-1 |
    !   v   o-------o-------o--
    !   2   |       |       |
    !
    u_2D_nc(@:,1,j) = &
      & EXTR_2( INTR_2(u_2D_cc(@:,1,j), u_2D_cc(@:,1,j-1) ), &
      &         INTR_2(u_2D_cc(@:,2,j), u_2D_cc(@:,2,j-1) ) )
    ! 
    !   |       |       |     m
    ! --o-------o-------o     ^
    !   | n-1,j |  n,j  |     |
    ! --o-------o-------o (n+1,j)
    !   |n-1,j-1| n,j-1 |     |
    ! --o-------o-------o     v
    !   |       |       |     2
    !
    u_2D_nc(@:,n+1,j) = &
      & EXTR_2( INTR_2(u_2D_cc(@:,n  ,j), u_2D_cc(@:,n  ,j-1)), &
      &         INTR_2(u_2D_cc(@:,n-1,j), u_2D_cc(@:,n-1,j-1)) )
  end do
  !$omp end parallel do
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(n, m, u_2D_nc, u_2D_cc)
  do i = 2, n
    !
    !   |       |       |
    ! --o-------o-------o--
    !   | i-1,2 |  i,2  |
    ! --o-------o-------o--
    !   | i-1,1 |  i,1  |
    ! --o-------o-------o--
    !     2<--(i,1)-->n
    !
    u_2D_nc(@:,i,1) = &
      & EXTR_2( INTR_2(u_2D_cc(@:,i-1,1), u_2D_cc(@:,i,1)), &
      &         INTR_2(u_2D_cc(@:,i-1,2), u_2D_cc(@:,i,2)) )
    !
    !    2<--(i,m+1)-->n
    ! --o-------o-------o--
    !   | i-1,m |  i,m  |
    ! --o-------o-------o--
    !   |i-1,m-1| i,m-1 |
    ! --o-------o-------o--
    !   |       |       |
    !
    u_2D_nc(@:,i,m+1) = &
      & EXTR_2( INTR_2(u_2D_cc(@:,i-1,m  ), u_2D_cc(@:,i,m  )), &
      &         INTR_2(u_2D_cc(@:,i-1,m-1), u_2D_cc(@:,i,m-1)) )
  end do
  !$omp end parallel do

  ! ----------------------
  ! Interior interpolation.
  ! ----------------------
  !
  !   |       |       |     m
  ! --o-------o-------o--   ^
  !   | i-1,j |  i,j  |     |
  ! --o-------o-------o-- (i,j)
  !   |i-1,j-1| i,j-1 |     |
  ! --o-------o-------o--   v
  !   |       |       |     2
  !     2<--(i,j)-->n
  !
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(n, m, u_2D_nc, u_2D_cc)
  do j = 2, m
    do i = 2, n
      u_2D_nc(@:,i,j) = &
        & INTR_2( INTR_2(u_2D_cc(@:,i-1,j-1), u_2D_cc(@:,i-1,j)), &
        &         INTR_2(u_2D_cc(@:,i  ,j-1), u_2D_cc(@:,i  ,j)) )
    end do
  end do
  !$omp end parallel do
end subroutine Interpolate2D_CellToNode$rank
#$end do

end module StormRuler_Helpers
