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
#@use 'StormRuler_Parameters.f90'

implicit none  

interface Inner
#@do rank = 0, NumRanks-1
  module procedure Inner$rank
#@end do
end interface Inner

interface Outer
#@do rank = 0, NumRanks-1
  module procedure Outer$rank
#@end do
end interface Outer

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
!! Convert an integer to string.
function IntToString(value) result(string)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer, intent(in) :: value
  character(len=20) :: string
  ! >>>>>>>>>>>>>>>>>>>>>>
  write(string, *) value
  string = adjustl(string)
end function IntToString

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
#@do rank = 1, NumRanks-1
pure function Inner$rank(v,w) result(u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: v(:), w(:,@:)
  real(dp) :: u(@{size(w, dim=$$+1)}@)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: i
  u(@:) = 0.0_dp
  do i = 1, size(v)
    u(@:) += v(i)*w(i,@:)
  end do
end function Inner$rank
#@end do
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
pure function Outer$0(v,w) result(u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: v(:),w
  real(dp) :: u(size(v))
  ! >>>>>>>>>>>>>>>>>>>>>>
  u = v*w
end function Outer$0
#@do rank = 1, NumRanks-1
pure function Outer$rank(v,w) result(u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: v(:), w(@:)
  real(dp) :: u(size(v), @{size(w, dim=$$+1)}@)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: i
  do i = 1, size(v)
    u(i,@:) = v(i)*w(@:)
  end do
end function Outer$rank
#@end do
!! -----------------------------------------------------------------  

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
pure function PixelToInt(colorChannels) result(int)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer, intent(in) :: colorChannels(3)
  integer :: int
  ! >>>>>>>>>>>>>>>>>>>>>>
  int = ior(iand(255,colorChannels(1)), &
            ior(ishft(iand(255,colorChannels(2)),8), &
                ishft(iand(255,colorChannels(3)),16)))
end function PixelToInt
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------
!! Load PPM image.
subroutine Load_PPM(file,pixels)
  ! <<<<<<<<<<<<<<<<<<<<<<
  character(len=*), intent(in) :: file
  integer, intent(out), allocatable :: pixels(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: unit,offset
  integer :: numRows,numColumns
  ! ----------------------
  ! Parse PPM header.
  block
    character(len=2) :: magic
    integer :: colorRange
    open(newunit=unit,file=file, &
         access='stream',form='formatted',status='old')
    read(unit,'(A2)') magic
    if (magic/='P6') &
      error stop 'unexpected PPM magic, P6 expected'
    read(unit,*) numRows,numColumns
    read(unit,*) colorRange
    if (colorRange/=255) &
      error stop 'unsupported PPM color range, 255 expected'
    inquire(unit, pos=offset)
    close(unit)
  end block
  ! ----------------------
  ! Allocate and read image pixels.
  block
    character :: byte
    integer :: row,column,colorChannel
    allocate(pixels(numRows,numColumns))
    open(newunit=unit,file=file, &
         access='stream',status='old')
    read(unit, pos=offset-1) byte
    do column = numColumns, 1, -1
      do row = 1, numRows
        do colorChannel = 0, 2
          read(unit) byte
          pixels(row,column) = &
            ior(pixels(row,column), ishft(iachar(byte),8*colorChannel))
        end do
      end do
    end do
    close(unit)
  end block
end subroutine Load_PPM
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------
!! Save PPM image.  
subroutine Save_PPM(file,pixels)
  ! <<<<<<<<<<<<<<<<<<<<<<
  character(len=*), intent(in) :: file
  integer, intent(in) :: pixels(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: unit
  integer :: numRows,numColumns
  integer :: row,column,colorChannel
  ! ----------------------
  ! Write PPM header.
  numRows = size(pixels,dim=1) 
  numColumns = size(pixels,dim=2)
  open(newunit=unit,file=file,status='replace')
  write(unit,'(A2)') 'P6'
  write(unit,'(I0," ",I0)') numRows,numColumns
  write(unit,'(I0)') 255
  ! ----------------------
  ! Write PPM image pixels.
  do column = numColumns, 1, -1
    do row = 1, numRows
      do colorChannel = 0, 2
        write(unit,'(A1)',advance='no') &
          achar(iand(255,ishft(pixels(row,column),-8*colorChannel)))
      end do
    end do
  end do
  close(unit)
end subroutine Save_PPM
!! -----------------------------------------------------------------  

end module StormRuler_Helpers
