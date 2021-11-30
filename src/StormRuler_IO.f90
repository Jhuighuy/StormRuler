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
module StormRuler_IO

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Helpers, only: I2S, RgbToInt, IntToRgb

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

type, abstract :: IOListItem
  character(len=:), allocatable :: name
  class(IOListItem), pointer :: next => null()
end type IOListItem

#$do rank = 0, 2
type, extends(IOListItem) :: IOListItem$rank
  real(dp), pointer :: values(@:,:) => null()
end type !IOListItem$rank
#$end do

type :: IOList
  class(IOListItem), pointer :: first => null()
contains
  generic :: Add => @{Add$$@|@0, 2}@
#$do rank = 0, 2
  procedure :: Add$rank => IOList_Add$rank
#$end do
end type IOList

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

#$do rank = 0, 2
subroutine IOList_Add$rank(list, name, values)
  class(IOList), intent(inout) :: list
  character(len=*), intent(in) :: name
  real(dp), target, intent(in) :: values(@:,:)

  class(IOListItem$rank), pointer :: item

  ! ----------------------
  ! Generate item.
  ! ----------------------
  allocate(item)
  item%name = name
  item%values => values
  ! ----------------------
  ! Append item.
  ! ----------------------
  item%next => list%first
  list%first => item
end subroutine IOList_Add$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Load PPM image.
!! ----------------------------------------------------------------- !!
subroutine Load_PPM(file, image)
  character(len=*), intent(in) :: file
  integer(ip), allocatable, intent(out) :: image(:,:)

  integer(ip) :: unit, offset
  integer(ip) :: numRows, numCols

  ! ----------------------
  ! Parse PPM header.
  ! ----------------------
  block
    character(len=2) :: magic
    integer(ip) :: colorRange
    open(newunit=unit, file=file, access='stream', form='formatted', status='old')
    read(unit, '(A2)') magic
    if (magic/='P6') call ErrorStop('unexpected PPM magic, `P6` expected')
    read(unit, *) numRows, numCols
    read(unit, *) colorRange
    if (colorRange/=255) call ErrorStop('unsupported PPM color range value, "255" expected')
    inquire(unit, pos=offset)
    close(unit)
  end block
  
  ! ----------------------
  ! Allocate and read image image.
  ! ----------------------
  block
    character :: byte, bytes(3)
    integer(ip) :: row, col, rgb
    allocate(image(0:(numRows-1),0:(numCols-1)))
    open(newunit=unit, file=file, access='stream', status='old')
    read(unit, pos=offset-1) byte
    do col = numCols - 1, 0, -1
      do row = 0, numRows - 1
        do rgb = 1, 3
          read(unit) byte; bytes(rgb) = byte
        end do
        image(row,col) = RgbToInt(iachar(bytes))
      end do
    end do
    close(unit)
  end block
end subroutine Load_PPM

!! ----------------------------------------------------------------- !!
!! Save PPM image.  
!! ----------------------------------------------------------------- !!
subroutine Save_PPM(file, image)
  character(len=*), intent(in) :: file
  integer(ip), intent(in) :: image(:,:)

  integer(ip) :: unit
  character :: bytes(3)
  integer(ip) :: row, numRows, col, numCols
  
  open(newunit=unit, file=file, access='stream', status='replace')

  ! ----------------------
  ! Write PPM header.
  ! ----------------------
  numRows = size(image, dim=1) 
  numCols = size(image, dim=2)
  write(unit) 'P6'//char(10)
  write(unit) I2S(numRows)//' '//I2S(numCols)//char(10)
  write(unit) '255'//char(10)
  
  ! ----------------------
  ! Write PPM image image.
  ! ----------------------
  do col = 1, numCols
    do row = 1, numRows
      bytes(:) = achar(IntToRgb(image(row,col)))
      write(unit) bytes(1)
      write(unit) bytes(2)
      write(unit) bytes(3)
    end do
  end do

  close(unit)

end subroutine Save_PPM

end module StormRuler_IO