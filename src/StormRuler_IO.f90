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

use StormRuler_Helpers
#$use 'StormRuler_Parameters.f90'

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

type, abstract :: IOListItem
  character(len=:), allocatable :: name
  class(IOListItem), pointer :: next => null()
end type IOListItem

#$do rank = 0, NUM_RANKS
type, extends(IOListItem) :: IOListItem$rank
  real(dp), pointer :: values(@:,:) => null()
end type !IOListItem$rank
#$end do

type :: IOList
  class(IOListItem), pointer :: first => null()
contains
  generic :: Add=>Add0,Add1,Add2
#$do rank = 0, NUM_RANKS
  procedure :: Add$rank=>IOList_Add$rank
#$end do
end type IOList

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

#$do rank = 0, NUM_RANKS
subroutine IOList_Add$rank(list,name,values)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(IOList), intent(inout) :: list
  character(len=*), intent(in) :: name
  real(dp), target, intent(in) :: values(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
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

pure function IntToPixel(int) result(colorChannels)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer, intent(in) :: int
  integer :: colorChannels(3)
  ! >>>>>>>>>>>>>>>>>>>>>>
  colorChannels(1) = iand(255,int)
  colorChannels(2) = iand(255,ishft(int,-8))
  colorChannels(3) = iand(255,ishft(int,-16))
end function IntToPixel

!! -----------------------------------------------------------------
!! Load PPM image.
subroutine Load_PPM(file,pixels)
  ! <<<<<<<<<<<<<<<<<<<<<<
  character(len=*), intent(in) :: file
  integer, allocatable, intent(out) :: pixels(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: unit,offset
  integer :: numRows,numColumns
  ! ----------------------
  ! Parse PPM header.
  ! ----------------------
  block
    character(len=2) :: magic
    integer :: colorRange
    open(newunit=unit,file=file, &
         access='stream',form='formatted',status='old')
    read(unit,'(A2)') magic
    if (magic/='P6') &
      error stop 'unexpected PPM magic, "P6" expected'
    read(unit,*) numRows,numColumns
    read(unit,*) colorRange
    if (colorRange/=255) &
      error stop 'unsupported PPM color range value, "255" expected'
    inquire(unit,pos=offset)
    close(unit)
  end block
  ! ----------------------
  ! Allocate and read image pixels.
  ! ----------------------
  block
    character :: byte,bytes(3)
    integer :: row,column,colorChannel
    allocate(pixels(0:numRows-1,0:numColumns-1))
    open(newunit=unit,file=file, &
         access='stream',status='old')
    read(unit,pos=offset-1) byte
    do column = numColumns-1, 0, -1
      do row = 0, numRows-1
        do colorChannel = 1, 3
          read(unit) byte; bytes(colorChannel) = byte
        end do
        pixels(row,column) = PixelToInt(iachar(bytes))
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
  ! ----------------------
  numRows = size(pixels,dim=1) 
  numColumns = size(pixels,dim=2)
  open(newunit=unit,file=file,status='replace')
  write(unit,'(A2)') 'P6'
  write(unit,'(I0," ",I0)') numRows,numColumns
  write(unit,'(I0)') 255
  ! ----------------------
  ! Write PPM image pixels.
  ! ----------------------
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

end module StormRuler_IO