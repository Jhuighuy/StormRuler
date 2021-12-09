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
module StormRuler_IO_Writer

use StormRuler_Consts, only: bp, ip, lp, sp, dp

use StormRuler_Helpers, only: ErrorStop, I2S, R2S

use StormRuler_IO_Stream, only: tOutputStream

use, intrinsic :: iso_c_binding, only: c_size_t, c_long, &
  & c_loc, c_f_pointer, c_sizeof

#$use 'StormRuler_Macros.fi'

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Abstract data writer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, abstract :: tWriter
contains
  generic :: Write => &
    & WriteByte, WriteByteArray, &
    & WriteInteger, WriteIntegerArray, &
    & WriteReal, WriteRealVector, WriteRealVectorArray
  procedure(tWriteByteToStreamFunc), deferred :: WriteByte
  procedure(tWriteIntegerToStreamFunc), deferred :: WriteInteger
  procedure(tWriteRealToStreamFunc), deferred :: WriteReal
  procedure :: WriteByteArray => WriteByteArrayToStream
  procedure :: WriteIntegerArray => WriteIntegerArrayToStream
  procedure :: WriteRealVector => WriteRealVectorToStream
  procedure :: WriteRealVectorArray => WriteRealVectorArrayToStream

end type tWriter

!! ----------------------------------------------------------------- !!
!! Encode and write the data.
!! ----------------------------------------------------------------- !!
abstract interface
  subroutine tWriteByteToStreamFunc(writer, stream, data)
    import :: tWriter, tOutputStream, bp
    class(tWriter), intent(inout) :: writer
    class(tOutputStream), intent(inout) :: stream
    integer(bp), intent(in) :: data
  end subroutine tWriteByteToStreamFunc
  subroutine tWriteIntegerToStreamFunc(writer, stream, data)
    import :: tWriter, tOutputStream, ip
    class(tWriter), intent(inout) :: writer
    class(tOutputStream), intent(inout) :: stream
    integer(ip), intent(in) :: data
  end subroutine tWriteIntegerToStreamFunc
  subroutine tWriteRealToStreamFunc(writer, stream, data)
    import :: tWriter, tOutputStream, dp
    class(tWriter), intent(inout) :: writer
    class(tOutputStream), intent(inout) :: stream
    real(dp), intent(in) :: data
  end subroutine tWriteRealToStreamFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Text encoded writer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tWriter) :: tTextWriter
  character(len=:), allocatable :: Separator

contains
  procedure, non_overridable :: WriteByte => WriteTextByteToStream
  procedure, non_overridable :: WriteInteger => WriteTextIntegerToStream
  procedure, non_overridable :: WriteReal => WriteTextRealToStream

end type tTextWriter

interface tTextWriter
  module procedure MakeTextWriter
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Binary encoded writer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tWriter) :: tBinaryWriter
  logical, private :: SwapEndian
  logical, private :: LongInts, SingleReals

contains
  procedure, non_overridable :: WriteByte => WriteBinaryByteToStream
  procedure, non_overridable :: WriteInteger => WriteBinaryIntegerToStream
  procedure, non_overridable :: WriteReal => WriteBinaryRealToStream

end type tBinaryWriter

interface tBinaryWriter
  module procedure MakeBinaryWriter
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Encode arrays.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine WriteByteArrayToStream(writer, stream, data)
  class(tWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(bp), intent(in) :: data(:)

  integer(ip) :: index

  do index = 1, size(data)
    call writer%WriteByte(stream, data(index))
  end do

end subroutine WriteByteArrayToStream
subroutine WriteIntegerArrayToStream(writer, stream, data)
  class(tWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(ip), intent(in) :: data(:)

  integer(ip) :: index

  do index = 1, size(data)
    call writer%WriteInteger(stream, data(index))
  end do

end subroutine WriteIntegerArrayToStream
subroutine WriteRealVectorToStream(writer, stream, data, paddedSize)
  class(tWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  real(dp), intent(in) :: data(:)
  integer(ip), intent(in), optional :: paddedSize

  integer(ip) :: index

  do index = 1, size(data)
    call writer%WriteReal(stream, data(index))
  end do
  if (present(paddedSize)) then
    do index = size(data) + 1, paddedSize
      call writer%WriteReal(stream, 0.0_dp)
    end do
  end if

end subroutine WriteRealVectorToStream
subroutine WriteRealVectorArrayToStream(writer, stream, data, paddedSize)
  class(tWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  real(dp), intent(in) :: data(:,:)
  integer(ip), intent(in), optional :: paddedSize

  integer(ip) :: index

  do index = 1, size(data, dim=2)
    call writer%WriteRealVector(stream, data(:,index), paddedSize)
  end do

end subroutine WriteRealVectorArrayToStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Create a text writer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
function MakeTextWriter(separator) result(writer)
  character(len=*), intent(in), optional :: separator
  type(tTextWriter) :: writer

  writer%Separator = ' '
  if (present(separator)) writer%Separator = separator 

end function MakeTextWriter

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Encode data as text.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine WriteStringToStream(stream, string)
  class(tOutputStream), intent(inout) :: stream
  character(len=*), intent(in), target :: string

  integer(bp), contiguous, pointer :: bytes(:)
  
  call c_f_pointer(cptr=c_loc(string), fptr=bytes, shape=[len(string)])

  call stream%Write(bytes)

end subroutine WriteStringToStream
subroutine WriteTextByteToStream(writer, stream, data)
  class(tTextWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(bp), intent(in) :: data

  call writer%WriteInteger(stream, int(data, kind=ip))

end subroutine WriteTextByteToStream
subroutine WriteTextIntegerToStream(writer, stream, data)
  class(tTextWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(ip), intent(in) :: data

  call WriteStringToStream(stream, string=I2S(data)//writer%Separator)

end subroutine WriteTextIntegerToStream
subroutine WriteTextRealToStream(writer, stream, data)
  class(tTextWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  real(dp), intent(in) :: data

  call WriteStringToStream(stream, string=R2S(data)//writer%Separator)

end subroutine WriteTextRealToStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Create a binary writer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
function MakeBinaryWriter(endian, longInts, singleReals) result(writer)
  character(len=*), intent(in), optional :: endian
  logical, intent(in), optional :: longInts, singleReals
  type(tBinaryWriter) :: writer

  ! ----------------------
  ! Setup the endianness.
  ! ----------------------
  writer%SwapEndian = .false.
  if (present(endian)) then
    if (endian == 'big') then
#$if not BIG_ENDIAN
      writer%SwapEndian = .true.
#$end if
    else if (endian == 'little') then
#$if BIG_ENDIAN
      writer%SwapEndian = .true.
#$end if
    else
      call ErrorStop('Invalid endianness.')
    end if
  end if

  ! ----------------------
  ! Set up the data casts.
  ! ----------------------
  writer%LongInts = .false.
  if (present(longInts)) writer%LongInts = longInts

  writer%SingleReals = .false.
  if (present(singleReals)) writer%SingleReals = singleReals

end function MakeBinaryWriter

!! ----------------------------------------------------------------- !!
!! Swap endianness of the bytes.
!! ----------------------------------------------------------------- !!
subroutine SwapEndian(bytes)
  integer(bp), intent(inout) :: bytes(:)

  integer(bp) :: tmp
  integer(ip) :: i, j, n

  n = size(bytes)
  do i = 1, n/2
    j = n - i + 1; tmp = bytes(i)
    bytes(i) = bytes(j); bytes(j) = tmp
  end do

end subroutine SwapEndian

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Write data as binary bytes.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine WriteBinaryDataToStream(writer, stream, data, size)
  class(tBinaryWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  type(*), intent(in), target :: data
  integer(c_size_t), intent(in) :: size

  integer(bp), contiguous, pointer :: bytes(:)

  call c_f_pointer(cptr=c_loc(data), fptr=bytes, shape=[size])

  !! TODO: inplace swapping is a bad idea.
  if (writer%SwapEndian) call SwapEndian(bytes)
  call stream%Write(bytes)
  if (writer%SwapEndian) call SwapEndian(bytes)

end subroutine WriteBinaryDataToStream
subroutine WriteBinaryByteToStream(writer, stream, data)
  class(tBinaryWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(bp), intent(in) :: data

  call WriteBinaryDataToStream(writer, stream, data, c_sizeof(data))

end subroutine WriteBinaryByteToStream
subroutine WriteBinaryIntegerToStream(writer, stream, data)
  class(tBinaryWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(ip), intent(in) :: data

  if (writer%LongInts) then
    associate(longData => int(data, kind=lp))
      call WriteBinaryDataToStream(writer, stream, longData, c_sizeof(longData))
    end associate
  else
    call WriteBinaryDataToStream(writer, stream, data, c_sizeof(data))
  end if

end subroutine WriteBinaryIntegerToStream
subroutine WriteBinaryRealToStream(writer, stream, data)
  class(tBinaryWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  real(dp), intent(in) :: data

  if (writer%SingleReals) then
    associate(singleData => real(data, kind=sp))
      call WriteBinaryDataToStream(writer, stream, singleData, c_sizeof(singleData))
    end associate
  else
    call WriteBinaryDataToStream(writer, stream, data, c_sizeof(data))
  end if

end subroutine WriteBinaryRealToStream

end module StormRuler_IO_Writer
