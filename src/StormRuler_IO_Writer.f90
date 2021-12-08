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
!! Abstract data encoded writer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, abstract :: tWriter
contains
  generic :: Write => &
    & WriteByte, WriteByteArray, &
    & WriteInteger, WriteIntegerArray, &
    & WriteReal, WriteRealVector, WriteRealVectorArray
  procedure(tWriteEncodedByteFunc), deferred :: WriteByte
  procedure(tWriteEncodedIntegerFunc), deferred :: WriteInteger
  procedure(tWriteEncodeRealFunc), deferred :: WriteReal
  procedure :: WriteByteArray, WriteIntegerArray
  procedure :: WriteRealVector, WriteRealVectorArray

end type tWriter

!! ----------------------------------------------------------------- !!
!! Encode and write the data.
!! ----------------------------------------------------------------- !!
abstract interface
  subroutine tWriteEncodedByteFunc(writer, stream, data)
    import :: tWriter, tOutputStream, bp
    class(tWriter), intent(inout) :: writer
    class(tOutputStream), intent(inout) :: stream
    integer(bp), intent(in) :: data
  end subroutine tWriteEncodedByteFunc
  subroutine tWriteEncodedIntegerFunc(writer, stream, data)
    import :: tWriter, tOutputStream, ip
    class(tWriter), intent(inout) :: writer
    class(tOutputStream), intent(inout) :: stream
    integer(ip), intent(in) :: data
  end subroutine tWriteEncodedIntegerFunc
  subroutine tWriteEncodeRealFunc(writer, stream, data)
    import :: tWriter, tOutputStream, dp
    class(tWriter), intent(inout) :: writer
    class(tOutputStream), intent(inout) :: stream
    real(dp), intent(in) :: data
  end subroutine tWriteEncodeRealFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! ASCII encoded writer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tWriter) :: tAsciiWriter
contains
  procedure, non_overridable :: WriteByte => WriteByteAscii
  procedure, non_overridable :: WriteInteger => WriteIntegerAscii
  procedure, non_overridable :: WriteReal => WriteRealAscii

end type tAsciiWriter

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Binary encoded writer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tWriter) :: tBinaryWriter
  logical, private :: SwapEndianness
  logical, private :: LongIntegers, SingleReals

contains
  procedure, non_overridable :: WriteByte => WriteByteBinary
  procedure, non_overridable :: WriteInteger => WriteIntegerBinary
  procedure, non_overridable :: WriteReal => WriteRealBinary

end type tBinaryWriter

interface tBinaryWriter
  module procedure MakeBinaryWriter
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Encode arrays.
!! ----------------------------------------------------------------- !!
subroutine WriteByteArray(writer, stream, data)
  class(tWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(bp), intent(in) :: data(:)

  integer(ip) :: index

  do index = 1, size(data)
    call writer%WriteByte(stream, data(index))
  end do

end subroutine WriteByteArray
subroutine WriteIntegerArray(writer, stream, data)
  class(tWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(ip), intent(in) :: data(:)

  integer(ip) :: index

  do index = 1, size(data)
    call writer%WriteInteger(stream, data(index))
  end do

end subroutine WriteIntegerArray
subroutine WriteRealVector(writer, stream, data, paddedSize)
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

end subroutine WriteRealVector
subroutine WriteRealVectorArray(writer, stream, data, paddedSize)
  class(tWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  real(dp), intent(in) :: data(:,:)
  integer(ip), intent(in), optional :: paddedSize

  integer(ip) :: index

  do index = 1, size(data, dim=2)
    call writer%WriteRealVector(stream, data(:,index), paddedSize)
  end do

end subroutine WriteRealVectorArray

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Encode data as ASCII text.
!! ----------------------------------------------------------------- !!
subroutine WriteStringAscii(stream, data)
  class(tOutputStream), intent(inout) :: stream
  character(len=*), intent(in) :: data

  integer(bp), contiguous, pointer :: bytes(:)
  
  call c_f_pointer(cptr=c_loc(data), &
    & fptr=bytes, shape=[len(data)*c_sizeof(data(1:1))])

  call stream%Write(bytes)

end subroutine WriteStringAscii
subroutine WriteByteAscii(writer, stream, data)
  class(tAsciiWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(bp), intent(in) :: data

  call writer%WriteInteger(stream, int(data, kind=ip))

end subroutine WriteByteAscii
subroutine WriteIntegerAscii(writer, stream, data)
  class(tAsciiWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(ip), intent(in) :: data

  call WriteStringAscii(stream, I2S(data)//' ')

end subroutine WriteIntegerAscii
subroutine WriteRealAscii(writer, stream, data)
  class(tAsciiWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  real(dp), intent(in) :: data

  call WriteStringAscii(stream, R2S(data)//' ')

end subroutine WriteRealAscii

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Create the binary writer.
!! ----------------------------------------------------------------- !!
function MakeBinaryWriter(endianness, longIntegers, singleReals) result(writer)
  type(tBinaryWriter) :: writer
  character(len=*), intent(in), optional :: endianness
  logical, intent(in), optional :: longIntegers, singleReals

  ! ----------------------
  ! Setup the endianness.
  ! ----------------------
  writer%SwapEndianness = .false.
  if (present(endianness)) then
    if (endianness == 'big') then
#$if not BIG_ENDIAN
      writer%SwapEndianness = .true.
#$end if
    else if (endianness == 'little') then
#$if BIG_ENDIAN
      writer%SwapEndianness = .true.
#$end if
    else
      call ErrorStop('Invalid endianness.')
    end if
  end if

  ! ----------------------
  ! Set up the data casts.
  ! ----------------------
  writer%LongIntegers = .false.
  if (present(longIntegers)) writer%LongIntegers = longIntegers

  writer%SingleReals = .false.
  if (present(singleReals)) writer%SingleReals = singleReals

end function MakeBinaryWriter

!! ----------------------------------------------------------------- !!
!! Swap endianness of the bytes.
!! ----------------------------------------------------------------- !!
subroutine SwapEndianness(bytes)
  integer(bp), intent(inout) :: bytes(:)

  integer(bp) :: tmp
  integer(ip) :: i, j, n

  n = size(bytes)
  do i = 1, n/2
    j = n - i + 1; tmp = bytes(i)
    bytes(i) = bytes(j); bytes(j) = tmp
  end do

end subroutine SwapEndianness

!! ----------------------------------------------------------------- !!
!! Write data as binary bytes.
!! ----------------------------------------------------------------- !!
subroutine WriteBinaryData(writer, stream, data, size)
  class(tBinaryWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  type(*), intent(in), target :: data
  integer(c_size_t), intent(in) :: size

  integer(bp), contiguous, pointer :: bytes(:)

  call c_f_pointer(cptr=c_loc(data), fptr=bytes, shape=[size])

  if (writer%SwapEndianness) call SwapEndianness(bytes)
  call stream%Write(bytes)
  if (writer%SwapEndianness) call SwapEndianness(bytes)

end subroutine WriteBinaryData
subroutine WriteByteBinary(writer, stream, data)
  class(tBinaryWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(bp), intent(in) :: data

  call WriteBinaryData(writer, stream, data, c_sizeof(data))

end subroutine WriteByteBinary
subroutine WriteIntegerBinary(writer, stream, data)
  class(tBinaryWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  integer(ip), intent(in) :: data

  if (writer%LongIntegers) then
    associate(longData => int(data, kind=lp))
      call WriteBinaryData(writer, stream, longData, c_sizeof(longData))
    end associate
  else
    call WriteBinaryData(writer, stream, data, c_sizeof(data))
  end if

end subroutine WriteIntegerBinary
subroutine WriteRealBinary(writer, stream, data)
  class(tBinaryWriter), intent(inout) :: writer
  class(tOutputStream), intent(inout) :: stream
  real(dp), intent(in) :: data

  if (writer%SingleReals) then
    associate(singleData => real(data, kind=sp))
      call WriteBinaryData(writer, stream, singleData, c_sizeof(singleData))
    end associate
  else
    call WriteBinaryData(writer, stream, data, c_sizeof(data))
  end if

end subroutine WriteRealBinary

end module StormRuler_IO_Writer
