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
module StormRuler_IO_Stream_Base64

use StormRuler_Consts, only: bp, ip

use StormRuler_IO_Stream, only: tOutputStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Output stream that encodes bytes with Base64.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tOutputStream) :: tBase64OutputStream
  class(tOutputStream), allocatable, private :: Inner
  integer(bp), private :: Buffer(3), BufferLength = 0

contains
  procedure, non_overridable :: Write => WriteToBase64Stream
  procedure, non_overridable :: EndWrite => EndWriteToBase64Stream

end type tBase64OutputStream

interface tBase64OutputStream
  module procedure MakeBase64OutputStream
end interface tBase64OutputStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Create a Base64 output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
function MakeBase64OutputStream(innerStream) result(stream)
  class(tOutputStream), intent(inout), target :: innerStream
  type(tBase64OutputStream) :: stream

  allocate(stream%Inner, source=innerStream)

end function MakeBase64OutputStream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Write bytes to the Base64 output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine WriteToBase64Stream(stream, bytes)
  class(tBase64OutputStream), intent(inout) :: stream
  integer(bp), intent(in), contiguous, target :: bytes(:)

  integer(ip) :: index

  do index = 1, size(bytes)
    stream%BufferLength = stream%BufferLength + 1
    stream%Buffer(stream%BufferLength) = bytes(index)
    if (stream%BufferLength == 3) then
      call EncodeBytesAndWriteToStream(stream)
    end if
  end do

end subroutine WriteToBase64Stream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! End writing to the Base64 output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine EndWriteToBase64Stream(stream)
  class(tBase64OutputStream), intent(inout) :: stream

    if (stream%BufferLength > 0) then
      call EncodeBytesAndWriteToStream(stream)
    end if
  call stream%Inner%EndWrite()

end subroutine EndWriteToBase64Stream

!! ----------------------------------------------------------------- !!
!! Encode the buffered triplet with Base64 and flush it.
!! ----------------------------------------------------------------- !!
subroutine EncodeBytesAndWriteToStream(stream)
  class(tBase64OutputStream), intent(inout) :: stream

  integer(bp) :: encoded(4)

  select case(stream%BufferLength)
    case(1)
      call EncodeSingle(stream%Buffer(1))
    case(2)
      call EncodeDuplet(stream%Buffer(1), stream%Buffer(2))
    case(3)
      call EncodeTriplet(stream%Buffer(1), stream%Buffer(2), stream%Buffer(3))
  end select

  stream%BufferLength = 0
  call stream%Inner%Write(encoded)

contains
  integer(bp) function EncodeChar(e)
    integer(bp), intent(in) :: e

    character, parameter :: encodeTable(0:63) = &
      [ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', & 
      & 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', &
      & 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', &
      & 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', &
      & 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', &
      & 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', &
      & 'w', 'x', 'y', 'z', '0', '1', '2', '3', &
      & '4', '5', '6', '7', '8', '9', '+', '/'  ]

    EncodeChar = iachar(encodeTable(e), kind=bp)

  end function EncodeChar
  subroutine EncodeSingle(a)
    integer(bp), intent(in) :: a

    encoded(1) = EncodeChar(iand(ishft(a, -2), 63_bp)) ! 0x3F
    encoded(2) = EncodeChar(iand(ishft(a, +4), 48_bp)) ! 0x30
    encoded(3:4) = iachar('=')

  end subroutine EncodeSingle
  subroutine EncodeDuplet(a, b)
    integer(bp), intent(in) :: a, b

    encoded(1) = EncodeChar(    iand(ishft(a, -2), 63_bp)  ) ! 0x3F
    encoded(2) = EncodeChar(ior(iand(ishft(a, +4), 48_bp), & ! 0x30
                              & iand(ishft(b, -4), 15_bp)) ) ! 0x0F
    encoded(3) = EncodeChar(    iand(ishft(b, +2), 60_bp)  ) ! 0x3C
    encoded(4) = iachar('=')

  end subroutine EncodeDuplet
  subroutine EncodeTriplet(a, b, c)
    integer(bp), intent(in) :: a, b, c

    encoded(1) = EncodeChar(    iand(ishft(a, -2), 63_bp)  ) ! 0x3F
    encoded(2) = EncodeChar(ior(iand(ishft(a, +4), 48_bp), & ! 0x30
                              & iand(ishft(b, -4), 15_bp)) ) ! 0x0F
    encoded(3) = EncodeChar(ior(iand(ishft(b, +2), 60_bp), & ! 0x3C
                              & iand(ishft(c, -6), 03_bp)) ) ! 0x03 
    encoded(4) = EncodeChar(    iand(ishft(c,  0), 63_bp)  ) ! 0x3F

  end subroutine EncodeTriplet
end subroutine EncodeBytesAndWriteToStream

end module StormRuler_IO_Stream_Base64
