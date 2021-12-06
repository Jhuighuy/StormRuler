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
module StormRuler_IO_Base64

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: ip, i8

use StormRuler_IO_Stream, only: tOutputStream, tBinaryOutputStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Output stream that encodes output bytes with Base64.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tBinaryOutputStream) :: tBase64OutputStream
  class(tOutputStream), pointer, private :: Inner
  integer(i8), private :: Buffer(3), BuffAddr = 0

contains
  procedure, non_overridable :: Close => Base64OStreamClose
  procedure, non_overridable :: Flush => Base64OStreamFlush
  procedure, non_overridable :: WriteBytes => Base64OStreamWriteBytes

end type tBase64OutputStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Flush the Base64 output stream.
!! ----------------------------------------------------------------- !!
subroutine Base64OStreamFlush(stream)
  class(tBase64OutputStream), intent(inout) :: stream

  call stream%Inner%Flush()

end subroutine Base64OStreamFlush

!! ----------------------------------------------------------------- !!
!! Close the Base64 output stream.
!! ----------------------------------------------------------------- !!
subroutine Base64OStreamClose(stream)
  class(tBase64OutputStream), intent(inout) :: stream

  call Base64OStreamEncodeBufferedAndFlush(stream)
  call stream%Inner%Close()

end subroutine Base64OStreamClose

!! ----------------------------------------------------------------- !!
!! Write bytes to the Base64 output stream.
!! ----------------------------------------------------------------- !!
subroutine Base64OStreamWriteBytes(stream, data)
  class(tBase64OutputStream), intent(inout) :: stream
  integer(i8), intent(in), contiguous, target :: data(:)

  integer(ip) :: index

  do index = 1, size(data)
    stream%BuffAddr = stream%BuffAddr + 1
    stream%Buffer(stream%BuffAddr) = data(index)
    if (stream%BuffAddr == 3) then
      call Base64OStreamEncodeBufferedAndFlush(stream)
    end if
  end do

end subroutine Base64OStreamWriteBytes

!! ----------------------------------------------------------------- !!
!! Encode the buffered triplet with Base64 and flush it.
!! Implementation was shamelessly stolen from here:
!! https://gist.github.com/tomykaira/f0fd86b6c73063283afe550bc5d77594
!! ----------------------------------------------------------------- !!
subroutine Base64OStreamEncodeBufferedAndFlush(stream)
  class(tBase64OutputStream), intent(inout) :: stream

  integer(ip) :: x, y
  character(len=4) :: chars
  character, parameter :: tbl(0:63) = &
    [ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', & 
    & 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', &
    & 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', &
    & 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', &
    & 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', &
    & 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', &
    & 'w', 'x', 'y', 'z', '0', '1', '2', '3', &
    & '4', '5', '6', '7', '8', '9', '+', '/'  ]

  associate(b => stream%Buffer, s => stream%BuffAddr)
    if (s == 0) return
    x = ishft(iand(b(1), z'03'), 4)
    chars(1:1) = tbl(iand(ishft(b(1), -2), z'3F'))
    if (s == 1) then
      chars(2:2) = tbl(x)
      chars(3:4) = '='
    else if (s >= 2) then
      chars(2:2) = tbl(ior(x, ishft(iand(b(2), z'F0'), -4)))
      y = ishft(iand(b(2), z'0F'), 2)
      if (s == 2) then
        chars(3:3) = tbl(y)
        chars(4:4) = '='
      else if (s == 3) then
        chars(3:3) = tbl(ior(y, ishft(iand(b(3), z'C0'), -6)))
        chars(4:4) = tbl(iand(b(3), z'3F'))
      end if
    end if
    s = 0
  end associate

  print *, chars
  !call stream%Inner%WriteString(chars)

end subroutine Base64OStreamEncodeBufferedAndFlush

end module StormRuler_IO_Base64
