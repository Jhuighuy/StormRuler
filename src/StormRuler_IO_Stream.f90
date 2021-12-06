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
module StormRuler_IO_Stream

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8

use StormRuler_Helpers, only: L2S, I2S, R2S

use, intrinsic :: iso_c_binding, only: c_size_t, &
  & c_loc, c_f_pointer, c_sizeof

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Abstract output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, abstract :: tOutputStream
contains
  procedure(tFlushOStreamFunc), deferred :: Flush
  procedure(tCloseOStreamFunc), deferred :: Close

  generic :: Write => WriteBytes, WriteString, &
    & WriteLogicals, WriteIntegers, WriteReals
  procedure(tOStreamWriteBytesFunc), deferred :: WriteBytes
  procedure(tOStreamWriteStringFunc), deferred :: WriteString
  procedure(tOStreamWriteLogicalsFunc), deferred :: WriteLogicals
  procedure(tOStreamWriteIntegersFunc), deferred :: WriteIntegers
  procedure(tOStreamWriteRealsFunc), deferred :: WriteReals

end type tOutputStream

!! ----------------------------------------------------------------- !!
!! Flush the output stream.
!! ----------------------------------------------------------------- !!
abstract interface
  subroutine tFlushOStreamFunc(stream)
    import :: tOutputStream
    class(tOutputStream), intent(inout) :: stream
  end subroutine tFlushOStreamFunc
end interface

!! ----------------------------------------------------------------- !!
!! Close the output stream.
!! ----------------------------------------------------------------- !!
abstract interface
  subroutine tCloseOStreamFunc(stream)
    import :: tOutputStream
    class(tOutputStream), intent(inout) :: stream
  end subroutine tCloseOStreamFunc
end interface

!! ----------------------------------------------------------------- !!
!! Write data to the output stream.
!! ----------------------------------------------------------------- !!
abstract interface
  subroutine tOStreamWriteBytesFunc(stream, data)
    import :: tOutputStream, i8
    class(tOutputStream), intent(inout) :: stream
    integer(i8), intent(in), contiguous, target :: data(:)
  end subroutine tOStreamWriteBytesFunc
  subroutine tOStreamWriteStringFunc(stream, data)
    import :: tOutputStream
    class(tOutputStream), intent(inout) :: stream
    character(len=*), intent(in), target :: data
  end subroutine tOStreamWriteStringFunc
  subroutine tOStreamWriteLogicalsFunc(stream, data, separator)
    import :: tOutputStream
    class(tOutputStream), intent(inout) :: stream
    logical, intent(in), contiguous, target :: data(:)
    character(len=*), intent(in), optional :: separator
  end subroutine tOStreamWriteLogicalsFunc
  subroutine tOStreamWriteIntegersFunc(stream, data, separator)
    import :: tOutputStream, ip
    class(tOutputStream), intent(inout) :: stream
    integer(ip), intent(in), contiguous, target :: data(:)
    character(len=*), intent(in), optional :: separator
  end subroutine tOStreamWriteIntegersFunc
  subroutine tOStreamWriteRealsFunc(stream, data, separator)
    import :: tOutputStream, dp
    class(tOutputStream), intent(inout) :: stream
    real(dp), intent(in), contiguous, target :: data(:)
    character(len=*), intent(in), optional :: separator
  end subroutine tOStreamWriteRealsFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Abstract output stream that converts output data into to the
!! textual representation and writes it as string.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tOutputStream), abstract :: tTextOutputStream
contains
  procedure, non_overridable :: WriteString => TextOStreamWriteString
  procedure, non_overridable :: WriteLogicals => TextOStreamWriteLogicals
  procedure, non_overridable :: WriteIntegers => TextOStreamWriteIntegers
  procedure, non_overridable :: WriteReals => TextOStreamWriteReals

end type tTextOutputStream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Abstract output stream that writes output data in the
!! binary representation.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tOutputStream), abstract :: tBinaryOutputStream
contains
  procedure, non_overridable :: WriteString => BinaryOStreamWriteString
  procedure, non_overridable :: WriteLogicals => BinaryOStreamWriteLogicals
  procedure, non_overridable :: WriteIntegers => BinaryOStreamWriteIntegers
  procedure, non_overridable :: WriteReals => BinaryOStreamWriteReals

end type tBinaryOutputStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Write data to the output stream as raw bytes.
!! ----------------------------------------------------------------- !!
subroutine OStreamWriteStringAsBytes(stream, data)
  class(tOutputStream), intent(inout) :: stream
  character(len=*), intent(in), target :: data

  integer(i8), contiguous, pointer :: dataBytes(:)

  associate(numBytes => c_sizeof(data(1:1))*len(data))
    call c_f_pointer(cptr=c_loc(data), fptr=dataBytes, shape=[numBytes])
  end associate
  call stream%WriteBytes(dataBytes)

end subroutine OStreamWriteStringAsBytes
subroutine OStreamWriteDataAsBytes(stream, data, sizeof)
  class(tOutputStream), intent(inout) :: stream
  type(*), intent(in), contiguous, target :: data(:)
  integer(c_size_t), intent(in) :: sizeof

  integer(i8), contiguous, pointer :: dataBytes(:)

  associate(numBytes => sizeof*size(data))
    call c_f_pointer(cptr=c_loc(data), fptr=dataBytes, shape=[numBytes])
  end associate
  call stream%WriteBytes(dataBytes)

end subroutine OStreamWriteDataAsBytes

!! ----------------------------------------------------------------- !!
!! Write data to the text output stream.
!! ----------------------------------------------------------------- !!
subroutine TextOStreamWriteString(stream, data)
  class(tTextOutputStream), intent(inout) :: stream
  character(len=*), intent(in), target :: data

  call OStreamWriteStringAsBytes(stream, data)

end subroutine TextOStreamWriteString
subroutine TextOStreamWriteLogicals(stream, data, separator)
  class(tTextOutputStream), intent(inout) :: stream
  logical, intent(in), contiguous, target :: data(:)
  character(len=*), intent(in), optional :: separator

  call stream%WriteString(L2S(data, separator))

end subroutine TextOStreamWriteLogicals
subroutine TextOStreamWriteIntegers(stream, data, separator)
  class(tTextOutputStream), intent(inout) :: stream
  integer(ip), intent(in), contiguous, target :: data(:)
  character(len=*), intent(in), optional :: separator

  call stream%WriteString(I2S(data, separator))

end subroutine TextOStreamWriteIntegers
subroutine TextOStreamWriteReals(stream, data, separator)
  class(tTextOutputStream), intent(inout) :: stream
  real(dp), intent(in), contiguous, target :: data(:)
  character(len=*), intent(in), optional :: separator

  call stream%WriteString(R2S(data, separator))

end subroutine TextOStreamWriteReals

!! ----------------------------------------------------------------- !!
!! Write data to the binary output stream.
!! ----------------------------------------------------------------- !!
subroutine BinaryOStreamWriteString(stream, data)
  class(tBinaryOutputStream), intent(inout) :: stream
  character(len=*), intent(in), target :: data

  call OStreamWriteStringAsBytes(stream, data)

end subroutine BinaryOStreamWriteString
subroutine BinaryOStreamWriteLogicals(stream, data, separator)
  class(tBinaryOutputStream), intent(inout) :: stream
  logical, intent(in), contiguous, target :: data(:)
  character(len=*), intent(in), optional :: separator

  call OStreamWriteDataAsBytes(stream, data, c_sizeof(data(1)))

end subroutine BinaryOStreamWriteLogicals
subroutine BinaryOStreamWriteIntegers(stream, data, separator)
  class(tBinaryOutputStream), intent(inout) :: stream
  integer(ip), intent(in), contiguous, target :: data(:)
  character(len=*), intent(in), optional :: separator

  call OStreamWriteDataAsBytes(stream, data, c_sizeof(data(1)))

end subroutine BinaryOStreamWriteIntegers
subroutine BinaryOStreamWriteReals(stream, data, separator)
  class(tBinaryOutputStream), intent(inout) :: stream
  real(dp), intent(in), contiguous, target :: data(:)
  character(len=*), intent(in), optional :: separator

  call OStreamWriteDataAsBytes(stream, data, c_sizeof(data(1)))

end subroutine BinaryOStreamWriteReals

end module StormRuler_IO_Stream
