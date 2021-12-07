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

use, intrinsic :: iso_fortran_env, only: output_unit
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
  generic :: Write => WriteBytes, WriteString, &
    & WriteLogicals, WriteIntegers, WriteReals
  procedure(tOStreamWriteBytesFunc), deferred :: WriteBytes
  procedure(tOStreamWriteStringFunc), deferred :: WriteString
  procedure(tOStreamWriteLogicalsFunc), deferred :: WriteLogicals
  procedure(tOStreamWriteIntegersFunc), deferred :: WriteIntegers
  procedure(tOStreamWriteRealsFunc), deferred :: WriteReals

  procedure(tFinalizeOStreamFunc), deferred :: Finalize

end type tOutputStream

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
  subroutine tOStreamWriteLogicalsFunc(stream, data)
    import :: tOutputStream
    class(tOutputStream), intent(inout) :: stream
    logical, intent(in), contiguous, target :: data(:)
  end subroutine tOStreamWriteLogicalsFunc
  subroutine tOStreamWriteIntegersFunc(stream, data)
    import :: tOutputStream, ip
    class(tOutputStream), intent(inout) :: stream
    integer(ip), intent(in), contiguous, target :: data(:)
  end subroutine tOStreamWriteIntegersFunc
  subroutine tOStreamWriteRealsFunc(stream, data)
    import :: tOutputStream, dp
    class(tOutputStream), intent(inout) :: stream
    real(dp), intent(in), contiguous, target :: data(:)
  end subroutine tOStreamWriteRealsFunc
end interface

!! ----------------------------------------------------------------- !!
!! Finalize the output stream.
!! ----------------------------------------------------------------- !!
abstract interface
  subroutine tFinalizeOStreamFunc(stream)
    import :: tOutputStream
    class(tOutputStream), intent(inout) :: stream
  end subroutine tFinalizeOStreamFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Abstract output stream that converts output data into to the
!! textual representation and writes it as string.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tOutputStream), abstract :: tTextOutputStream
contains
  procedure, non_overridable :: WriteBytes => TextOStreamWriteBytes
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

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Text file output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tTextOutputStream) :: tTextFileOutputStream
  integer(ip), private :: Unit

contains
  procedure :: Init => InitTextFileOStream
  procedure, non_overridable :: WriteString => TextFileOStreamWriteString
  procedure, non_overridable :: Finalize => FinalizeTextFileOStream

end type tTextFileOutputStream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Binary file output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tBinaryOutputStream) :: tBinaryFileOutputStream
  integer(ip), private :: Unit

contains
  procedure :: Init => InitBinaryFileOStream
  procedure, non_overridable :: WriteBytes => BinaryFileOStreamWriteBytes
  procedure, non_overridable :: Finalize => FinalizeBinaryFileOStream

end type tBinaryFileOutputStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Write data to the text output stream.
!! ----------------------------------------------------------------- !!
subroutine TextOStreamWriteBytes(stream, data)
  class(tTextOutputStream), intent(inout) :: stream
  integer(i8), intent(in), contiguous, target :: data(:)

  call stream%WriteIntegers(int(data, kind=ip))

end subroutine TextOStreamWriteBytes
subroutine TextOStreamWriteLogicals(stream, data)
  class(tTextOutputStream), intent(inout) :: stream
  logical, intent(in), contiguous, target :: data(:)

  call stream%WriteString(L2S(data)//' ')

end subroutine TextOStreamWriteLogicals
subroutine TextOStreamWriteIntegers(stream, data)
  class(tTextOutputStream), intent(inout) :: stream
  integer(ip), intent(in), contiguous, target :: data(:)

  call stream%WriteString(I2S(data)//' ')

end subroutine TextOStreamWriteIntegers
subroutine TextOStreamWriteReals(stream, data)
  class(tTextOutputStream), intent(inout) :: stream
  real(dp), intent(in), contiguous, target :: data(:)

  call stream%WriteString(R2S(data)//' ')

end subroutine TextOStreamWriteReals

!! ----------------------------------------------------------------- !!
!! Write data to the binary output stream.
!! ----------------------------------------------------------------- !!
subroutine BinaryOStreamWriteString(stream, data)
  class(tBinaryOutputStream), intent(inout) :: stream
  character(len=*), intent(in), target :: data

  character :: char
  integer(i8), contiguous, pointer :: dataBytes(:)

  associate(numBytes => c_sizeof(char)*len(data))
    call c_f_pointer(cptr=c_loc(data), fptr=dataBytes, shape=[numBytes])
  end associate
  call stream%WriteBytes(dataBytes)

end subroutine BinaryOStreamWriteString
subroutine BinaryOStreamWriteAsBytes(stream, data, sizeof)
  class(tBinaryOutputStream), intent(inout) :: stream
  type(*), intent(in), contiguous, target :: data(:)
  integer(c_size_t), intent(in) :: sizeof

  integer(i8), contiguous, pointer :: dataBytes(:)

  associate(numBytes => sizeof*size(data))
    call c_f_pointer(cptr=c_loc(data), fptr=dataBytes, shape=[numBytes])
  end associate
  call stream%WriteBytes(dataBytes)

end subroutine BinaryOStreamWriteAsBytes
subroutine BinaryOStreamWriteLogicals(stream, data)
  class(tBinaryOutputStream), intent(inout) :: stream
  logical, intent(in), contiguous, target :: data(:)

  call BinaryOStreamWriteAsBytes(stream, data, c_sizeof(data(1)))

end subroutine BinaryOStreamWriteLogicals
subroutine BinaryOStreamWriteIntegers(stream, data)
  class(tBinaryOutputStream), intent(inout) :: stream
  integer(ip), intent(in), contiguous, target :: data(:)

  call BinaryOStreamWriteAsBytes(stream, data, c_sizeof(data(1)))

end subroutine BinaryOStreamWriteIntegers
subroutine BinaryOStreamWriteReals(stream, data)
  class(tBinaryOutputStream), intent(inout) :: stream
  real(dp), intent(in), contiguous, target :: data(:)

  call BinaryOStreamWriteAsBytes(stream, data, c_sizeof(data(1)))

end subroutine BinaryOStreamWriteReals

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Initialize the file output stream.
!! ----------------------------------------------------------------- !!
subroutine InitTextFileOStream(stream, unit)
  class(tTextFileOutputStream), intent(inout) :: stream
  integer(ip), intent(in), optional :: unit

  if (present(unit)) then
    stream%Unit = unit
  else
    stream%Unit = output_unit
  end if

end subroutine InitTextFileOStream
subroutine InitBinaryFileOStream(stream, unit)
  class(tBinaryFileOutputStream), intent(inout) :: stream
  integer(ip), intent(in) :: unit

  stream%Unit = unit

end subroutine InitBinaryFileOStream

!! ----------------------------------------------------------------- !!
!! Write a data to the file output stream.
!! ----------------------------------------------------------------- !!
subroutine TextFileOStreamWriteString(stream, data)
  class(tTextFileOutputStream), intent(inout) :: stream
  character(len=*), intent(in), target :: data

  write(unit=stream%Unit) data

end subroutine TextFileOStreamWriteString
subroutine BinaryFileOStreamWriteBytes(stream, data)
  class(tBinaryFileOutputStream), intent(inout) :: stream
  integer(i8), intent(in), contiguous, target :: data(:)

  write(unit=stream%Unit) data

end subroutine BinaryFileOStreamWriteBytes

!! ----------------------------------------------------------------- !!
!! Finalize the file output stream.
!! ----------------------------------------------------------------- !!
subroutine FinalizeTextFileOStream(stream)
  class(tTextFileOutputStream), intent(inout) :: stream

end subroutine FinalizeTextFileOStream
subroutine FinalizeBinaryFileOStream(stream)
  class(tBinaryFileOutputStream), intent(inout) :: stream

end subroutine FinalizeBinaryFileOStream

end module StormRuler_IO_Stream
