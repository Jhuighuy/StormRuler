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
module StormRuler_IO_Stream_ZLib

#$use 'StormRuler_Macros.fi'
#$if HAS_ZLIB

use StormRuler_Consts, only: bp, ip, lp
  
use StormRuler_Helpers, only: ErrorStop, I2S

use StormRuler_IO_Writer, only: tBinaryWriter
use StormRuler_IO_Stream, only: tOutputStream

use, intrinsic :: iso_c_binding, only: c_int8_t, c_int, c_long

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Output stream that compresses bytes with ZLib.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tOutputStream) :: tZLibOutputStream
  class(tOutputStream), allocatable, private :: Inner
  integer(bp), allocatable, private :: Buffer(:) 
  integer(lp), private :: BufferSize = 0

contains
  procedure, non_overridable :: Write => WriteToZLibStream
  procedure, non_overridable :: EndWrite => EndWriteToZLibStream

end type tZLibOutputStream

interface tZLibOutputStream
  module procedure MakeZLibOutputStream
end interface tZLibOutputStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Create a ZLib output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
function MakeZLibOutputStream(innerStream) result(stream)
  class(tOutputStream), intent(inout), target :: innerStream
  type(tZLibOutputStream) :: stream

  allocate(stream%Inner, source=innerStream)

end function MakeZLibOutputStream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Write bytes to the ZLib output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine WriteToZLibStream(stream, bytes)
  class(tZLibOutputStream), intent(inout) :: stream
  integer(bp), intent(in), contiguous, target :: bytes(:)

  !! TODO: implement appending with preallocation.
  stream%BufferSize = stream%BufferSize + size(bytes)
  if (allocated(stream%Buffer)) then
    stream%Buffer = [stream%Buffer, bytes]
  else
    stream%Buffer = bytes
  end if

end subroutine WriteToZLibStream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! End writing to the ZLib output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine EndWriteToZLibStream(stream)
  class(tZLibOutputStream), intent(inout) :: stream

  interface
    function ZLibCompress2(dest, destSize, &
        & source, sourceSize, level) result(result) bind(C, name='compress2')
      import :: c_int8_t, c_int, c_long
      integer(c_int8_t), intent(inout) :: dest(*)
      integer(c_long), intent(inout) :: destSize
      integer(c_int8_t), intent(in) :: source(*)
      integer(c_long), intent(in), value :: sourceSize
      integer(c_int), intent(in), value :: level
      integer(c_int) :: result
    end function ZLibCompress2
  end interface

  integer(ip) :: result
  integer(lp) :: destBufferSize
  integer(bp), allocatable :: destBuffer(:)
  class(tBinaryWriter), allocatable :: writer

  ! ----------------------
  ! Allocate the destination buffer.
  ! ----------------------
  destBufferSize = stream%BufferSize + &
    & (stream%BufferSize + 999)/1000 + 12
  allocate(destBuffer(destBufferSize))
  writer = tBinaryWriter()

  ! ----------------------
  ! Compress the data.
  ! ----------------------
  result = ZLibCompress2(destBuffer, destBufferSize, &
    & stream%Buffer, stream%BufferSize, 7_lp)
  if (result /= 0) then
    call ErrorStop('ZLib `compress2` has failed, result='//I2S(result))
  end if

  ! ----------------------
  ! Write compressed data to the inner stream. 
  ! ----------------------
  call stream%Inner%BeginWrite()
  call writer%Write(stream%Inner, int(1, kind=ip))
  call writer%Write(stream%Inner, int(stream%BufferSize, kind=ip))
  call writer%Write(stream%Inner, int(stream%BufferSize, kind=ip))
  call writer%Write(stream%Inner, int(destBufferSize, kind=ip))
  call stream%Inner%EndWrite()
  call stream%Inner%BeginWrite()
  call stream%Inner%Write(destBuffer(1:destBufferSize))
  call stream%Inner%EndWrite()
  stream%BufferSize = 0
  deallocate(stream%Buffer)

end subroutine EndWriteToZLibStream

#$end if

end module StormRuler_IO_Stream_ZLib
