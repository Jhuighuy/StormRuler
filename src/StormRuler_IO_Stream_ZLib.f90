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

use StormRuler_Libs, only: compress2

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Output stream that compresses bytes with ZLib.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tOutputStream) :: tZLibOutputStream
  class(tOutputStream), allocatable, private :: InnerStream
  integer(ip), private :: CompressionLevel
  integer(bp), allocatable, private :: Bytes(:) 
  integer(lp), private :: NumBytes

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
function MakeZLibOutputStream(innerStream, compressionLevel) result(stream)
  class(tOutputStream), intent(inout), allocatable :: innerStream
  integer(ip), intent(in), optional :: compressionLevel
  type(tZLibOutputStream) :: stream

  !! BUG: Intel Fortran bug, use `move_alloc` after it is fixed. 
  !call move_alloc(from=innerStream, to=stream%InnerStream)
  stream%InnerStream = innerStream; deallocate(innerStream)

  stream%CompressionLevel = 6
  if (present(compressionLevel)) then
    stream%CompressionLevel = compressionLevel
  end if

  stream%NumBytes = 0
  allocate(stream%Bytes(1024))

end function MakeZLibOutputStream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Write bytes to the ZLib output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine WriteToZLibStream(stream, bytes)
  class(tZLibOutputStream), intent(inout) :: stream
  integer(bp), intent(in), contiguous, target :: bytes(:)

  integer(lp) :: capacity

  ! ----------------------
  ! Accumulate bytes into the internal buffer.
  ! ----------------------
  capacity = size(stream%Bytes)
  if (capacity < (stream%NumBytes + size(bytes))) then
    stream%Bytes = [stream%Bytes, spread(1_bp, 1, capacity)]
  end if

  stream%Bytes((stream%NumBytes + 1):(stream%NumBytes + size(bytes))) = bytes
  stream%NumBytes = stream%NumBytes + size(bytes)

end subroutine WriteToZLibStream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! End writing to the ZLib output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine EndWriteToZLibStream(stream)
  class(tZLibOutputStream), intent(inout) :: stream

  integer(ip) :: errorCode
  integer(lp) :: numCompressedBytes
  integer(bp), allocatable :: bytesCompressed(:)

  ! ----------------------
  ! Allocate the destination buffer.
  ! ----------------------
  numCompressedBytes = stream%NumBytes + (stream%NumBytes + 999)/1000 + 12
  allocate(bytesCompressed(numCompressedBytes))

  ! ----------------------
  ! Compress the data.
  ! ----------------------
  errorCode = compress2(bytesCompressed, numCompressedBytes, &
    & stream%Bytes(1:stream%NumBytes), stream%NumBytes, &
    & stream%CompressionLevel)
  if (errorCode /= 0) then
    call ErrorStop('ZLib `compress2` has failed, errorCode='//I2S(errorCode))
  end if

  ! ----------------------
  ! Write compressed data to the inner stream. 
  ! ----------------------
  block
    !! TODO: header should be written externally somehow.
    class(tBinaryWriter), allocatable :: writer
    writer = tBinaryWriter()
    call stream%InnerStream%BeginWrite()
    call writer%Write(stream%InnerStream, int(1, kind=ip))
    call writer%Write(stream%InnerStream, int(stream%NumBytes, kind=ip))
    call writer%Write(stream%InnerStream, int(stream%NumBytes, kind=ip))
    call writer%Write(stream%InnerStream, int(numCompressedBytes, kind=ip))
    call stream%InnerStream%EndWrite()
  end block
  call stream%InnerStream%BeginWrite()
  call stream%InnerStream%Write(bytesCompressed(1:numCompressedBytes))
  call stream%InnerStream%EndWrite()

  ! ----------------------
  ! Reset the buffer.
  ! ----------------------
  stream%NumBytes = 0

end subroutine EndWriteToZLibStream

#$end if

end module StormRuler_IO_Stream_ZLib
