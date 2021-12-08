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
module StormRuler_IO_Encoder

use StormRuler_Consts, only: bp, ip, lp, sp, dp

use StormRuler_Helpers, only: ErrorStop, L2S, I2S, R2S

use, intrinsic :: iso_fortran_env, only: output_unit
use, intrinsic :: iso_c_binding, only: c_size_t, c_long, &
  & c_loc, c_f_pointer, c_sizeof

#$use 'StormRuler_Macros.fi'

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Abstract output stream.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, abstract :: tEncoder
contains
  generic :: Encode => &
    & EncodeByte, EncodeByteArray, &
    & EncodeInteger, EncodeIntegerArray, &
    & EncodeReal, EncodeRealVector, EncodeRealVectorArray
  procedure(tEncodeByteFunc), deferred :: EncodeByte
  procedure(tEncodeIntegerFunc), deferred :: EncodeInteger
  procedure(tEncodeRealFunc), deferred :: EncodeReal
  procedure :: EncodeByteArray, EncodeIntegerArray
  procedure :: EncodeRealVector, EncodeRealVectorArray
  procedure :: Leftover => EncodingLeftover

end type tEncoder

!! ----------------------------------------------------------------- !!
!! Encode the data.
!! ----------------------------------------------------------------- !!
abstract interface
  function tEncodeByteFunc(encoder, data) result(bytes)
    import :: tEncoder, bp
    class(tEncoder), intent(inout) :: encoder
    integer(bp), intent(in) :: data
    character(len=:), allocatable :: bytes
  end function tEncodeByteFunc
  function tEncodeIntegerFunc(encoder, data) result(bytes)
    import :: tEncoder, ip
    class(tEncoder), intent(inout) :: encoder
    integer(ip), intent(in) :: data
    character(len=:), allocatable :: bytes
  end function tEncodeIntegerFunc
  function tEncodeRealFunc(encoder, data) result(bytes)
    import :: tEncoder, dp
    class(tEncoder), intent(inout) :: encoder
    real(dp), intent(in) :: data
    character(len=:), allocatable :: bytes
  end function tEncodeRealFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! ASCII text encoder. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tEncoder) :: tAsciiEncoder
contains
  procedure, non_overridable :: EncodeByte => AsciiEncodeByte
  procedure, non_overridable :: EncodeInteger => AsciiEncodeInteger
  procedure, non_overridable :: EncodeReal => AsciiEncodeReal

end type tAsciiEncoder

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Data as bytes encoder with a support for endianness swap and casts
!! from double precision reals to the single precision real and 
!! from integers to long integers.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, extends(tEncoder) :: tBinaryEncoder
  logical, private :: SwapEndianness
  logical, private :: longIntegers, SingleReals

contains
  procedure, non_overridable :: EncodeByte => BinaryEncodeByte
  procedure, non_overridable :: EncodeInteger => BinaryEncodeInteger
  procedure, non_overridable :: EncodeReal => BinaryEncodeReal

end type tBinaryEncoder

interface tBinaryEncoder
  module procedure MakeBinaryEncoder
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Encode arrays.
!! ----------------------------------------------------------------- !!
function EncodeByteArray(encoder, data) result(bytes)
  class(tEncoder), intent(inout) :: encoder
  integer(bp), intent(in) :: data(:)
  character(len=:), allocatable :: bytes

  integer(ip) :: index

  bytes = encoder%EncodeByte(data(1))
  do index = 2, size(data)
    bytes = bytes//encoder%EncodeByte(data(index))
  end do

end function EncodeByteArray
function EncodeIntegerArray(encoder, data) result(bytes)
  class(tEncoder), intent(inout) :: encoder
  integer(ip), intent(in) :: data(:)
  character(len=:), allocatable :: bytes

  integer(ip) :: index

  bytes = encoder%EncodeInteger(data(1))
  do index = 2, size(data)
    bytes = bytes//encoder%EncodeInteger(data(index))
  end do

end function EncodeIntegerArray
function EncodeRealVector(encoder, data, paddedSize) result(bytes)
  class(tEncoder), intent(inout) :: encoder
  real(dp), intent(in) :: data(:)
  integer(ip), intent(in), optional :: paddedSize
  character(len=:), allocatable :: bytes

  integer(ip) :: index

  bytes = encoder%EncodeReal(data(1))
  do index = 2, size(data)
    bytes = bytes//encoder%EncodeReal(data(index))
  end do
  if (present(paddedSize)) then
    do index = size(data) + 1, paddedSize
      bytes = bytes//encoder%EncodeReal(0.0_dp)
    end do
  end if

end function EncodeRealVector
function EncodeRealVectorArray(encoder, data, paddedSize) result(bytes)
  class(tEncoder), intent(inout) :: encoder
  real(dp), intent(in) :: data(:,:)
  integer(ip), intent(in), optional :: paddedSize
  character(len=:), allocatable :: bytes

  integer(ip) :: index

  bytes = encoder%EncodeRealVector(data(:,1), paddedSize)
  do index = 2, size(data, dim=2)
    bytes = bytes//encoder%EncodeRealVector(data(:,index), paddedSize)
  end do

end function EncodeRealVectorArray

!! ----------------------------------------------------------------- !!
!! Output the encoding leftover.
!! ----------------------------------------------------------------- !!
function EncodingLeftover(encoder) result(bytes)
  class(tEncoder), intent(inout) :: encoder
  character(len=:), allocatable :: bytes

  bytes = '' ! no leftover by default.

end function EncodingLeftover

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Encode data as ASCII text.
!! ----------------------------------------------------------------- !!
function AsciiEncodeByte(encoder, data) result(bytes)
  class(tAsciiEncoder), intent(inout) :: encoder
  integer(bp), intent(in) :: data
  character(len=:), allocatable :: bytes

  bytes = encoder%EncodeInteger(int(data, kind=ip))

end function AsciiEncodeByte
function AsciiEncodeInteger(encoder, data) result(bytes)
  class(tAsciiEncoder), intent(inout) :: encoder
  integer(ip), intent(in) :: data
  character(len=:), allocatable :: bytes

  bytes = I2S(data)//' '

end function AsciiEncodeInteger
function AsciiEncodeReal(encoder, data) result(bytes)
  class(tAsciiEncoder), intent(inout) :: encoder
  real(dp), intent(in) :: data
  character(len=:), allocatable :: bytes

  bytes = R2S(data)//' '

end function AsciiEncodeReal

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Create the binary encoder.
!! ----------------------------------------------------------------- !!
function MakeBinaryEncoder(endianness, longIntegers, singleReals) result(encoder)
  type(tBinaryEncoder) :: encoder
  character(len=*), intent(in), optional :: endianness
  logical, intent(in), optional :: longIntegers, singleReals

  ! ----------------------
  ! Setup the endianness.
  ! ----------------------
  encoder%SwapEndianness = .false.
  if (present(endianness)) then
    if (endianness == 'big') then
#$if not BIG_ENDIAN
      encoder%SwapEndianness = .true.
#$end if
    else if (endianness == 'little') then
#$if BIG_ENDIAN
      encoder%SwapEndianness = .true.
#$end if
    else
      call ErrorStop('Invalid endianness.')
    end if
  end if

  ! ----------------------
  ! Set up the data casts.
  ! ----------------------
  encoder%longIntegers = .false.
  if (present(longIntegers)) encoder%longIntegers = longIntegers

  encoder%SingleReals = .false.
  if (present(singleReals)) encoder%SingleReals = singleReals

end function MakeBinaryEncoder

!! ----------------------------------------------------------------- !!
!! Swap endianness of the bytes.
!! ----------------------------------------------------------------- !!
subroutine SwapEndianness(bytes)
  character(len=*), intent(inout) :: bytes

  character :: tmp
  integer(ip) :: i, j, n

  n = len(bytes)
  do i = 1, n/2
    j = n - i + 1; tmp = bytes(i:i)
    bytes(i:i) = bytes(j:j); bytes(j:j) = tmp
  end do

end subroutine SwapEndianness

!! ----------------------------------------------------------------- !!
!! Encode data as raw bytes.
!! ----------------------------------------------------------------- !!
function BinaryEncodeData(encoder, data, size) result(bytes)
  class(tBinaryEncoder), intent(inout) :: encoder
  type(*), intent(in), target :: data
  integer(c_size_t), intent(in) :: size
  character(len=:), allocatable :: bytes

  character(len=size), pointer :: dataAsBytes

  call c_f_pointer(cptr=c_loc(data), fptr=dataAsBytes)

  bytes = dataAsBytes
  if (encoder%SwapEndianness) call SwapEndianness(bytes)

end function BinaryEncodeData
function BinaryEncodeByte(encoder, data) result(bytes)
  class(tBinaryEncoder), intent(inout) :: encoder
  integer(bp), intent(in) :: data
  character(len=:), allocatable :: bytes

  bytes = BinaryEncodeData(encoder, data, c_sizeof(data))

end function BinaryEncodeByte
function BinaryEncodeInteger(encoder, data) result(bytes)
  class(tBinaryEncoder), intent(inout) :: encoder
  integer(ip), intent(in) :: data
  character(len=:), allocatable :: bytes

  if (encoder%longIntegers) then
    associate(longData => int(data, kind=lp))
      bytes = BinaryEncodeData(encoder, longData, c_sizeof(longData))
    end associate
  else
    bytes = BinaryEncodeData(encoder, data, c_sizeof(data))
  end if

end function BinaryEncodeInteger
function BinaryEncodeReal(encoder, data) result(bytes)
  class(tBinaryEncoder), intent(inout) :: encoder
  real(dp), intent(in) :: data
  character(len=:), allocatable :: bytes

  if (encoder%SingleReals) then
    associate(singleData => real(data, kind=sp))
      bytes = BinaryEncodeData(encoder, singleData, c_sizeof(singleData))
    end associate
  else
    bytes = BinaryEncodeData(encoder, data, c_sizeof(data))
  end if

end function BinaryEncodeReal

end module StormRuler_IO_Encoder
