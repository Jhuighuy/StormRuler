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
module StormRuler_Array

use StormRuler_Consts, only: ip, dp

use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

#$use 'StormRuler_Macros.fi'
#$let NUM_ARGUMENTS = 10

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Fortran-style array with inplace reshapes.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
type :: tArray
  ! ----------------------
  ! Contiguous data contained in the array.
  ! ----------------------
  real(dp), contiguous, pointer :: mData(:) => null()
  ! ----------------------
  ! Shape of the array.
  ! ----------------------
  integer(ip), contiguous, pointer :: mShape(:) => null()

contains
  procedure :: Free => FreeArray$1

  procedure :: Rank => ArrayRank

  procedure :: Slice => ArraySlice

#$do N = 1, NUM_RANKS
  generic :: Get => Get$N
  procedure :: Get$N => GetArrayData$N
#$end do
end type tArray

interface AllocArray
  module procedure AllocArrayShape
#$do N = 1, NUM_ARGUMENTS
  module procedure AllocArrayMold$N
#$end do
end interface AllocArray

interface FreeArray
#$do N = 1, NUM_ARGUMENTS
  module procedure FreeArray$N
#$end do
end interface FreeArray

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Allocate a contiguous array with specified shape.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine AllocArrayShape(array, shape)
  class(tArray), intent(inout) :: array
  integer(ip), intent(in) :: shape(:)

  allocate(array%mShape, mold=shape); array%mShape = shape
  allocate(array%mData(product(array%mShape)))

end subroutine AllocArrayShape

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Allocate a contiguous array a mold.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do N = 1, NUM_ARGUMENTS
subroutine AllocArrayMold$N(@{array$$}@, mold)
  class(tArray), intent(inout) :: @{array$$}@
  class(tArray), intent(in) :: mold

#$do I = 1, N
  call AllocArrayShape(array$I, mold%mShape) 
#$end do

end subroutine AllocArrayMold$N
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Free an array. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do N = 1, NUM_ARGUMENTS
subroutine FreeArray$N(@{array$$}@)
  class(tArray), intent(inout) :: @{array$$}@

#$do I = 1, N
  deallocate(array$I%mShape)
  deallocate(array$I%mData)
#$end do

end subroutine FreeArray$N
#$end do

pure logical function ArrayAllocated(array)
  class(tArray), intent(in) :: array

  ArrayAllocated = associated(array%mShape)
end function ArrayAllocated

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Get array rank.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
pure integer(ip) function ArrayRank(array)
  class(tArray), intent(in) :: array

  ArrayRank = size(array%mShape)

end function ArrayRank

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Get a slice of the array.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type(tArray) function ArraySlice(array, index) result(slice)
  class(tArray), intent(in) :: array
  integer(ip) :: index

  associate(rank => ArrayRank(array))

    slice%mShape => array%mShape(1:(rank - 1))

    associate(block => product(slice%mShape))
      slice%mData => array%mData((block*(index - 1) + 1):(block*index))
    end associate

  end associate

end function ArraySlice

!! ----------------------------------------------------------------- !!
!! Transform the shape.
!! ----------------------------------------------------------------- !!
pure function Rerank(shape, rank)
  integer(ip), intent(in) :: shape(:), rank
  integer(ip) :: Rerank(rank)

  associate(delta => rank - size(shape))
    if (delta >= 0) then
      Rerank(:delta) = 1
      Rerank(delta + 1:) = shape(:)
    else
      Rerank(1) = product(shape(:1 - delta))
      Rerank(2:) = shape(2 - delta:)
    end if
  end associate

end function Rerank

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Get array data as a Fortran array pointer.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do fRank = 1, NUM_RANKS
subroutine GetArrayData$fRank(array, fData)
  class(tArray), intent(in) :: array
  real(dp), intent(out), pointer :: fData(@:)

  integer(ip) :: rank
  integer(ip) :: fShape($fRank)

  rank = size(array%mShape)

  fShape = Rerank(array%mShape, $fRank)
  call c_f_pointer(cptr=c_loc(array%mData), fptr=fData, shape=fShape)

end subroutine GetArrayData$fRank
#$end do

end module StormRuler_Array
