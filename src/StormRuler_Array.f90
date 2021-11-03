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

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Fortran-style array with inplace reshapes.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
type :: tArrayR
  ! ----------------------
  ! Contiguous data contained in the array.
  ! ----------------------
  real(dp), contiguous, pointer :: mData(:)
  ! ----------------------
  ! Shape of the array.
  ! ----------------------
  integer(ip), contiguous, pointer :: mShape(:)

contains
  generic :: Alloc => AllocShape, AllocMold  
  procedure :: AllocShape => AllocArray
  procedure :: AllocMold => AllocArrayMold$1

  procedure :: Free => FreeArray$1

  procedure :: Rank => ArrayRank

#$do N = 1, NUM_RANKS
  generic :: Get => Get$N 
  procedure :: Get$N => GetArrayData$N
#$end do
end type tArrayR

interface AllocArray
  module procedure AllocArray
#$do N = 1, NUM_RANKS
  module procedure AllocArrayMold$N
#$end do
end interface AllocArray

interface FreeArray
#$do N = 1, NUM_RANKS
  module procedure FreeArray$N
#$end do
end interface FreeArray

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Allocate a contiguous array with specified shape.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine AllocArray(array, shape)
  class(tArrayR), intent(inout) :: array
  integer(ip), intent(in) :: shape(:)

  allocate(array%mShape, mold=shape); array%mShape = shape
  allocate(array%mData(product(array%mShape)))

end subroutine AllocArray

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Allocate a contiguous array a mold.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do N = 1, NUM_RANKS
subroutine AllocArrayMold$N(@{array$$}@, mold)
  class(tArrayR), intent(inout) :: @{array$$}@
  class(tArrayR), intent(in) :: mold

#$do I = 1, N
  allocate(array$I%mShape, source=mold%mShape) 
  allocate(array$I%mData, mold=mold%mData)
#$end do

end subroutine AllocArrayMold$N
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Free an array. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do N = 1, NUM_RANKS
subroutine FreeArray$N(@{array$$}@)
  class(tArrayR), intent(inout) :: @{array$$}@

#$do I = 1, N
  deallocate(array$I%mShape)
  deallocate(array$I%mData)
#$end do

end subroutine FreeArray$N
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Get array rank.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
pure integer(ip) function ArrayRank(array)
  class(tArrayR), intent(in) :: array

  ArrayRank = size(array%mShape)

end function ArrayRank

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
  class(tArrayR), intent(in) :: array
  real(dp), intent(out), pointer :: fData(@:)

  integer(ip) :: rank
  integer(ip) :: fShape($fRank)

  rank = size(array%mShape)

  fShape = Rerank(array%mShape, $fRank)
  call c_f_pointer(cptr=c_loc(array%mData), fptr=fData, shape=fShape)

end subroutine GetArrayData$fRank
#$end do

end module StormRuler_Array
