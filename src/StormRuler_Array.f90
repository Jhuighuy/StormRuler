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
  procedure :: Alloc => AllocArray
  procedure :: Free => FreeArray

  procedure :: Rank => ArrayRank

  procedure :: Get => GetArrayData
end type tArrayR

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
!! Free an array. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine FreeArray(array)
  class(tArrayR), intent(inout) :: array

  deallocate(array%mShape)
  deallocate(array%mData)

end subroutine FreeArray

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
subroutine GetArrayData(array, fData)
  class(tArrayR), intent(in) :: array
  real(dp), intent(out), pointer :: fData(..)

  integer(ip) :: rank

  rank = size(array%mShape)

  select rank(fData)
    rank(0)
      error stop 'array argument `data` expected.'
#$do fRank = 1, 15
    rank($fRank)
      call c_f_pointer(cptr=c_loc(array%mData), fptr=fData, shape=Rerank(array%mShape, $fRank))
#$end do
  end select

end subroutine GetArrayData

end module StormRuler_Array
