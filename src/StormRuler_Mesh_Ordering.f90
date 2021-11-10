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
module StormRuler_Mesh_Ordering

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Helpers, only: I2S

use StormRuler_Mesh, only: tMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute mesh ordering quality, lower is better.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Mesh_Ordering_Quality(mesh, iperm) result(quality)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in), optional :: iperm(:)

  quality = mesh%RunCellKernel_Sum(Quality_Kernel)

contains
  real(dp) function Quality_Kernel(iCell) result(q)
    integer(ip), intent(in) :: iCell

    integer(ip) :: iCellFace

    ! ----------------------
    ! For each cell face:
    ! ----------------------
    q = 0.0_dp
    do iCellFace = 1, mesh%NumCellFaces
      ! ----------------------
      ! Quality is measured as the absolute value of the
      ! difference of the current index and adjacent cell index.
      ! (or permuted indices).
      ! ----------------------
      associate(jCell => mesh%CellToCell(iCellFace, iCell))
        if (jCell <= mesh%NumCells) then
          if (present(iperm)) then
            q = q + abs(iperm(iCell) - iperm(jCell))
          else
            q = q + abs(iCell - jCell)
          end if
        end if
      end associate
    end do

  end function Quality_Kernel
end function Mesh_Ordering_Quality

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Dump mesh ordering.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_Dump(mesh, file, iperm)
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: file
  integer(ip), intent(in), optional :: iperm(:)

  integer(ip) :: unit
  integer(ip) :: iCell

  open(newunit=unit, file=file, status='replace')

  do iCell = 1, mesh%NumCells
    if (present(iperm)) then
      write(unit, *) I2S(mesh%CellMDIndex(1,iperm(iCell)))//' '//I2S(mesh%CellMDIndex(2,iperm(iCell)))
    else
      write(unit, *) I2S(mesh%CellMDIndex(1,iCell))//' '//I2S(mesh%CellMDIndex(2,iCell))
    end if
  end do

  close(unit)

end subroutine Mesh_Ordering_Dump

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Swap function.
!! ----------------------------------------------------------------- !!
recursive subroutine Swap(i, j)
  integer(ip), intent(inout) :: i, j

  integer(ip) :: k

  k = i; i = j; j = k

end subroutine Swap

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a 2D mesh ordering along the Hilbert curve. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_Hilbert2D(mesh, iperm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(inout) :: iperm(:)

  integer(ip) :: iCell

  ! ----------------------
  ! Generate identity inverse permutation.
  ! ----------------------
  do iCell = 1, mesh%NumCells
    iperm(iCell) = iCell
  end do

  call Sort(1, mesh%NumCells + 1, &
    & 0.0_dp*mesh%MDIndexBounds + 1.0_dp, &
    & 1.0_dp*mesh%MDIndexBounds, 2, 0)

contains
  recursive subroutine Separate(iCell, iCellEnd, centerCoords, dim, sign)
    integer(ip), intent(inout) :: iCell, iCellEnd
    real(dp), intent(in) :: centerCoords(:)
    integer(ip), intent(in) :: dim, sign

    do while(iCell /= iCellEnd)
      if (sign*(mesh%CellMDIndex(dim,iperm(iCell)) - centerCoords(dim)) > 0) then
        iCellEnd = iCellEnd - 1
        call Swap(iperm(iCell), iperm(iCellEnd))
      else
        iCell = iCell + 1
      end if
    end do
  end subroutine Separate
  recursive subroutine Sort(iCell, iCellEnd, &
      & lowerCoords,upperCoords, orientation, threads)
    integer(ip), intent(in) :: iCell, iCellEnd
    real(dp), intent(in) :: lowerCoords(:),upperCoords(:)
    integer(ip), intent(in) :: orientation, threads

    integer(ip) :: iCellPiv1, iCellEnd1
    integer(ip) :: iCellPiv2, iCellEnd2
    integer(ip) :: iCellPiv3, iCellEnd3
    real(dp) :: centerCoords(2)

    ! ----------------------
    ! Check if recursion terminates.
    ! ----------------------
    if (iCell >= iCellEnd - 1) return
    if (any(lowerCoords >= upperCoords)) return

    ! ----------------------
    ! Separate quoters based on the orientation.
    ! ----------------------
    centerCoords = 0.5_dp*(lowerCoords + upperCoords)
    select case(orientation)
      
      ! ----------------------
      ! 4-|3
      ! 1-|2 orienation.
      ! ----------------------
      case(1)
        ! Separate lower and upper halves.
        iCellPiv2 = iCell; iCellEnd2 = iCellEnd
        call Separate(iCellPiv2, iCellEnd2, centerCoords, 2, +1)

        ! Separate left and right quadrants of the lower half.
        iCellPiv1 = iCell; iCellEnd1 = iCellPiv2
        call Separate(iCellPiv1, iCellEnd1, centerCoords, 1, +1)

        ! Separate right and left quadrants of the upper half.
        iCellPiv3 = iCellPiv2; iCellEnd3 = iCellEnd
        call Separate(iCellPiv3, iCellEnd3, centerCoords, 1, -1)

        ! Recursively process the quadrants.
        call Sort(iCell, iCellPiv1, &
          & [ lowerCoords(1),  lowerCoords(2)], &
          & [centerCoords(1), centerCoords(2)], &
          & 2, 0)
        call Sort(iCellPiv1, iCellPiv2, &
          & [centerCoords(1),  lowerCoords(2)], &
          & [ upperCoords(1), centerCoords(2)], &
          & 1, 0)
        call Sort(iCellPiv2, iCellPiv3, &
          & [centerCoords(1), centerCoords(2)], &
          & [ upperCoords(1),  upperCoords(2)], &
          & 1, 0)
        call Sort(iCellPiv3, iCellEnd, &
          & [ lowerCoords(1), centerCoords(2)], &
          & [centerCoords(1),  upperCoords(2)], &
          & 3, 0)

      ! ----------------------
      ! 2--3
      ! 1||4 orientation.
      ! ----------------------
      case(2)
        ! Separate left and right halves.
        iCellPiv2 = iCell; iCellEnd2 = iCellEnd
        call Separate(iCellPiv2, iCellEnd2, centerCoords, 1, +1)

        ! Separate lower and upper quadrants of the left half.
        iCellPiv1 = iCell; iCellEnd1 = iCellPiv2
        call Separate(iCellPiv1, iCellEnd1, centerCoords, 2, +1)

        ! Separate upper and lower quadrants of the right half.
        iCellPiv3 = iCellPiv2; iCellEnd3 = iCellEnd
        call Separate(iCellPiv3, iCellEnd3, centerCoords, 2, -1)

        ! Recursively process the quadrants.
        call Sort(iCell, iCellPiv1, &
          & [ lowerCoords(1),  lowerCoords(2)], &
          & [centerCoords(1), centerCoords(2)], &
          & 1, 0)
        call Sort(iCellPiv1, iCellPiv2, &
          & [ lowerCoords(1), centerCoords(2)], &
          & [centerCoords(1),  upperCoords(2)], &
          & 2, 0)
        call Sort(iCellPiv2, iCellPiv3, &
          & [centerCoords(1), centerCoords(2)], &
          & [ upperCoords(1),  upperCoords(2)], &
          & 2, 0)
        call Sort(iCellPiv3, iCellEnd, &
          & [centerCoords(1),  lowerCoords(2)], &
          & [ upperCoords(1), centerCoords(2)], &
          & 4, 0)

      ! ----------------------
      ! 4||1
      ! 3--2 orientation.
      ! ----------------------
      case(3)
        ! Separate right and left halves.
        iCellPiv2 = iCell; iCellEnd2 = iCellEnd
        call Separate(iCellPiv2, iCellEnd2, centerCoords, 1, -1)

        ! Separate upper and lower quadrants of the right half.
        iCellPiv1 = iCell; iCellEnd1 = iCellPiv2
        call Separate(iCellPiv1, iCellEnd1, centerCoords, 2, -1)

        ! Separate lower and upper quadrants of the left half.
        iCellPiv3 = iCellPiv2; iCellEnd3 = iCellEnd
        call Separate(iCellPiv3, iCellEnd3, centerCoords, 2, +1)

        ! Recursively process the quadrants.
        call Sort(iCell, iCellPiv1, &
          & [centerCoords(1), centerCoords(2)], &
          & [ upperCoords(1),  upperCoords(2)], &
          & 4, 0)
        call Sort(iCellPiv1, iCellPiv2, &
          & [centerCoords(1) , lowerCoords(2)], &
          & [ upperCoords(1), centerCoords(2)], &
          & 3, 0)
        call Sort(iCellPiv2, iCellPiv3, &
          & [ lowerCoords(1),  lowerCoords(2)], &
          & [centerCoords(1), centerCoords(2)], &
          & 3, 0)
        call Sort(iCellPiv3, iCellEnd, &
          & [ lowerCoords(1), centerCoords(2)], &
          & [centerCoords(1),  upperCoords(2)], &
          & 1, 0)

      ! ----------------------
      ! 2|-1
      ! 3|-4 orientation.
      ! ----------------------
      case(4)
        ! Separate upper and lower halves.
        iCellPiv2 = iCell; iCellEnd2 = iCellEnd
        call Separate(iCellPiv2, iCellEnd2, centerCoords, 2, -1)

        ! Separate right and left quadrants of the lower half.
        iCellPiv1 = iCell; iCellEnd1 = iCellPiv2
        call Separate(iCellPiv1, iCellEnd1, centerCoords, 1, -1)

        ! Separate left and right quadrants of the upper half.
        iCellPiv3 = iCellPiv2; iCellEnd3 = iCellEnd
        call Separate(iCellPiv3, iCellEnd3, centerCoords, 1, +1)

        ! Recursively process the quadrants.
        call Sort(iCell, iCellPiv1, &
          & [centerCoords(1), centerCoords(2)], &
          & [ upperCoords(1),  upperCoords(2)], &
          & 3, 0)
        call Sort(iCellPiv1, iCellPiv2, &
          & [ lowerCoords(1), centerCoords(2)], &
          & [centerCoords(1),  upperCoords(2)], &
          & 4, 0)
        call Sort(iCellPiv2, iCellPiv3, &
          & [ lowerCoords(1),  lowerCoords(2)], &
          & [centerCoords(1), centerCoords(2)], &
          & 4, 0)
        call Sort(iCellPiv3, iCellEnd, &
          & [centerCoords(1),  lowerCoords(2)], &
          & [ upperCoords(1), centerCoords(2)], &
          & 2, 0)

    end select

  end subroutine Sort
end subroutine Mesh_Ordering_Hilbert2D

end module StormRuler_Mesh_Ordering
