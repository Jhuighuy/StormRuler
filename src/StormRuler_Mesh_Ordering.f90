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
use StormRuler_Helpers, only: Swap, I2S

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

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute an identity ordering.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_Identity(mesh, iperm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(inout) :: iperm(:)

  integer(ip) :: iCell

  do iCell = 1, mesh%NumCells
    iperm(iCell) = iCell
  end do

end subroutine Mesh_Ordering_Identity

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a 2D/3D mesh ordering along the Hilbert curve. 
!!
!! References:
!! [1] MFEM source code: 
!!     https://github.com/mfem/mfem/blob/master/mesh/mesh.cpp#L1807
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_HilbertCurve(mesh, iperm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(inout) :: iperm(:)

  ! ----------------------
  ! Generate identity inverse permutation.
  ! ----------------------
  call Mesh_Ordering_Identity(mesh, iperm)

  ! ----------------------
  ! Generate the ordering!
  ! ----------------------
  if (mesh%NumDims == 2) then
    call HilbertSort(1, mesh%NumCells + 1, 1, -1, -1)
  else if (mesh%NumDims == 3) then
    call HilbertSort(1, mesh%NumCells + 1, 1, -1, -1, -1)
  else
    error stop 'Hilbert ordering supports 2D/3D meshes only.'
  end if

contains
  recursive pure logical function Condition(iCell, centerCoords, dim, sign)
    integer(ip), intent(in) :: iCell
    real(dp), intent(in) :: centerCoords(:)
    integer(ip), intent(in) :: dim, sign

    Condition = sign*(mesh%CellMDIndex(dim,iperm(iCell)) - centerCoords(dim)) < 0

  end function Condition
  recursive integer(ip) function Partition(iCell, iCellEnd, centerCoords, dim, sign)
    integer(ip), intent(in) :: iCell, iCellEnd
    real(dp), intent(in) :: centerCoords(:)
    integer(ip), intent(in) :: dim, sign

    integer(ip) :: jCell, kCell

    integer(ip), allocatable :: ipermTemp(:)
    logical, allocatable :: matchTemp(:)

    allocate(ipermTemp(iCell:(iCellEnd - 1)))
    allocate(matchTemp(iCell:(iCellEnd - 1)))

    ! ----------------------
    ! Precompute the condition.
    ! ----------------------
    do jCell = iCell, iCellEnd - 1
      matchTemp(jCell) = Condition(jCell, centerCoords, dim, sign)
    end do

    kCell = iCell

    ! ----------------------
    ! Copy the matching cells first.
    ! ----------------------
    do jCell = iCell, iCellEnd - 1
      if (matchTemp(jCell)) then
        ipermTemp(kCell) = iperm(jCell); kCell = kCell + 1
      end if
    end do

    Partition = kCell

    ! ----------------------
    ! Copy the non-matching cells next.
    ! ----------------------
    do jCell = iCell, iCellEnd - 1
      if (.not.matchTemp(jCell)) then
        ipermTemp(kCell) = iperm(jCell); kCell = kCell + 1
      end if
    end do

    iperm(iCell:(iCellEnd - 1)) = ipermTemp(iCell:(iCellEnd - 1))

  end function Partition
  recursive subroutine HilbertSort(iCell, iCellEnd, dim1, sign1, sign2, sign3)
    integer(ip), intent(in) :: iCell, iCellEnd
    integer(ip), intent(in) :: dim1, sign1, sign2
    integer(ip), intent(in), optional :: sign3

    real(dp) :: centerCoords(mesh%NumDims)
    integer(ip) :: lowerCoords(mesh%NumDims), upperCoords(mesh%NumDims)
    integer(ip) :: iCellPiv0
    integer(ip) :: iCellPiv1, iCellPiv2, iCellPiv3, iCellPiv4
    integer(ip) :: iCellPiv5, iCellPiv6, iCellPiv7, iCellPiv8

    ! ----------------------
    ! Check if recursion terminates.
    ! ----------------------
    if (iCell >= iCellEnd - 1) return

    ! ----------------------
    ! Compute the bounding coordinates.
    ! ----------------------
    lowerCoords = mesh%CellMDIndex(:,iperm(iCell))
    upperCoords = mesh%CellMDIndex(:,iperm(iCell))
    do iCellPiv1 = iCell + 1, iCellEnd - 1
      lowerCoords = min(lowerCoords, mesh%CellMDIndex(:,iperm(iCellPiv1)))
      upperCoords = max(upperCoords, mesh%CellMDIndex(:,iperm(iCellPiv1)))
    end do
    centerCoords = 0.5_dp*(lowerCoords + upperCoords)

    ! ----------------------
    ! Terminate if odd block is detected. 
    ! ----------------------
    associate(blockShape => upperCoords - lowerCoords + 1)
      if (any(mod(abs(blockShape), 2) == 1)) then
        iperm(iCell:(iCellEnd - 1)) = iperm((iCellEnd - 1):iCell:-1)
        return
      end if
    end associate

    ! ----------------------
    ! Continue the recursive partition.
    ! ----------------------
    if (mesh%NumDims == 2) then

      iCellPiv0 = iCell; iCellPiv4 = iCellEnd
      associate(dim2 => mod(dim1, 2) + 1)

        ! ----------------------
        ! Partition block into the quadrants.
        ! ----------------------
        iCellPiv2 = Partition(iCellPiv0, iCellPiv4, centerCoords, dim1, +sign1)
        iCellPiv1 = Partition(iCellPiv0, iCellPiv2, centerCoords, dim2, +sign2)
        iCellPiv3 = Partition(iCellPiv2, iCellPiv4, centerCoords, dim2, -sign2)

        ! ----------------------
        ! Recursively sort the quadrants.
        ! ----------------------
        if (iCellPiv1 /= iCellPiv4) then
          call HilbertSort(iCellPiv0, iCellPiv1, dim2, +sign2, +sign1)
        end if
        if (iCellPiv1 /= iCellPiv0.or.iCellPiv2 /= iCellPiv4) then
          call HilbertSort(iCellPiv1, iCellPiv2, dim1, +sign1, +sign2)
        end if
        if (iCellPiv2 /= iCellPiv0.or.iCellPiv3 /= iCellPiv4) then
          call HilbertSort(iCellPiv2, iCellPiv3, dim1, +sign1, +sign2)
        end if
        if (iCellPiv3 /= iCellPiv0) then
          call HilbertSort(iCellPiv3, iCellPiv4, dim2, -sign2, -sign1)
        end if

      end associate

    else if (mesh%NumDims == 3) then

      iCellPiv0 = iCell; iCellPiv8 = iCellEnd
      associate(dim2 => mod(dim1, 3) + 1, dim3 => mod(dim1 + 1, 3) + 1)

        ! ----------------------
        ! Partition block into the octants.
        ! ----------------------
        iCellPiv4 = Partition(iCellPiv0, iCellPiv8, centerCoords, dim1, +sign1)
        iCellPiv2 = Partition(iCellPiv0, iCellPiv4, centerCoords, dim2, +sign2)
        iCellPiv6 = Partition(iCellPiv4, iCellPiv8, centerCoords, dim2, -sign2)
        iCellPiv1 = Partition(iCellPiv0, iCellPiv2, centerCoords, dim3, +sign3)
        iCellPiv3 = Partition(iCellPiv2, iCellPiv4, centerCoords, dim3, -sign3)
        iCellPiv5 = Partition(iCellPiv4, iCellPiv6, centerCoords, dim3, +sign3)
        iCellPiv7 = Partition(iCellPiv6, iCellPiv8, centerCoords, dim3, -sign3)

        ! ----------------------
        ! Recursively sort the octants.
        ! ----------------------
        if (iCellPiv1 /= iCellPiv8) then
          call HilbertSort(iCellPiv0, iCellPiv1, dim3, +sign3, +sign1, +sign2)
        end if
        if (iCellPiv1 /= iCellPiv0.or.iCellPiv2 /= iCellPiv8) then
          call HilbertSort(iCellPiv1, iCellPiv2, dim2, +sign2, +sign3, +sign1)
        end if
        if (iCellPiv2 /= iCellPiv0.or.iCellPiv3 /= iCellPiv8) then
          call HilbertSort(iCellPiv2, iCellPiv3, dim2, +sign2, +sign3, +sign1)
        end if
        if (iCellPiv3 /= iCellPiv0.or.iCellPiv4 /= iCellPiv8) then
          call HilbertSort(iCellPiv3, iCellPiv4, dim1, +sign1, -sign2, -sign3)
        end if
        if (iCellPiv4 /= iCellPiv0.or.iCellPiv5 /= iCellPiv8) then
          call HilbertSort(iCellPiv4, iCellPiv5, dim1, +sign1, -sign2, -sign3)
        end if
        if (iCellPiv5 /= iCellPiv0.or.iCellPiv6 /= iCellPiv8) then
          call HilbertSort(iCellPiv5, iCellPiv6, dim2, -sign2, +sign3, -sign1)
        end if
        if (iCellPiv6 /= iCellPiv0.or.iCellPiv7 /= iCellPiv8) then
          call HilbertSort(iCellPiv6, iCellPiv7, dim2, -sign2, +sign3, -sign1)
        end if
        if (iCellPiv7 /= iCellPiv0) then
          call HilbertSort(iCellPiv7, iCellPiv8, dim3, -sign3, -sign1, +sign2)
        end if

      end associate

    end if

  end subroutine HilbertSort
end subroutine Mesh_Ordering_HilbertCurve

end module StormRuler_Mesh_Ordering
