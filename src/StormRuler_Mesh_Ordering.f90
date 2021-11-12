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

use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_null_ptr

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

integer(ip) pure function iPermAt(iCell, iPerm)
  integer(ip), intent(in) :: iCell
  integer(ip), intent(in), optional :: iPerm(:)

  if (present(iPerm)) then
    iPermAt = iPerm(iCell)
  else
    iPermAt = iCell
  end if

end function iPermAt

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute mesh ordering quality, lower is better.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Mesh_Ordering_Quality(mesh, iPerm) result(quality)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in), optional :: iPerm(:)

  quality = mesh%RunCellKernel_Sum(Quality_Kernel)

contains
  real(dp) function Quality_Kernel(iCell) result(q)
    integer(ip), intent(in) :: iCell

    integer(ip) :: iCellFace

    q = 0.0_dp

    do iCellFace = 1, mesh%NumCellFaces
      associate(jCell => mesh%CellToCell(iCellFace, iCell))

        if (jCell <= mesh%NumCells) then
          q = q + log(abs(iPermAt(iCell, iPerm) - iPermAt(jCell, iPerm)) + 0.0_dp)
        end if

      end associate
    end do

  end function Quality_Kernel
end function Mesh_Ordering_Quality

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Dump mesh ordering.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_Dump(mesh, file, iPerm)
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: file
  integer(ip), intent(in), optional :: iPerm(:)

  integer(ip) :: unit
  integer(ip) :: iCell

  open(newunit=unit, file=file, status='replace')

  do iCell = 1, mesh%NumCells
    write(unit, *) &
      & I2S(mesh%CellMDIndex(1,iPermAt(iCell, iPerm)))//' '// &
      & I2S(mesh%CellMDIndex(2,iPermAt(iCell, iPerm)))
  end do

  close(unit)

end subroutine Mesh_Ordering_Dump

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute an identity ordering.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_Identity(mesh, iPerm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(inout) :: iPerm(:)

  integer(ip) :: iCell

  !$omp parallel do default(none) shared(mesh, iPerm)
  do iCell = 1, mesh%NumCells
    iPerm(iCell) = iCell
  end do
  !$omp end parallel do

end subroutine Mesh_Ordering_Identity

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a 2D/3D mesh ordering along the Hilbert curve. 
!!
!! References:
!! [1] MFEM source code: 
!!     https://github.com/mfem/mfem/blob/master/mesh/mesh.cpp#L1807
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_HilbertCurve(mesh, iPerm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(inout) :: iPerm(:)

  integer(ip), allocatable :: iPermTemp(:)
  logical, allocatable :: matchTemp(:)

  allocate(iPermTemp, mold=iPerm)
  allocate(matchTemp, mold=(iPerm == iPerm))

  ! ----------------------
  ! Generate identity inverse permutation.
  ! ----------------------
  call Mesh_Ordering_Identity(mesh, iPerm)

  ! ----------------------
  ! Generate the ordering!
  ! ----------------------
  if (2 <= mesh%NumDims.or.mesh%NumDims <= 3) then
    call HilbertSort(1, mesh%NumCells + 1, 1, -1, -1, -1)
  else
    error stop 'Hilbert ordering supports 2D/3D meshes only.'
  end if

contains
  recursive pure logical function Condition(iCell, centerCoords, dim, sign)
    integer(ip), intent(in) :: iCell
    real(dp), intent(in) :: centerCoords(:)
    integer(ip), intent(in) :: dim, sign

    Condition = sign*(mesh%CellMDIndex(dim,iPerm(iCell)) - centerCoords(dim)) < 0

  end function Condition
  recursive integer(ip) function Partition(iCell, iCellEnd, centerCoords, dim, sign)
    integer(ip), intent(in) :: iCell, iCellEnd
    real(dp), intent(in) :: centerCoords(:)
    integer(ip), intent(in) :: dim, sign

    integer(ip) :: jCell, kCell

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
        iPermTemp(kCell) = iPerm(jCell); kCell = kCell + 1
      end if
    end do

    Partition = kCell

    ! ----------------------
    ! Copy the non-matching cells next.
    ! ----------------------
    do jCell = iCell, iCellEnd - 1
      if (.not.matchTemp(jCell)) then
        iPermTemp(kCell) = iPerm(jCell); kCell = kCell + 1
      end if
    end do

    iPerm(iCell:(iCellEnd - 1)) = iPermTemp(iCell:(iCellEnd - 1))

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
    integer(ip) :: dim2, dim3

    ! ----------------------
    ! Check if recursion terminates.
    ! ----------------------
    if (iCell >= iCellEnd - 1) return

    ! ----------------------
    ! Compute the bounding coordinates.
    ! ----------------------
    lowerCoords = mesh%CellMDIndex(:,iPerm(iCell))
    upperCoords = mesh%CellMDIndex(:,iPerm(iCell))
    do iCellPiv1 = iCell + 1, iCellEnd - 1
      lowerCoords = min(lowerCoords, mesh%CellMDIndex(:,iPerm(iCellPiv1)))
      upperCoords = max(upperCoords, mesh%CellMDIndex(:,iPerm(iCellPiv1)))
    end do
    centerCoords = 0.5_dp*(lowerCoords + upperCoords)

    ! ----------------------
    ! Terminate if odd block is detected. 
    ! ----------------------
    associate(blockShape => upperCoords - lowerCoords + 1)
    
      if (any(mod(blockShape, 2) == 1)) return
    
    end associate

    ! ----------------------
    ! Continue the recursive partition.
    ! ----------------------
    if (mesh%NumDims == 2) then

      iCellPiv0 = iCell; iCellPiv4 = iCellEnd
      dim2 = mod(dim1, 2) + 1

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

    else if (mesh%NumDims == 3) then

      iCellPiv0 = iCell; iCellPiv8 = iCellEnd
      dim2 = mod(dim1, 3) + 1; dim3 = mod(dim1 + 1, 3) + 1

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

    end if

  end subroutine HilbertSort
end subroutine Mesh_Ordering_HilbertCurve

#$if HAS_METIS
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a 2D/3D mesh ordering with METIS.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_METIS(mesh, iPerm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(inout) :: iPerm(:)

  interface
    function cMETIS_NodeND(nvtxs, xadj, adjncy, &
        & vwgt, options, perm, iPerm) bind(C, name='METIS_NodeND')
      import :: c_int, c_int32_t
      integer(c_int32_t), intent(in) :: nvtxs, xadj(*), adjncy(*)
      integer(c_int32_t), intent(in), optional :: vwgt(*), options(*)
      integer(c_int32_t), intent(inout) :: perm(*), iPerm(*)
      integer(c_int) :: cMETIS_NodeND
    end function cMETIS_NodeND
  end interface

  integer(c_int) :: ierror
  integer(ip) :: iCell, iCellFace
  integer(ip), allocatable :: perm(:), xadj(:), adjncy(:)

  allocate(perm, mold=iPerm)

  xadj = [0]
  do iCell = 1, mesh%NumCells
    do iCellFace = 1, mesh%NumCellFaces
      associate(jCell => mesh%CellToCell(iCellFace, iCell))
        if (jCell <= mesh%NumCells) adjncy = [adjncy, jCell - 1]
      end associate
    end do
    xadj = [xadj, size(adjncy)]
  end do 

  print *, 'METIS'
  ierror = cMETIS_NodeND(mesh%NumCells, xadj, adjncy, perm=perm, iPerm=iPerm)

  perm(:) = perm(:) + 1
  iPerm(:) = iPerm(:) + 1
  print *, 'ierror=', ierror

end subroutine Mesh_Ordering_METIS
#$end if

end module StormRuler_Mesh_Ordering
