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

use StormRuler_Consts, only: ip, dp

use StormRuler_Helpers, only: Swap, I2S

use StormRuler_Mesh, only: tMesh

use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_null_ptr

#$use 'StormRuler_Macros.fi'

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

integer(ip) pure function ipermAt(cell, iperm)
  integer(ip), intent(in) :: cell
  integer(ip), intent(in), optional :: iperm(:)

  if (present(iperm)) then
    ipermAt = iperm(cell)
  else
    ipermAt = cell
  end if

end function ipermAt

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute mesh ordering quality, lower is better.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Mesh_Ordering_Quality(mesh, iperm) result(quality)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in), optional :: iperm(:)

  quality = mesh%RunCellKernel_Sum(Quality_Kernel)

contains
  real(dp) function Quality_Kernel(cell) result(q)
    integer(ip), intent(in) :: cell

    integer(ip) :: cellFace

    q = 0.0_dp

    do cellFace = 1, mesh%NumCellFaces
      associate(cellCell => mesh%CellToCell(cellFace, cell))

        if (cellCell <= mesh%NumCells) then
          q = q + log(abs(ipermAt(cell, iperm) - ipermAt(cellCell, iperm)) + 0.0_dp)
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
  integer(ip) :: cell

  open(newunit=unit, file=file, status='replace')

  do cell = 1, mesh%NumCells
    write(unit, *) &
      & I2S(mesh%CellToIndex(:,ipermAt(cell, iperm)))
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

  integer(ip) :: cell

  !$omp parallel do default(none) shared(mesh, iperm)
  do cell = 1, mesh%NumCells
    iperm(cell) = cell
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
subroutine Mesh_Ordering_HilbertCurve(mesh, iperm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(inout) :: iperm(:)

  integer(ip), allocatable :: ipermTemp(:)
  logical, allocatable :: matchTemp(:)

  allocate(ipermTemp, mold=iperm)
  allocate(matchTemp, mold=(iperm == iperm))

  ! ----------------------
  ! Generate identity inverse permutation.
  ! ----------------------
  call Mesh_Ordering_Identity(mesh, iperm)

  ! ----------------------
  ! Generate the ordering!
  ! ----------------------
  if (2 <= mesh%NumDims.or.mesh%NumDims <= 3) then
    call HilbertSort(1, mesh%NumCells + 1, 1, -1, -1, -1)
  else
    error stop 'Hilbert ordering supports 2D/3D meshes only.'
  end if

contains
  recursive pure logical function Condition(cell, centerCoords, dim, sign)
    integer(ip), intent(in) :: cell
    real(dp), intent(in) :: centerCoords(:)
    integer(ip), intent(in) :: dim, sign

    Condition = sign*(mesh%CellToIndex(dim,iperm(cell)) - centerCoords(dim)) < 0

  end function Condition
  recursive integer(ip) function Partition(cell, cellEnd, centerCoords, dim, sign)
    integer(ip), intent(in) :: cell, cellEnd
    real(dp), intent(in) :: centerCoords(:)
    integer(ip), intent(in) :: dim, sign

    integer(ip) :: cellCell, cellCellCell

    ! ----------------------
    ! Precompute the condition.
    ! ----------------------
    do cellCell = cell, cellEnd - 1
      matchTemp(cellCell) = Condition(cellCell, centerCoords, dim, sign)
    end do

    cellCellCell = cell

    ! ----------------------
    ! Copy the matching cells first.
    ! ----------------------
    do cellCell = cell, cellEnd - 1
      if (matchTemp(cellCell)) then
        ipermTemp(cellCellCell) = iperm(cellCell); cellCellCell = cellCellCell + 1
      end if
    end do

    Partition = cellCellCell

    ! ----------------------
    ! Copy the non-matching cells next.
    ! ----------------------
    do cellCell = cell, cellEnd - 1
      if (.not.matchTemp(cellCell)) then
        ipermTemp(cellCellCell) = iperm(cellCell); cellCellCell = cellCellCell + 1
      end if
    end do

    iperm(cell:(cellEnd - 1)) = ipermTemp(cell:(cellEnd - 1))

  end function Partition
  recursive subroutine HilbertSort(cell, cellEnd, dim1, sign1, sign2, sign3)
    integer(ip), intent(in) :: cell, cellEnd
    integer(ip), intent(in) :: dim1, sign1, sign2
    integer(ip), intent(in), optional :: sign3

    real(dp) :: centerCoords(mesh%NumDims)
    integer(ip) :: lowerCoords(mesh%NumDims), upperCoords(mesh%NumDims)
    integer(ip) :: cellPiv0
    integer(ip) :: cellPiv1, cellPiv2, cellPiv3, cellPiv4
    integer(ip) :: cellPiv5, cellPiv6, cellPiv7, cellPiv8
    integer(ip) :: dim2, dim3

    ! ----------------------
    ! Check if recursion terminates.
    ! ----------------------
    if (cell >= cellEnd - 1) return

    ! ----------------------
    ! Compute the bounding coordinates.
    ! ----------------------
    lowerCoords = mesh%CellToIndex(:,iperm(cell))
    upperCoords = mesh%CellToIndex(:,iperm(cell))
    do cellPiv1 = cell + 1, cellEnd - 1
      lowerCoords = min(lowerCoords, mesh%CellToIndex(:,iperm(cellPiv1)))
      upperCoords = max(upperCoords, mesh%CellToIndex(:,iperm(cellPiv1)))
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

      cellPiv0 = cell; cellPiv4 = cellEnd
      dim2 = mod(dim1, 2) + 1

      ! ----------------------
      ! Partition block into the quadrants.
      ! ----------------------
      cellPiv2 = Partition(cellPiv0, cellPiv4, centerCoords, dim1, +sign1)
      cellPiv1 = Partition(cellPiv0, cellPiv2, centerCoords, dim2, +sign2)
      cellPiv3 = Partition(cellPiv2, cellPiv4, centerCoords, dim2, -sign2)

      ! ----------------------
      ! Recursively sort the quadrants.
      ! ----------------------
      if (cellPiv1 /= cellPiv4) then
        call HilbertSort(cellPiv0, cellPiv1, dim2, +sign2, +sign1)
      end if
      if (cellPiv1 /= cellPiv0.or.cellPiv2 /= cellPiv4) then
        call HilbertSort(cellPiv1, cellPiv2, dim1, +sign1, +sign2)
      end if
      if (cellPiv2 /= cellPiv0.or.cellPiv3 /= cellPiv4) then
        call HilbertSort(cellPiv2, cellPiv3, dim1, +sign1, +sign2)
      end if
      if (cellPiv3 /= cellPiv0) then
        call HilbertSort(cellPiv3, cellPiv4, dim2, -sign2, -sign1)
      end if

    else if (mesh%NumDims == 3) then

      cellPiv0 = cell; cellPiv8 = cellEnd
      dim2 = mod(dim1, 3) + 1; dim3 = mod(dim1 + 1, 3) + 1

      ! ----------------------
      ! Partition block into the octants.
      ! ----------------------
      cellPiv4 = Partition(cellPiv0, cellPiv8, centerCoords, dim1, +sign1)
      cellPiv2 = Partition(cellPiv0, cellPiv4, centerCoords, dim2, +sign2)
      cellPiv6 = Partition(cellPiv4, cellPiv8, centerCoords, dim2, -sign2)
      cellPiv1 = Partition(cellPiv0, cellPiv2, centerCoords, dim3, +sign3)
      cellPiv3 = Partition(cellPiv2, cellPiv4, centerCoords, dim3, -sign3)
      cellPiv5 = Partition(cellPiv4, cellPiv6, centerCoords, dim3, +sign3)
      cellPiv7 = Partition(cellPiv6, cellPiv8, centerCoords, dim3, -sign3)

      ! ----------------------
      ! Recursively sort the octants.
      ! ----------------------
      if (cellPiv1 /= cellPiv8) then
        call HilbertSort(cellPiv0, cellPiv1, dim3, +sign3, +sign1, +sign2)
      end if
      if (cellPiv1 /= cellPiv0.or.cellPiv2 /= cellPiv8) then
        call HilbertSort(cellPiv1, cellPiv2, dim2, +sign2, +sign3, +sign1)
      end if
      if (cellPiv2 /= cellPiv0.or.cellPiv3 /= cellPiv8) then
        call HilbertSort(cellPiv2, cellPiv3, dim2, +sign2, +sign3, +sign1)
      end if
      if (cellPiv3 /= cellPiv0.or.cellPiv4 /= cellPiv8) then
        call HilbertSort(cellPiv3, cellPiv4, dim1, +sign1, -sign2, -sign3)
      end if
      if (cellPiv4 /= cellPiv0.or.cellPiv5 /= cellPiv8) then
        call HilbertSort(cellPiv4, cellPiv5, dim1, +sign1, -sign2, -sign3)
      end if
      if (cellPiv5 /= cellPiv0.or.cellPiv6 /= cellPiv8) then
        call HilbertSort(cellPiv5, cellPiv6, dim2, -sign2, +sign3, -sign1)
      end if
      if (cellPiv6 /= cellPiv0.or.cellPiv7 /= cellPiv8) then
        call HilbertSort(cellPiv6, cellPiv7, dim2, -sign2, +sign3, -sign1)
      end if
      if (cellPiv7 /= cellPiv0) then
        call HilbertSort(cellPiv7, cellPiv8, dim3, -sign3, -sign1, +sign2)
      end if

    end if

  end subroutine HilbertSort
end subroutine Mesh_Ordering_HilbertCurve

#$if HAS_METIS
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a 2D/3D mesh ordering with METIS.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh_Ordering_METIS(mesh, iperm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(inout) :: iperm(:)

  interface
    function cMETIS_NodeND(nvtxs, xadj, adjncy, &
        & vwgt, options, perm, iperm) bind(C, name='METIS_NodeND')
      import :: c_int, c_int32_t
      integer(c_int32_t), intent(in) :: nvtxs, xadj(*), adjncy(*)
      integer(c_int32_t), intent(in), optional :: vwgt(*), options(*)
      integer(c_int32_t), intent(inout) :: perm(*), iperm(*)
      integer(c_int) :: cMETIS_NodeND
    end function cMETIS_NodeND
  end interface

  integer(c_int) :: ierror
  integer(ip) :: cell, cellFace
  integer(ip), allocatable :: perm(:), xadj(:), adjncy(:)

  allocate(perm, mold=iperm)

  xadj = [0]
  do cell = 1, mesh%NumCells
    do cellFace = 1, mesh%NumCellFaces
      associate(cellCell => mesh%CellToCell(cellFace, cell))
        if (cellCell <= mesh%NumCells) adjncy = [adjncy, cellCell - 1]
      end associate
    end do
    xadj = [xadj, size(adjncy)]
  end do 

  print *, 'METIS'
  ierror = cMETIS_NodeND(mesh%NumCells, xadj, adjncy, perm=perm, iperm=iperm)

  perm(:) = perm(:) + 1
  iperm(:) = iperm(:) + 1
  print *, 'ierror=', ierror

end subroutine Mesh_Ordering_METIS
#$end if

end module StormRuler_Mesh_Ordering
