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
module StormRuler_Mesh

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Helpers, only: Flip, IndexOf, I2S, R2S, IntToPixel
use StormRuler_IO, only: IOList, IOListItem, @{IOListItem$$@|@0, 2}@

use, intrinsic :: iso_fortran_env, only: error_unit

#$if HAS_OpenMP
use :: omp_lib
#$endif

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Semi-structured multidimensional mesh.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
type :: tMesh
  ! ----------------------
  ! Number of spatial dimension.
  ! ----------------------
  integer(ip) :: NumDims

  ! ----------------------
  ! Number of faces (edges in 2D) per each cell.
  ! NumCellFaces = 4 for 2D meshes, NumCellFaces = 6 for 3D meshes.
  ! ----------------------
  integer(ip) :: NumCellFaces
  ! ----------------------
  ! Number of the connection per each cell.
  ! Typically, NumConns = NumCellFaces.
  ! For LBM simulations, extended connectivity is required, here 
  ! NumConns = 9 for 2D meshes, NumConns = 19 for 3D meshes (D2Q9 and D3Q19).
  ! ----------------------
  integer(ip) :: NumConns

  ! ----------------------
  ! Number of interior cells.
  ! This value should be used for field operations.
  ! ----------------------
  integer(ip) :: NumCells
  ! ----------------------
  ! Number of boundary cells 
  ! (domain-exterior cells, first in the ghost layers).
  ! Cells may be counted multiple times..
  ! ----------------------
  integer(ip) :: NumBCCells
  ! ----------------------
  ! Total number of cells (including interior and ghost cells).
  ! This value should be used for field allocation.
  ! ----------------------
  integer(ip) :: NumAllCells
  ! ----------------------
  ! Cell connectivity table.
  ! Shape is [1, NumConns]×[1, NumAllCells].
  ! ----------------------
  integer(ip), allocatable :: CellToCell(:,:)
  ! ----------------------
  ! Logical table for the periodic face connections.
  ! Shape is [1, NumConns]×[1, NumAllCells].
  ! ----------------------
  logical, allocatable, private :: mIsCellFacePeriodic(:,:)

  ! ----------------------
  ! Multidimensional index bounds table.
  ! Shape is [1, NumDims].
  ! ----------------------
  integer(ip), allocatable :: MDIndexBounds(:)
  ! ----------------------
  ! Cell multidimensional index table.
  ! Shape is [1, NumDims]×[1, NumAllCells].
  ! ----------------------
  integer(ip), allocatable :: CellMDIndex(:,:)

  ! ----------------------
  ! Number of boundary condition marks.
  ! ----------------------
  integer(ip) :: NumBCMs
  ! ----------------------
  ! BC mark to boundary cell index (in CSR format).
  ! Shape is [1, NumBCMs].
  ! ----------------------
  integer(ip), allocatable :: BCMs(:)
  ! ----------------------
  ! Shape is [1, NumBCCells].
  ! ----------------------
  integer(ip), allocatable :: BCMToCell(:)
  ! ----------------------
  ! Shape is [1, NumBCCells].
  ! ----------------------
  integer(ip), allocatable :: BCMToCellFace(:)

  ! ----------------------
  ! Distance between centers of the adjacent cells per face.
  ! Shape is [1, NumConns].
  ! ----------------------
  real(dp), allocatable :: dl(:)
  ! ----------------------
  ! Difference between centers of the adjacent cells per face.
  ! Shape is [1, NumDims]×[1, NumConns].
  ! ----------------------
  real(dp), allocatable :: dr(:,:)

contains

  ! ----------------------
  ! Initializers.
  ! ----------------------
  procedure :: InitRect => tMesh_InitRect
  procedure :: InitFromImage2D => tMesh_InitFromImage2D
  procedure :: InitFromImage3D => tMesh_InitFromImage3D
  procedure, private :: GenerateBCCells => tMesh_GenerateBCCells

  ! ----------------------
  ! Parallel mesh walkthough subroutines.
  ! ----------------------
  procedure :: RunCellKernel => RunMeshCellSimpleKernel
  procedure :: RunCellKernel_Block => RunMeshBlockCellKernel
  procedure :: RunCellKernel_Sum => RunMeshSumCellKernel
  procedure :: RunCellKernel_Min => RunMeshMinCellKernel
  procedure :: RunCellKernel_Max => RunMeshMaxCellKernel

  ! ----------------------
  ! Field wrappers.
  ! ----------------------
  generic :: CellCenter => CellCenterVec, CellCenterCoord
  procedure :: CellCenterVec => GetMeshCellCenterVec
  procedure :: CellCenterCoord => GetMeshCellCenterCoord
  procedure :: IsCellFacePeriodic => IsMeshCellFacePeriodic

  ! ----------------------
  ! Ordering routines. 
  ! ----------------------
  procedure :: ApplyOrdering => ApplyMeshOrdering
  procedure :: GenerateExtraConnectivity => GenerateMeshExtraConnectivity 

  ! ----------------------
  ! Printers.
  ! ----------------------
  procedure :: PrintTo_Neato => PrintMeshToNeato
  procedure :: PrintTo_LegacyVTK => PrintMeshToLegacyVTK
  
end type tMesh

abstract interface
  subroutine tKernelFunc(cell)
    import ip
    integer(ip), intent(in) :: cell
  end subroutine tKernelFunc
end interface

abstract interface
  subroutine tBlockKernelFunc(firstCell, lastCell)
    import ip
    integer(ip), intent(in) :: firstCell, lastCell
  end subroutine tBlockKernelFunc
end interface

abstract interface
  function tReduceKernelFunc(cell) result(r)
    import ip, dp
    integer(ip), intent(in) :: cell
    real(dp) :: r
  end function tReduceKernelFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Launch a cell kernel.
!! ----------------------------------------------------------------- !!
subroutine RunMeshCellSimpleKernel(mesh, Kernel)
  class(tMesh), intent(in) :: mesh
  procedure(tKernelFunc) :: Kernel

  integer :: cell

  !$omp parallel do default(none) shared(mesh) schedule(static)
  do cell = 1, mesh%NumCells
    call Kernel(cell)
  end do
  !$omp end parallel do

end subroutine RunMeshCellSimpleKernel

!! ----------------------------------------------------------------- !!
!! Launch a block-kernel.
!! ----------------------------------------------------------------- !!
subroutine RunMeshBlockCellKernel(mesh, BlockKernel)
  class(tMesh), intent(in) :: mesh
  procedure(tBlockKernelFunc) :: BlockKernel

#$if not HAS_OpenMP

  call BlockKernel(1, mesh%NumCells)

#$else

  integer(ip) :: i, thread
  integer(ip), allocatable :: ranges(:)

  ! ----------------------
  ! Compute block range for each thread.
  ! ----------------------
  allocate(ranges(omp_get_max_threads()))
  associate(size => mesh%NumCells)
    ranges(:) = size/omp_get_max_threads()
    associate(remainder => ranges(:mod(size, omp_get_max_threads())))
      remainder(:) = remainder(:) + 1
    end associate
  end associate
  ranges = [1, ranges]
  do i = 1, omp_get_max_threads()
    ranges(i + 1) = ranges(i) + ranges(i + 1) 
  end do

  ! ----------------------
  ! Lauch threads.
  ! ----------------------
  !$omp parallel default(none) shared(ranges) private(thread)
  thread = omp_get_thread_num() + 1
  call BlockKernel(ranges(thread), ranges(thread + 1) - 1)
  !$omp end parallel

#$endif
end subroutine RunMeshBlockCellKernel

!! ----------------------------------------------------------------- !!
!! Launch a SUM-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function RunMeshSumCellKernel(mesh, Kernel) result(sum)
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: sum

  integer :: cell

  sum = 0.0_dp
  !$omp parallel do default(none) shared(mesh) &
  !$omp & schedule(static) reduction(+:sum)
  do cell = 1, mesh%NumCells
    sum = sum + Kernel(cell)
  end do
  !$omp end parallel do

end function RunMeshSumCellKernel

!! ----------------------------------------------------------------- !!
!! Launch a min-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function RunMeshMinCellKernel(mesh, Kernel) result(minValue)
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: minValue

  integer :: cell

  minValue = +huge(minValue)
  !$omp parallel do default(none) shared(mesh) &
  !$omp & schedule(static) reduction(min:minValue)
  do cell = 1, mesh%NumCells
    minValue = min(minValue, Kernel(cell))
  end do
  !$omp end parallel do

end function RunMeshMinCellKernel

!! ----------------------------------------------------------------- !!
!! Launch a max-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function RunMeshMaxCellKernel(mesh, Kernel) result(maxValue)
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: maxValue

  integer :: cell

  maxValue = -huge(maxValue)
  !$omp parallel do default(none) shared(mesh) &
  !$omp & schedule(static) reduction(max:maxValue)
  do cell = 1, mesh%NumCells
    maxValue = max(maxValue, Kernel(cell))
  end do
  !$omp end parallel do

end function RunMeshMaxCellKernel

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Get cell center.
!! ----------------------------------------------------------------- !!
pure function GetMeshCellCenterVec(mesh, cell) result(cellCenter)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: cell
  real(dp) :: cellCenter(mesh%NumDims)

  !! TODO: add initial offset as mesh class field.
  associate(cellIndex => mesh%CellMDIndex(:,cell))
    cellCenter = [0.0_dp,0.0_dp] + mesh%dl(:2*mesh%NumDims:2)*(cellIndex - 0.5_dp)
  end associate

end function GetMeshCellCenterVec
pure function GetMeshCellCenterCoord(mesh, dim, cell) result(cellCenterCoord)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: dim, cell
  real(dp) :: cellCenterCoord

  associate(cellCenter => mesh%CellCenterVec(cell))
    cellCenterCoord = cellCenter(dim)
  end associate

end function GetMeshCellCenterCoord

!! ----------------------------------------------------------------- !!
!! Check if cell face connects it to the BC cell or 
!! to the interioir periodic. 
!! ----------------------------------------------------------------- !!
logical pure function IsMeshCellFacePeriodic(mesh, cellFace, cell) result(isPeriodicFace)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: cell, cellFace
  
  isPeriodicFace = allocated(mesh%mIsCellFacePeriodic)
  if (isPeriodicFace) then
    isPeriodicFace = mesh%mIsCellFacePeriodic(cellFace, cell)
  end if

end function IsMeshCellFacePeriodic

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

#$let STENCIL = { &
  & 2:{ 1:[+1, 0], 2:[-1, 0],   &
  &     3:[0, +1], 4:[0, -1] }, &
  & 3:{ 1:[+1, 0, 0], 2:[-1, 0, 0], &
  &     3:[0, +1, 0], 4:[0, -1, 0], &
  &     5:[0, 0, +1], 6:[0, 0, -1] } }

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize a mesh from image.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do dim = 2, 3
subroutine tMesh_InitFromImage${dim}$D(mesh, image, fluidColor, colorToBCM, numBCLayers)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in) :: image(@{0:}@), fluidColor, colorToBCM(:), numBCLayers
  ! >>>>>>>>>>>>>>>>>>>>>>

#$define iCellMD @{iCellMD$$}@
  integer(ip) :: iCell, $iCellMD, iBCM
  integer(ip), allocatable :: cache(@:)
  allocate(cache, mold=image); cache(@:) = 0

  ! ----------------------
  ! Process image data in order to:
  ! 1. Generate initial "Cell MD Index" ⟺ "Cell (Plain Index)" table 
  !    for interior cells and compute number of interior and boundary cells.
  ! 2. Generate and pre-fill "BCM" ⟺ "Num BC Cells per BCM" CSR-style table.
  ! ----------------------
  mesh%NumDims = $dim
  mesh%NumCellFaces = ${2*dim}$
  mesh%NumConns = merge(9, 19, mesh%NumDims == 2)
  mesh%MDIndexBounds = shape(image)-2
  mesh%NumBCMs = size(colorToBCM, dim=1)
  allocate(mesh%BCMs(mesh%NumBCMs+1)); mesh%BCMs(:) = 0

  ! ----------------------
  iCell = 0; mesh%BCMs(:) = 0
#$do rank = dim, 1, -1
  do iCellMD$rank = 1, mesh%MDIndexBounds($rank) 
#$end do
    if (image($iCellMD) == fluidColor) then
      iCell = iCell + 1
      cache($iCellMD) = iCell
#$for iCellFace, iCellFaceMD in STENCIL[dim].items()
#$define iiCellMD @{iCellMD$$+${iCellFaceMD[$$-1]}$}@
      if (image($iiCellMD) /= fluidColor) then
        iBCM = IndexOf(image($iiCellMD), colorToBCM)
        if (iBCM == 0) then
          write(error_unit, *) 'INVALID IMAGE COLOR', &
            & IntToPixel(image($iiCellMD)), 'AT', $iiCellMD
          error stop
        end if
        mesh%BCMs(iBCM+1) = mesh%BCMs(iBCM+1) + 1
      end if
#$end for
    end if
#$do rank = dim, 1, -1
  end do
#$end do
  do iBCM = 2, mesh%NumBCMs
    mesh%BCMs(iBCM+1) = mesh%BCMs(iBCM+1) + mesh%BCMs(iBCM)
  end do

  ! ----------------------
  ! TODO: apply plain indices sorting.
  ! ----------------------
  mesh%NumCells = iCell
  mesh%NumBCCells = mesh%BCMs(mesh%NumBCMs+1)
  mesh%NumAllCells = mesh%NumCells+mesh%NumBCCells*numBCLayers
  print *, 'NC:', mesh%NumCells, 'NBC:', mesh%NumBCCells, 'NAC:', mesh%NumAllCells
  print *, 'BCMs:', mesh%BCMs

  ! ----------------------
  ! 4. Allocate the connectivity tables and fill them with data,
  !    inverting the "Cell MD Index" ⟺ "Cell (Plain Index)" table.
  ! ----------------------
  associate(nac => mesh%NumAllCells, ncc => mesh%NumConns)
    allocate(mesh%CellMDIndex(mesh%NumDims, nac))
    allocate(mesh%CellToCell(ncc, nac)); mesh%CellToCell(:,:) = 0
  end associate
  ! ----------------------
#$do rank = dim, 1, -1
  do iCellMD$rank = 1, mesh%MDIndexBounds($rank) 
#$end do
    if (image($iCellMD) == fluidColor) then
      iCell = cache($iCellMD)
      mesh%CellMDIndex(:,iCell) = [$iCellMD]
#$for iCellFace, iCellFaceMD in STENCIL[dim].items()
#$define iiCellMD @{iCellMD$$+${iCellFaceMD[$$-1]}$}@
      if (image($iiCellMD) == fluidColor) then
        mesh%CellToCell($iCellFace, iCell) = cache($iiCellMD)
      else
        iBCM = IndexOf(image($iiCellMD), colorToBCM)
        mesh%CellToCell($iCellFace, iCell) = -iBCM
      end if
#$end for
    end if
#$do rank = dim, 1, -1
  end do
#$end do
  
  ! ----------------------
  ! Clean-up.
#$del iCellMD, iiCellMD
  deallocate(cache)
  
  ! ----------------------
  ! Proceed to the next steps..
  ! ----------------------
  call mesh%GenerateBCCells(numBCLayers)
end subroutine tMesh_InitFromImage${dim}$D
#$end do

#$del STENCIL

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Generate BC/ghost cells.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine tMesh_GenerateBCCells(mesh, numBCLayers)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in) :: numBCLayers
  
  integer(ip) :: iCell, iCellFace, iBCM, iBCCell, iGCell
  
  ! ----------------------
  ! 5. Allocate and fill the "BCM" ⟺ "BC Cell" 
  !    and "BCM" ⟺ "BC Cell Face" CSR-style tables and fill it. 
  ! ----------------------
  associate(nbcc => mesh%NumBCCells)
    allocate(mesh%BCMToCell(nbcc))
    allocate(mesh%BCMToCellFace(nbcc))
  end associate
  ! ----------------------
  iBCCell = mesh%NumCells + 1
  do iCell = 1, mesh%NumCells
    do iCellFace = 1, mesh%NumCellFaces
      iBCM = -mesh%CellToCell(iCellFace, iCell) 
      if (iBCM > 0) then
        ! Classic CSR insertion procedure here.
        associate(ptr => mesh%BCMs(iBCM))
          mesh%BCMToCell(ptr+1) = iBCCell
          mesh%BCMToCellFace(ptr+1) = iCellFace; ptr = ptr + 1
        end associate
        ! Generate connectivity for BC cell and ghost cells.
        mesh%CellToCell(iCellFace, iCell) = iBCCell
        do iGCell = iBCCell, iBCCell+numBCLayers-1
          mesh%CellToCell(iCellFace, iGCell) = iGCell+1
          mesh%CellToCell(Flip(iCellFace), iGCell) = iGCell-1
          ! Fill MD index information by linear extrapolation.
          associate(iCellMD => mesh%CellMDIndex(:,iCell), &
            &      iiCellMD => mesh%CellMDIndex(:,mesh%CellToCell(Flip(iCellFace), iCell)))
            mesh%CellMDIndex(:,iGCell) = &
              & iCellMD + (iGCell-iBCCell+1)*(iCellMD-iiCellMD)
          end associate
        end do
        mesh%CellToCell(iCellFace, iBCCell+numBCLayers-1) = 0
        mesh%CellToCell(Flip(iCellFace), iBCCell) = iCell
        ! Increment BC cell index.
        iBCCell = iBCCell+numBCLayers
      end if
    end do
  end do
  
  ! ----------------------
  ! 6. Classic CSR pointer correction procedure.
  ! ----------------------
  mesh%BCMs = eoshift(mesh%BCMs, shift=-1) + 1
  print *, 'ERROR PRONE CHECK:', iBCCell-1
  print *, 'AGAIN BCMs:', mesh%BCMs
end subroutine tMesh_GenerateBCCells

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the specified cell ordering.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ApplyMeshOrdering(mesh, iperm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in) :: iperm(:)

  integer(ip) :: cell, cellFace
  integer(ip), allocatable :: cellToCellTemp(:,:)
  integer(ip), allocatable :: cellIndexTemp(:,:)

  ! ----------------------
  ! Order cell-cell graph columns.
  ! ----------------------
  do cell = 1, mesh%NumAllCells
    do cellFace = 1, mesh%NumCellFaces
      associate(cellCell => mesh%CellToCell(cellFace, cell))

        if (0 < cellCell.and.cellCell <= mesh%NumCells) cellCell = iperm(cellCell)

      end associate
    end do
  end do

  ! ----------------------
  ! Order cell-cell graph rows.
  ! ----------------------
  cellToCellTemp = mesh%CellToCell(:,:mesh%NumCells)
  cellIndexTemp = mesh%CellMDIndex(:,:mesh%NumCells)
  do cell = 1, mesh%NumCells
    associate(cellCell => iperm(cell))

      mesh%CellToCell(:,cellCell) = cellToCellTemp(:,cell)
      mesh%CellMDIndex(:,cellCell) = cellIndexTemp(:,cell)

    end associate
  end do

end subroutine ApplyMeshOrdering

subroutine GenerateMeshExtraConnectivity(mesh)
  class(tMesh), intent(inout) :: mesh

  integer(ip) :: dim
  integer(ip) :: bcMark, bcMarkAddr, bcCell, bcCellFace, numExtraBcCells
  integer(ip) :: cell, cellCell

  numExtraBcCells = 0

  ! ----------------------
  ! Generate extra connectivity.
  ! ----------------------
  do cell = 1, mesh%NumCells
    if (mesh%NumDims == 2) then
      
      ! ----------------------
      ! D2Q9 mesh:
      ! ----------------------
      !  8  3  5
      !   \ | /
      ! 2 --9-- 1
      !   / | \
      !  6  4  7
      ! ----------------------

      ! ----------------------
      ! Find cell "5":
      ! either 1 -> 3 or 3 -> 1.
      ! ----------------------
      associate(cellCellCell => mesh%CellToCell(5,cell))
        cellCell = mesh%CellToCell(1,cell)
        cellCellCell = mesh%CellToCell(3,cellCell)
        if (cellCellCell == 0) then
          cellCell = mesh%CellToCell(3,cell)
          cellCellCell = mesh%CellToCell(1,cellCell)
          if (cellCellCell == 0) then
            cellCellCell = mesh%NumAllCells + 1
            mesh%NumAllCells = cellCellCell
            mesh%BCMToCell = [mesh%BCMToCell, cellCellCell]
            mesh%BCMToCellFace = [mesh%BCMToCellFace, 6]
            numExtraBcCells = numExtraBcCells + 1
          end if
        end if
      end associate

      ! ----------------------
      ! Find cell "6":
      ! either 2 -> 4 or 4 -> 2.
      ! ----------------------
      associate(cellCellCell => mesh%CellToCell(6,cell))
        cellCell = mesh%CellToCell(2,cell)
        cellCellCell = mesh%CellToCell(4,cellCell)
        if (cellCellCell == 0) then
          cellCell = mesh%CellToCell(4,cell)
          cellCellCell = mesh%CellToCell(2,cellCell)
          if (cellCellCell == 0) then
            cellCellCell = mesh%NumAllCells + 1
            mesh%NumAllCells = cellCellCell
            mesh%BCMToCell = [mesh%BCMToCell, cellCellCell]
            mesh%BCMToCellFace = [mesh%BCMToCellFace, 5]
            numExtraBcCells = numExtraBcCells + 1
          end if
        end if
      end associate

      ! ----------------------
      ! Find cell "7":
      ! either 1 -> 4 or 4 -> 1.
      ! ----------------------
      associate(cellCellCell => mesh%CellToCell(7,cell))
        cellCell = mesh%CellToCell(1,cell)
        cellCellCell = mesh%CellToCell(4,cellCell)
        if (cellCellCell == 0) then
          cellCell = mesh%CellToCell(4,cell)
          cellCellCell = mesh%CellToCell(1,cellCell)
          if (cellCellCell == 0) then
            cellCellCell = mesh%NumAllCells + 1
            mesh%NumAllCells = cellCellCell
            mesh%BCMToCell = [mesh%BCMToCell, cellCellCell]
            mesh%BCMToCellFace = [mesh%BCMToCellFace, 8]
            numExtraBcCells = numExtraBcCells + 1
          end if
        end if
      end associate

      ! ----------------------
      ! Find cell "8":
      ! either 2 -> 3 or 3 -> 2.
      ! ----------------------
      associate(cellCellCell => mesh%CellToCell(8,cell))
        cellCell = mesh%CellToCell(2,cell)
        cellCellCell = mesh%CellToCell(3,cellCell)
        if (cellCellCell == 0) then
          cellCell = mesh%CellToCell(3,cell)
          cellCellCell = mesh%CellToCell(2,cellCell)
          if (cellCellCell == 0) then
            cellCellCell = mesh%NumAllCells + 1
            mesh%NumAllCells = cellCellCell
            mesh%BCMToCell = [mesh%BCMToCell, cellCellCell]
            mesh%BCMToCellFace = [mesh%BCMToCellFace, 7]
            numExtraBcCells = numExtraBcCells + 1
          end if
        end if
      end associate

      ! ----------------------
      ! "9" is the easy one.
      ! ----------------------
      mesh%CellToCell(9,cell) = cell

    else if (mesh%NumDims == 3) then

      ! ----------------------
      ! D3Q19 mesh.
      ! ----------------------
      error stop 'not implemented'

    end if
  end do

  if (numExtraBcCells /= 0) then

    !mesh%NumBCMs = mesh%NumBCMs + 1
    mesh%BCMs = [mesh%BCMs, size(mesh%BCMToCell) + 1]

  end if

end subroutine GenerateMeshExtraConnectivity 

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh connectivity in '.dot' format.
!! Output may be visualized with: 'neato -n -Tpng c2c.dot > c2c.png' 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine PrintMeshToNeato(mesh, file)
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: file

  integer(ip) :: unit
  integer(ip) :: cell, cellCell, cellFace
  integer(ip) :: bcMark, bcMarkAddr, bcCell, bcCellCell, bcCellFace

  integer(ip), parameter :: DPI = 72
  character(len=10), parameter :: SHAPES(*) = &
    & [character(len=10) :: 'box', 'diamond']
  character(len=10), parameter :: PALETTE(*) = &
    & [character(len=10) :: 'red', 'green', 'blue', 'fuchsia', 'yellow', 'wheat']
  
  if (mesh%NumDims /= 2) then
    error stop 'Only 2D meshes can be printed to Neato'
  end if

  print *, ''
  print *, '-=-=-=-=-=-=-=-'
  print *, 'Print to Neato: ', file
  print *, '-=-=-=-=-=-=-=-'
  print *, ''

  open(newunit=unit, file=file, status='replace')
  
  write(unit, "(A)") '# Generated by StormRuler/Mesh2Neato'
  write(unit, "(A)") 'digraph Mesh {'
  
  ! ----------------------
  ! Write cell MD indices of the cells as graph nodes.
  ! ----------------------
  ! Difference node shapes are use for different face
  ! directions of the BC/ghost cells faces.
  do cell = 1, mesh%NumCells
    associate(cellPosition => DPI*mesh%CellMDIndex(:,cell))

      write(unit, "('  ', A, '[label=', A, ' ', 'pos=', A, ']')") &
        & 'C'//I2S(cell), '"C"', & ! ID and label
        & '"'//I2S(cellPosition(1))//','//I2S(cellPosition(2))//'!"' ! pos

    end associate
  end do

  do bcMark = 1, mesh%NumBCMs
    do bcMarkAddr = mesh%BCMs(bcMark), mesh%BCMs(bcMark+1)-1
      bcCell = mesh%BCMToCell(bcMarkAddr)
      bcCellFace = mesh%BCMToCellFace(bcMarkAddr)
  
      do while(bcCell /= 0)
        associate(iBCCellPos => DPI*mesh%CellMDIndex(:,bcCell))

          write(unit, "('  ', A, '[label=', A, ' pos=', A, ' shape=', A, ']')") &
            & 'C'//I2S(bcCell), & ! ID
            & '"B'//I2S(bcMark)//'"', & ! label
            & '"'//I2S(iBCCellPos(1))//','//I2S(iBCCellPos(2))//'!"', & ! pos
            & trim(SHAPES((bcCellFace+1)/2)) ! shape
        
        end associate
        bcCell = mesh%CellToCell(bcCellFace, bcCell)
      end do
  
    end do
  end do
  
  ! ----------------------
  ! Write cell to cell connectivity as graph edges.
  ! ----------------------
  ! BC/ghost cells are colored by the BC mark.
  do cell = 1, mesh%NumCells
    do cellFace = 1, mesh%NumCellFaces
      ! Do not dump the periodic faces.
      if (mesh%IsCellFacePeriodic(cellFace, cell)) cycle
      
      cellCell = mesh%CellToCell(cellFace, cell)
      write(unit, "('  ', A, '->', A)") &
        & 'C'//I2S(cell), 'C'//I2S(cellCell) ! source and dest ID

    end do
  end do

  do bcMark = 1, mesh%NumBCMs
    do bcMarkAddr = mesh%BCMs(bcMark), mesh%BCMs(bcMark+1)-1
      bcCell = mesh%BCMToCell(bcMarkAddr)
      bcCellFace = mesh%BCMToCellFace(bcMarkAddr)

      do
        associate(color => PALETTE(bcMark))

          bcCellCell = mesh%CellToCell(Flip(bcCellFace), bcCell) 
          
          write(unit, "('  ', A, '->', A, ' [color=', A, ']')") &
            & 'C'//I2S(bcCell), & ! source ID
            & 'C'//I2S(bcCellCell), color ! dest ID and color

          bcCellCell = mesh%CellToCell(bcCellFace, bcCell) 
          if (bcCellCell == 0) exit

          write(unit, "('  ', A, '->', A, ' [color=', A, ']')") &
            & 'C'//I2S(bcCell), & ! source ID
            & 'C'//I2S(bcCellCell), color ! dest ID and color
            
        end associate

        bcCell = bcCellCell
      end do

    end do
  end do

  write(unit, "(A)") '}'
  
  ! ----------------------
  ! Close file and exit.
  ! ----------------------
  close(unit)

end subroutine PrintMeshToNeato

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh in Legacy VTK '.vtk' format.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine PrintMeshToLegacyVTK(mesh, file, fields)
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: file
  type(IOList), intent(in), optional :: fields
  
  integer(ip) :: unit
  integer(ip) :: numVtkCells, vtkCell
  integer(ip) :: cell, cellCell, cellFace, cellCellCell, cellCellFace
  logical, allocatable :: cellIsFullInternal(:)
  class(IOListItem), pointer :: item

  !! TODO: it feels like this thing can be seriously optimized..

  if ((mesh%NumDims /= 2).and.(mesh%NumDims /= 3)) then
    error stop 'Only 2D/3D meshes can be printed to Legacy VTK'
  end if

  print *, 'Print to Legacy VTK: ', file
  
  open(newunit=unit, file=file, status='replace')

  write(unit, "(A)") '# vtk DataFile Version 2.0'
  write(unit, "(A)") '# Generated by StormRuler/Mesh2VTK'
  write(unit, "(A)") 'ASCII'
  write(unit, "(A)") 'DATASET UNSTRUCTURED_GRID'
  write(unit, "(A)") ''
  
  ! ----------------------
  ! Write cell centers as VTK nodes.
  ! ----------------------
  write(unit, "('POINTS ', A, ' double')") I2S(mesh%NumCells)
  do cell = 1, mesh%NumCells
    associate(r => mesh%CellCenter(cell))
      if (mesh%NumDims == 2) then
        write(unit, "(A, ' ', A, ' ', A)") R2S(r(1)), R2S(r(2)), '0.0'
      else
        write(unit, "(A, ' ', A, ' ', A)") R2S(r(1)), R2S(r(2)), R2S(r(3))
      end if
    end associate
  end do
  write(unit, "(A)") ''
  
  ! ----------------------
  ! Precompute if cell is far away from boundaries.
  ! ----------------------
  allocate(cellIsFullInternal(mesh%NumCells))
  !$omp parallel do
  do cell = 1, mesh%NumCells
    cellIsFullInternal(cell) = &
      & all(mesh%CellToCell(:,cell) <= mesh%NumCells)
  end do
  !$omp end parallel do
  
  ! ----------------------
  ! Write VTK cells.
  ! VTK cell connectivity is generated on the fly via parity algorithm.
  ! ----------------------
#$do writePass = 1, 2

  numVtkCells = 0
  !$omp parallel do reduction(+:numVtkCells) &
  !$omp private(cellFace, cellCell, cellCellFace, cellCellCell) if($writePass == 1)
  do cell = 1, mesh%NumCells

    ! ----------------------
    ! Locate the second cell (VTK cell node).
    ! ----------------------
    do cellFace = 1, mesh%NumCellFaces
      if (mesh%IsCellFacePeriodic(cellFace, cell)) cycle
      cellCell = mesh%CellToCell(cellFace, cell)
      if (cellCell > mesh%NumCells) cycle

      ! ----------------------
      ! Locate the third cell (VTK cell node).
      ! ----------------------
      do cellCellFace = 1, mesh%NumCellFaces
        if (cellCellFace == cellFace) cycle
        if (mesh%IsCellFacePeriodic(cellCellFace, cell)) cycle
        cellCellCell = mesh%CellToCell(cellCellFace, cell)
        if (cellCellCell > mesh%NumCells) cycle

        ! ----------------------
        ! Generate the VTK cell (triangle).
        ! ----------------------
        ! Denote cells with same index components parity as primary.
        ! Primary cells are used to generate VTK cells around them.
        ! Secondary cells are used to generate VTK cells near boundaries.
        associate(cellIndex => mesh%CellMDIndex(:,cell))
          if (mod(cellIndex(1) - cellIndex(2), 2) /= 0) then
            ! Skip secondary cells far away from boundaries.
            if ( cellIsFullInternal(cellCell).or. &
              &  cellIsFullInternal(cellCellCell) ) cycle
          end if
        end associate

        numVtkCells = numVtkCells + 1
#$if writePass == 2
        write(unit, "('3 ', A, ' ', A, ' ', A)") &
          & I2S(cell-1), I2S(cellCell-1), I2S(cellCellCell-1)
#$end if

      end do
    end do
  end do
  !$omp end parallel do

  ! ----------------------
  ! Write a number of the VTK cells as the result of the first pass.
  ! ----------------------
#$if writePass == 1
  associate(numVTKCellNodes => numVtkCells*(mesh%NumDims + 2))
    write(unit, "('CELLS ', A, ' ', A)") I2S(numVtkCells), I2S(numVTKCellNodes)
  end associate
#$end if

#$end do

  ! ----------------------
  ! Write VTK cell types.
  ! ----------------------
  write(unit, "('CELL_TYPES ', A)") I2S(numVtkCells)
  do vtkCell = 1, numVtkCells
    write(unit, "(A)") '5' ! 5 is triangle.
  end do
  write(unit, "(A)") ''

  ! ----------------------
  ! Write fields.
  ! ----------------------
  write(unit, "('POINT_DATA ', A)") I2S(mesh%NumCells)
  if (present(fields)) then
    item => fields%first
    do while(associated(item))
      select type(item)
        ! ----------------------
        ! Scalar field.
        ! ----------------------
        class is(IOListItem$0)
          write(unit, "('SCALARS ', A, ' double 1')") item%name 
          write(unit, "(A)") 'LOOKUP_TABLE default'
          do cell = 1, mesh%NumCells
            write(unit, "(A)") R2S(item%values(cell))
          end do

        ! ----------------------
        ! Vector field.
        ! ----------------------
        class is(IOListItem$1)
          write(unit, "('VECTORS ', A, ' double')") item%name 
          do cell = 1, mesh%NumCells
            associate(v => item%values(:,cell))
              write(unit, '(A, " ", A, " ", A)') R2S(v(1)), R2S(v(2)), '0.0'
            end associate
          end do

        end select
      write(unit, "(A)") ''
      item => item%next
    end do
  end if

  ! ----------------------
  ! Close file and exit.
  ! ----------------------
  close(unit)

end subroutine PrintMeshToLegacyVTK

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -----------------------------------------------------------------  
!! Initialize a 2D rectangular mesh.
subroutine tMesh_InitRect(mesh, xDelta, xNumCells, xPeriodic &
                               ,yDelta, yNumCells, yPeriodic, numLayers)
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: xDelta, yDelta
  integer(ip), intent(in) :: xNumCells, yNumCells
  logical, intent(in) :: xPeriodic, yPeriodic
  integer(ip), intent(in), optional :: numLayers
  
  integer(ip) :: xNumLayers, yNumLayers
  integer(ip), allocatable :: cellToIndex(:,:)
  ! ----------------------
  ! Set up amount of layers.
  mesh%NumDims = 2
  block
    if (present(numLayers)) then
      xNumLayers = numLayers
      yNumLayers = numLayers
    end if
    if (.not.present(numLayers).or.xPeriodic) xNumLayers = 1
    if (.not.present(numLayers).or.yPeriodic) yNumLayers = 1
  end block
  ! ----------------------
  ! Build a 2D to 1D index mapping.
  block
    integer(ip) :: xCell, yCell
    allocate(cellToIndex(1-xNumLayers:xNumCells+xNumLayers &
                       , 1-yNumLayers:yNumCells+yNumLayers))
    ! ----------------------
    ! Generate indices for the interior cells.
    associate(iCell => mesh%NumCells)
      iCell = 0
      do yCell = 1, yNumCells
        do xCell = 1, xNumCells
          iCell = iCell + 1
          cellToIndex(xCell, yCell) = iCell
        end do
      end do
    end associate
    ! ----------------------
    ! Generate indices for periodicicty or ghost cells.
    associate(iCell => mesh%NumAllCells)
      iCell = mesh%NumCells
      if (xPeriodic) then
        xNumLayers = 0
        !$omp parallel do private(yCell)
        do yCell = 1, yNumCells
          cellToIndex(0, yCell) = cellToIndex(xNumCells, yCell)
          cellToIndex(xNumCells+1, yCell) = cellToIndex(1, yCell)
        end do
        !$omp end parallel do
      else
        do yCell = 1, yNumCells
          do xCell = 1-xNumLayers, 0
            iCell = iCell + 1
            cellToIndex(xCell, yCell) = iCell
          end do
          do xCell = xNumCells+1, xNumCells+xNumLayers
            iCell = iCell + 1
            cellToIndex(xCell, yCell) = iCell
          end do
        end do
      end if
      if (yPeriodic) then
        yNumLayers = 0
        !$omp parallel do private(xCell)
        do xCell = 1, xNumCells
          cellToIndex(xCell, 0) = cellToIndex(xCell, yNumCells)
          cellToIndex(xCell, yNumCells+1) = cellToIndex(xCell, 1)
        end do
        !$omp end parallel do
      else
        do xCell = 1, xNumCells
          do yCell = 1-yNumLayers, 0
            iCell = iCell + 1
            cellToIndex(xCell, yCell) = iCell
          end do
          do yCell = yNumCells+1, yNumCells+yNumLayers
            iCell = iCell + 1
            cellToIndex(xCell, yCell) = iCell
          end do
        end do
      end if
    end associate
  end block
  ! ----------------------
  ! Fill the cell geometry and connectivity.
  block
    integer(ip) :: iCell, xCell, yCell
    allocate(mesh%dl(1:4))
    mesh%dl(:) = [xDelta, xDelta, yDelta, yDelta]
    mesh%MDIndexBounds = [xNumCells, yNumCells]
    allocate(mesh%dr(1:2, 1:9))
    mesh%dr(:,1) = [+1.0_dp, 0.0_dp]
    mesh%dr(:,2) = [-1.0_dp, 0.0_dp]
    mesh%dr(:,3) = [0.0_dp, +1.0_dp]
    mesh%dr(:,4) = [0.0_dp, -1.0_dp]
    mesh%dr(:,5) = [+1.0_dp, +1.0_dp]
    mesh%dr(:,6) = [-1.0_dp, -1.0_dp]
    mesh%dr(:,7) = [+1.0_dp, -1.0_dp]
    mesh%dr(:,8) = [-1.0_dp, +1.0_dp]
    mesh%dr(:,9) = [0.0_dp, 0.0_dp]
    mesh%NumConns = 9
    mesh%NumCellFaces = 4
    allocate(mesh%CellToCell(9, mesh%NumAllCells))
    allocate(mesh%mIsCellFacePeriodic(4, mesh%NumAllCells))
    allocate(mesh%CellMDIndex(2, mesh%NumAllCells))
    ! ----------------------
    ! Fill the cell information for the interior cells.
    !$omp parallel do private(iCell, yCell, xCell) collapse(2)
    do yCell = 1, yNumCells
      do xCell = 1, xNumCells
        iCell = cellToIndex(xCell, yCell)
        mesh%CellMDIndex(:,iCell) = [xCell, yCell]
        !---
        mesh%CellToCell(:,iCell) = &
          & [cellToIndex(xCell+1, yCell), &
          &  cellToIndex(xCell-1, yCell), &
          &  cellToIndex(xCell, yCell+1), &
          &  cellToIndex(xCell, yCell-1)]
        !---
        mesh%mIsCellFacePeriodic(:,iCell) = .false.
        if (xCell==1) mesh%mIsCellFacePeriodic(2, iCell) = .true.
        if (yCell==1) mesh%mIsCellFacePeriodic(4, iCell) = .true.
        if (xCell==xNumCells) mesh%mIsCellFacePeriodic(1, iCell) = .true.
        if (yCell==yNumCells) mesh%mIsCellFacePeriodic(3, iCell) = .true.
      end do
    end do
    !$omp end parallel do
    ! ----------------------
    ! Fill the cell information for the ghost cells.
    ! (Remember, the the ghost cells have a single connection to the interior cell.)
    !$omp parallel do private(iCell, yCell, xCell)
    do yCell = 1, yNumCells
      do xCell = 1-xNumLayers, 0
        iCell = cellToIndex(xCell, yCell)
        mesh%CellToCell(iCell, :) &
          = [ cellToIndex(1, yCell), 0, 0, 0 ]
      end do
      do xCell = xNumCells+1, xNumCells+xNumLayers 
        iCell = cellToIndex(xCell, yCell)
        mesh%CellToCell(iCell, :) &
          = [ 0, cellToIndex(xNumCells, yCell), 0, 0 ]
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(iCell, yCell, xCell)
    do xCell = 1, xNumCells
      do yCell = 1-yNumLayers, 0
        iCell = cellToIndex(xCell, yCell)
        mesh%CellToCell(iCell, :) &
          = [ 0, 0, cellToIndex(xCell, 1), 0 ]
      end do
      do yCell = yNumCells+1, yNumCells+yNumLayers
        iCell = cellToIndex(xCell, yCell)
        mesh%CellToCell(iCell, :) &
          = [ 0, 0, 0, cellToIndex(xCell, yNumCells) ]
      end do
    end do
    !$omp end parallel do
  end block

  block
    integer :: iCell, xCell
    do iCell = 1, mesh%NumCells

      mesh%CellToCell(9,iCell) = iCell

      xCell = mesh%CellToCell(1,iCell)
      mesh%CellToCell(5,iCell) = mesh%CellToCell(3,xCell)
      mesh%CellToCell(7,iCell) = mesh%CellToCell(4,xCell)

      xCell = mesh%CellToCell(2,iCell)
      mesh%CellToCell(6,iCell) = mesh%CellToCell(4,xCell)
      mesh%CellToCell(8,iCell) = mesh%CellToCell(3,xCell)

    end do
  end block

end subroutine tMesh_InitRect

end module StormRuler_Mesh
