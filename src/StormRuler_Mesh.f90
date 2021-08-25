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
use StormRuler_IO, only: IOList, IOListItem, @{IOListItem$$@|@0, NUM_RANKS}@

use, intrinsic :: iso_fortran_env, only: error_unit
#$if HAS_OpenMP
use :: omp_lib
#$endif

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

abstract interface
  subroutine tKernelFunc(iCell)
    import ip
    integer(ip), intent(in) :: iCell
  end subroutine tKernelFunc
  subroutine tBlockKernelFunc(firstCell, lastCell)
    import ip
    integer(ip), intent(in) :: firstCell, lastCell
  end subroutine tBlockKernelFunc
  function tReduceKernelFunc(iCell) result(r)
    import ip, dp
    integer(ip), intent(in) :: iCell
    real(dp) :: r
  end function tReduceKernelFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Semi-structured multidimensional mesh.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
type :: tMesh
  ! ----------------------
  ! Mesh dimension.
  ! ----------------------
  integer(ip) :: Dim

  ! ----------------------
  ! Number of faces (edges in 2D) per each cell.
  ! ----------------------
  integer(ip) :: NumCellFaces

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
  ! Shape is [1, NumFaces]×[1, NumAllCells].
  ! ----------------------
  integer(ip), allocatable :: CellToCell(:,:)
  ! ----------------------
  ! Logical table for the periodic face connections.
  ! Shape is [1, NumFaces]×[1, NumAllCells].
  ! ----------------------
  logical, allocatable, private :: mCellFacePeriodic(:,:)

  ! ----------------------
  ! Multidimensional index bounds table.
  ! Shape is [1, Dim].
  ! ----------------------
  integer(ip), allocatable :: MDIndexBounds(:)
  ! ----------------------
  ! Cell multidimensional index table.
  ! Shape is [1, Dim]×[1, NumCells].
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
  ! Temporal step value.
  ! ----------------------
  real(dp) :: dt
  ! ----------------------
  ! Distance between centers of the adjacent cells per face.
  ! Shape is [1, NumFaces].
  ! ----------------------
  real(dp), allocatable :: dl(:)
  ! ----------------------
  ! Difference between centers of the adjacent cells per face.
  ! Shape is [1, Dim]×[1, NumFaces].
  ! ----------------------
  real(dp), allocatable :: dr(:,:)
  ! ----------------------
  ! Cell center coordinates.
  ! Shape is [1, Dim]×[1, NumCells].
  ! ----------------------
  real(dp), allocatable :: mCellCenter(:,:)

  ! ----------------------
  ! Cell range indicators.
  ! ----------------------
  logical, private :: mParallel = .true.
#$if HAS_OpenMP
  integer(ip), allocatable, private :: mCellRange(:,:)
#$else
  integer(ip), private :: mCellRange(2)
#$end if

contains

  ! ----------------------
  ! Mesh walkthough subroutines.
  ! ----------------------
  procedure :: SetRange => tMesh_SetRange
  procedure :: FirstCell => tMesh_FirstCell
  procedure :: LastCell => tMesh_LastCell
  procedure :: RunCellKernel => tMesh_RunCellKernel
  procedure :: RunCellKernel_Sum => tMesh_RunCellKernel_Sum
  procedure :: RunCellKernel_Min => tMesh_RunCellKernel_Min
  procedure :: RunCellKernel_Max => tMesh_RunCellKernel_Max
  procedure :: RunCellKernel_Block => tMesh_RunCellKernel_Block
  procedure :: RunCellKernel_Forward => tMesh_RunCellKernel_Forward
  procedure :: RunCellKernel_Backward => tMesh_RunCellKernel_Backward

  ! ----------------------
  ! Field wrappers.
  ! ----------------------
  procedure :: CellFacePeriodic => tMesh_CellFacePeriodic
  procedure :: CellCenter => tMesh_CellCenter

  ! ----------------------
  ! Initializers.
  ! ----------------------
  procedure :: InitRect => tMesh_InitRect
  procedure :: InitFromImage2D => tMesh_InitFromImage2D
  procedure :: InitFromImage3D => tMesh_InitFromImage3D
  procedure, private :: GenerateBCCells => tMesh_GenerateBCCells

  ! ----------------------
  ! Printers.
  ! ----------------------
  procedure :: PrintTo_Neato => tMesh_PrintTo_Neato
  procedure :: PrintTo_LegacyVTK => tMesh_PrintTo_LegacyVTK
end type tMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! ----------------------------------------------------------------- !!
subroutine tMesh_SetRange(mesh, firstCell, lastCell, parallel)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in), optional :: firstCell, lastCell
  logical, optional :: parallel
  ! >>>>>>>>>>>>>>>>>>>>>>

#$if HAS_OpenMP
  ! Preallocate and reset the ranges.
  if (.not.allocated(mesh%mCellRange)) then
    !$omp single
    allocate(mesh%mCellRange(2, omp_get_max_threads()))
    mesh%mCellRange(1,:) = 1
    mesh%mCellRange(2,:) = mesh%NumCells
    !$omp end single
  end if
  if (present(firstCell)) then
    mesh%mParallel = .false.
    mesh%mCellRange(:,omp_get_thread_num()+1) = firstCell
    if (present(lastCell)) then
      mesh%mCellRange(2,omp_get_thread_num()+1) = lastCell
    end if
  else
    if (omp_in_parallel()) then
      write(error_unit, *) &
        & 'SetRange cannot be called without ', &
        & 'arguments inside of the parallel section.'
      error stop 2
    end if
    mesh%mParallel = .true.
    if (present(parallel)) then
      mesh%mParallel = parallel
    end if
    mesh%mCellRange(1,:) = 1
    mesh%mCellRange(2,:) = mesh%NumCells
  end if
#$else
  if (present(firstCell)) then
    mesh%mCellRange(:) = firstCell
    if (present(lastCell)) then
      mesh%mCellRange(2) = lastCell
    end if
  else
    mesh%mCellRange(1) = 1
    mesh%mCellRange(2) = mesh%NumCells
  end if
#$end if
end subroutine tMesh_SetRange

!! ----------------------------------------------------------------- !!
!! Indicate if the computation should be performed in parallel.
!! ----------------------------------------------------------------- !!
pure logical function tMesh_Parallel(mesh)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  ! >>>>>>>>>>>>>>>>>>>>>>

  tMesh_Parallel = mesh%mParallel
end function tMesh_Parallel

!! ----------------------------------------------------------------- !!
!! Get the lower cell bound for the current computation.
!! ----------------------------------------------------------------- !!
integer(ip) function tMesh_FirstCell(mesh)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  ! >>>>>>>>>>>>>>>>>>>>>>

#$if HAS_OpenMP
  tMesh_FirstCell = mesh%mCellRange(1, omp_get_thread_num()+1)
#$else
  tMesh_FirstCell = mesh%mCellRange(1)
#$end if
end function tMesh_FirstCell

!! ----------------------------------------------------------------- !!
!! Get the upper cell bound for the current computation.
!! ----------------------------------------------------------------- !!
integer(ip) function tMesh_LastCell(mesh)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  ! >>>>>>>>>>>>>>>>>>>>>>

#$if HAS_OpenMP
  tMesh_LastCell = mesh%mCellRange(2, omp_get_thread_num()+1)
#$else
  tMesh_LastCell = mesh%mCellRange(2)
#$end if
end function tMesh_LastCell

!! ----------------------------------------------------------------- !!
!! Launch a cell kernel.
!! ----------------------------------------------------------------- !!
subroutine tMesh_RunCellKernel(mesh, Kernel)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  procedure(tKernelFunc) :: Kernel
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer :: iCell

  if (mesh%mParallel) then
    !$omp parallel do schedule(static)
    do iCell = mesh%FirstCell(), mesh%LastCell()
      call Kernel(iCell)
    end do
    !$omp end parallel do
  else
    do iCell = mesh%FirstCell(), mesh%LastCell()
      call Kernel(iCell)
    end do
  end if
end subroutine tMesh_RunCellKernel

!! ----------------------------------------------------------------- !!
!! Launch a SUM-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function tMesh_RunCellKernel_Sum(mesh, Kernel) result(sum)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: sum
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer :: iCell

  sum = 0.0_dp
  if (mesh%mParallel) then
    !$omp parallel do schedule(static) reduction(+:sum)
    do iCell = mesh%FirstCell(), mesh%LastCell()
      sum = sum + Kernel(iCell)
    end do
    !$omp end parallel do
  else
    do iCell = mesh%FirstCell(), mesh%LastCell()
      sum = sum + Kernel(iCell)
    end do
  end if
end function tMesh_RunCellKernel_Sum

!! ----------------------------------------------------------------- !!
!! Launch a MIN-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function tMesh_RunCellKernel_Min(mesh, Kernel) result(mMin)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: mMin
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer :: iCell

  mMin = +huge(mMin)
  if (mesh%mParallel) then
    !$omp parallel do schedule(static) reduction(min:mMin)
    do iCell = mesh%FirstCell(), mesh%LastCell()
      mMin = min(mMin, Kernel(iCell))
    end do
    !$omp end parallel do
  else
    do iCell = mesh%FirstCell(), mesh%LastCell()
      mMin = min(mMin, Kernel(iCell))
    end do
  end if
end function tMesh_RunCellKernel_Min

!! ----------------------------------------------------------------- !!
!! Launch a MAX-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function tMesh_RunCellKernel_Max(mesh, Kernel) result(mMax)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: mMax
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer :: iCell

  mMax = -huge(mMax)
  if (mesh%mParallel) then
    !$omp parallel do schedule(static) reduction(max:mMax)
    do iCell = mesh%FirstCell(), mesh%LastCell()
      mMax = max(mMax, Kernel(iCell))
    end do
    !$omp end parallel do
  else
    do iCell = mesh%FirstCell(), mesh%LastCell()
      mMax = max(mMax, Kernel(iCell))
    end do
  end if
end function tMesh_RunCellKernel_Max

!! ----------------------------------------------------------------- !!
!! Launch a block-kernel.
!! ----------------------------------------------------------------- !!
subroutine tMesh_RunCellKernel_Block(mesh, BlockKernel)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  procedure(tBlockKernelFunc) :: BlockKernel
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: i
  integer(ip), allocatable :: ranges(:)

#$if HAS_OpenMP
  ! ----------------------
  ! Launch the whole range as a block if sequential mode is requested.
  ! ----------------------
  if (.not.mesh%mParallel) then
    call BlockKernel(mesh%FirstCell(), mesh%LastCell())
    return
  end if

  ! ----------------------
  ! Compute block range for each thread.
  ! ----------------------
  allocate(ranges(omp_get_max_threads()))
  associate(size => ( mesh%LastCell() - mesh%FirstCell() + 1 ) )
    ranges(:) = size/omp_get_max_threads()
    associate(remainder => ranges(:mod(size, omp_get_max_threads())))
      remainder(:) = remainder(:) + 1
    end associate
  end associate
  ranges = [mesh%FirstCell(), ranges]
  do i = 1, omp_get_max_threads()
    ranges(i + 1) = ranges(i) + ranges(i + 1) 
  end do

  ! ----------------------
  ! Lauch threads.
  ! ----------------------
  !$omp parallel default(none) shared(ranges)
  associate(iThread => omp_get_thread_num() + 1)
    call BlockKernel(ranges(iThread), ranges(iThread + 1) - 1)
  end associate
  !$omp end parallel
#$else
  call BlockKernel(mesh%FirstCell(), mesh%LastCell())
#$endif
end subroutine tMesh_RunCellKernel_Block

subroutine tMesh_RunCellKernel_Forward(mesh, Kernel)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  procedure(tKernelFunc) :: Kernel
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer :: iCell

  do iCell = mesh%FirstCell(), mesh%LastCell()
    call Kernel(iCell)
  end do
end subroutine tMesh_RunCellKernel_Forward

subroutine tMesh_RunCellKernel_Backward(mesh, Kernel)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  procedure(tKernelFunc) :: Kernel
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer :: iCell

  do iCell = mesh%LastCell(), mesh%FirstCell(), -1
    call Kernel(iCell)
  end do
end subroutine tMesh_RunCellKernel_Backward

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! ----------------------------------------------------------------- !!
logical pure function tMesh_CellFacePeriodic(mesh, iCellFace, iCell)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: iCell, iCellFace
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  tMesh_CellFacePeriodic = &
    & allocated(mesh%mCellFacePeriodic).and. &
    & mesh%mCellFacePeriodic(iCellFace, iCell)
end function tMesh_CellFacePeriodic

!! ----------------------------------------------------------------- !!
!! Get cell center.
!! ----------------------------------------------------------------- !!
pure function tMesh_CellCenter(mesh, iCell)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: iCell
  real(dp) :: tMesh_CellCenter(mesh%Dim)
  ! >>>>>>>>>>>>>>>>>>>>>>

  tMesh_CellCenter = &
    & mesh%dl(::2)*mesh%CellMDIndex(:,iCell)
end function tMesh_CellCenter

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

  !! TODO:
  allocate(mesh%dl(1:4))
  mesh%dl(:) = [0.01_dp, 0.01_dp, 0.01_dp, 0.01_dp]

  ! ----------------------
  ! Process image data in order to:
  ! 1. Generate initial "Cell MD Index" ⟺ "Cell (Plain Index)" table 
  !    for interior cells and compute number of interior and boundary cells.
  ! 2. Generate and pre-fill "BCM" ⟺ "Num BC Cells per BCM" CSR-style table.
  ! ----------------------
  mesh%Dim = $dim
  mesh%NumCellFaces = ${2*dim}$
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
          write(error_unit, *) 'INVALID IMAGE COLOR', IntToPixel(image($iiCellMD))
          error stop
        end if
        mesh%BCMs(iBCM+1) = mesh%BCMs(iBCM+1) + 1
      end if
#$end for
    end if
#$do rank = dim, 1, -1
  end do
#$end do
  
  ! ----------------------
  ! 3. TODO: apply plain indices sorting.
  ! ----------------------
  mesh%NumCells = iCell
  mesh%NumBCCells = sum(mesh%BCMs)
  mesh%NumAllCells = mesh%NumCells+mesh%NumBCCells*numBCLayers
  print *, 'NC:', mesh%NumCells, 'NBC:', mesh%NumBCCells, 'NAC:', mesh%NumAllCells

  ! ----------------------
  ! 4. Allocate the connectivity tables and fill them with data,
  !    inverting the "Cell MD Index" ⟺ "Cell (Plain Index)" table.
  ! ----------------------
  associate(nac => mesh%NumAllCells, ncf => mesh%NumCellFaces)
    allocate(mesh%CellMDIndex(mesh%Dim, nac))
    allocate(mesh%CellToCell(ncf, nac)); mesh%CellToCell(:,:) = 0
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
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in) :: numBCLayers
  ! >>>>>>>>>>>>>>>>>>>>>>
  
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
end subroutine tMesh_GenerateBCCells

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh connectivity in '.dot' format.
!! Output may be visualized with: 'neato -n -Tpng c2c.dot > c2c.png' 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine tMesh_PrintTo_Neato(mesh, file)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: file
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: unit
  integer(ip) :: iCell, iiCell, iCellFace
  integer(ip) :: iBCM, iBCMPtr, iBCCell, iiBCCell, iBCCellFace
  ! ----------------------
  integer(ip), parameter :: DPI = 72
  character(len=10), parameter :: SHAPES(2) = &
    & [character(len=10) :: 'box', 'diamond']
  character(len=10), parameter :: PALETTE(3) = &
    & [character(len=10) :: 'red', 'green', 'blue']
  
  ! ----------------------
  if (mesh%Dim /= 2) then
    error stop 'Only 2D meshes can be printed to Neato'
  end if
  print *, 'Print to Neato: ', file
  
  ! ----------------------
  open(newunit=unit, file=file, status='replace')
  write(unit, "(A)") '# Generated by StormRuler/Mesh2Neato'
  write(unit, "(A)") 'digraph Mesh {'
  
  ! ----------------------
  ! Write cell MD indices of the cells as graph nodes.
  ! Difference node shapes are use for different face 
  ! directions of the BC/ghost cells faces.
  ! ----------------------
  do iCell = 1, mesh%NumCells
    associate(iCellPos => DPI*mesh%CellMDIndex(:,iCell))
      100 format('  ', A, '[label=', A, ' ', 'pos=', A, ']')
      write(unit, 100) 'C'//I2S(iCell), &
        & '"C"', '"'//I2S(iCellPos(1))//', '//I2S(iCellPos(2))//'!"'
    end associate
  end do
  do iBCM = 1, mesh%NumBCMs
    do iBCMPtr = mesh%BCMs(iBCM), mesh%BCMs(iBCM+1)-1
      iBCCell = mesh%BCMToCell(iBCMPtr)
      iBCCellFace = mesh%BCMToCellFace(iBCMPtr)
      do while(iBCCell /= 0)
        associate(iBCCellPos => DPI*mesh%CellMDIndex(:,iBCCell))
          200 format('  ', A, '[label=', A, ' pos=', A, ' shape=', A, ']')
          write(unit, 200) 'C'//I2S(iBCCell), '"B'//I2S(iBCM)//'"', &
            & '"'//I2S(iBCCellPos(1))//', '//I2S(iBCCellPos(2))//'!"', &
            & trim(SHAPES((iBCCellFace+1)/2))
        end associate
        iBCCell = mesh%CellToCell(iBCCellFace, iBCCell)
      end do
    end do
  end do
  
  ! ----------------------
  ! Write cell to cell connectivity as graph edges.
  ! BC/ghost cells are colored by the BC mark.
  ! ----------------------
  do iCell = 1, mesh%NumCells
    do iCellFace = 1, mesh%NumCellFaces
      if (mesh%CellFacePeriodic(iCellFace, iCell)) cycle
      iiCell = mesh%CellToCell(iCellFace, iCell)
      write(unit, "('  ', A, '->', A)") 'C'//I2S(iCell), 'C'//I2S(iiCell)
    end do
  end do
  do iBCM = 1, mesh%NumBCMs
    do iBCMPtr = mesh%BCMs(iBCM), mesh%BCMs(iBCM+1)-1
      iBCCell = mesh%BCMToCell(iBCMPtr)
      iBCCellFace = mesh%BCMToCellFace(iBCMPtr)
      do
        associate(color => PALETTE(iBCM))
          iiBCCell = mesh%CellToCell(Flip(iBCCellFace), iBCCell) 
          300 format('  ', A, '->', A, ' [color=', A, ']')
          write(unit, 300) 'C'//I2S(iBCCell), 'C'//I2S(iiBCCell), trim(color)
          iiBCCell = mesh%CellToCell(iBCCellFace, iBCCell); if (iiBCCell == 0) exit
          write(unit, 300) 'C'//I2S(iBCCell), 'C'//I2S(iiBCCell), trim(color)
        end associate
        iBCCell = iiBCCell
      end do
    end do
  end do
  
  ! ----------------------
  ! Close file and exit.
  ! ----------------------
  write(unit, "(A)") '}'
  close(unit)
end subroutine tMesh_PrintTo_Neato

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh in Legacy VTK '.vtk' format.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine tMesh_PrintTo_LegacyVTK(mesh, file, fields)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: file
  type(IOList), intent(in), optional :: fields
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: unit
  integer(ip) :: numVTKCells, iVTKCell
  integer(ip) :: iCell, iiCell, iCellFace, jiCell, jCellFace
  logical, allocatable :: cellToCellIsAllInternal(:)
  class(IOListItem), pointer :: item
  
  ! ----------------------
  if ((mesh%Dim /= 2).and.(mesh%Dim /= 3)) then
    error stop 'Only 2D/3D meshes can be printed to Legacy VTK'
  end if
  print *, 'Print to Legacy VTK: ', file
  
  ! ----------------------
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
  do iCell = 1, mesh%NumCells
    associate(r => mesh%CellCenter(iCell))
      if (mesh%Dim == 2) then
        write(unit, "(A, ' ', A, ' ', A)") &
          & R2S(r(1)), R2S(r(2)), '0.0'
      else
        write(unit, "(A, ' ', A, ' ', A)") &
          & R2S(r(1)), R2S(r(2)), R2S(r(3))
      end if
    end associate
  end do
  write(unit, "(A)") ''
  
  ! ----------------------
  ! Precompute if cell is far away from boundaries.
  ! ----------------------
  allocate(cellToCellIsAllInternal(mesh%NumCells))
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    cellToCellIsAllInternal(iCell) = &
      & all(mesh%CellToCell(:,iCell) <= mesh%NumCells)
  end do
  !$omp end parallel do
  
  ! ----------------------
  ! Write VTK cells.
  ! VTK cell connectivity is generated on the fly via parity algorithm.
  ! ----------------------
#$do writePass = 1, 2
#$if writePass == 1
  numVTKCells = 0
  !$omp parallel do reduction(+:numVTKCells) &
  !$omp private(iCellFace, iiCell, jCellFace, jiCell)
#$end if
  do iCell = 1, mesh%NumCells
    ! ----------------------
    ! Locate the first cell (VTK cell node).
    ! ----------------------
    do iCellFace = 1, mesh%NumCellFaces
      if (mesh%CellFacePeriodic(iCellFace, iCell)) cycle
      iiCell = mesh%CellToCell(iCellFace, iCell)
      if (iiCell > mesh%NumCells) cycle

      ! ----------------------
      ! Locate the second cell (VTK cell node).
      ! ----------------------
      do jCellFace = 1, mesh%NumCellFaces
        if (jCellFace == iCellFace) cycle
        if (mesh%CellFacePeriodic(jCellFace, iCell)) cycle
        jiCell = mesh%CellToCell(jCellFace, iCell)
        if (jiCell > mesh%NumCells) cycle
    
        ! ----------------------
        ! Generate the VTK cell.
        ! ----------------------
        ! Denote cells with same MD index components parity as primary.
        ! Primary cells are used to generate VTK cells around them.
        ! Secondary cells are used to generate VTK cells near boundaries.
        associate(iCellMD => mesh%CellMDIndex(:,iCell))
          if (mod(iCellMD(1)-iCellMD(2), 2) /= 0) then
            ! Skip secondary cells far away from boundaries.
            if ( cellToCellIsAllInternal(iiCell).or. &
              &  cellToCellIsAllInternal(jiCell) ) cycle
          end if
        end associate
#$if writePass == 1
        numVTKCells = numVTKCells + 1
#$else
        write(unit, "('3 ', A, ' ', A, ' ', A)") I2S(iCell-1), I2S(iiCell-1), I2S(jiCell-1)
#$end if
      end do
    end do
  end do
#$if writePass == 1
  !$omp end parallel do
  associate(numVTKCellNodes => numVTKCells*(mesh%Dim+2))
    write(unit, "('CELLS ', A, ' ', A)") I2S(numVTKCells), I2S(numVTKCellNodes)
  end associate
#$end if
#$end do
  write(unit, "('CELL_TYPES ', A)") I2S(numVTKCells)
  do iVTKCell = 1, numVTKCells
    write(unit, "(A)") '5'
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
        ! Scalar field.
        class is(IOListItem$0)
          write(unit, "('SCALARS ', A, ' double 1')") item%name 
          write(unit, "(A)") 'LOOKUP_TABLE default'
          do iCell = 1, mesh%NumCells
            write(unit, "(A)") R2S(item%values(iCell))
          end do
        ! ----------------------
        ! Vector field.
        class is(IOListItem$1)
          write(unit, "('VECTORS ', A, ' double')") item%name 
          do iCell = 1, mesh%NumCells
            associate(v => item%values(:,iCell))
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
end subroutine tMesh_PrintTo_LegacyVTK

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
  mesh%Dim = 2
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
    allocate(mesh%dr(1:2, 1:4))
    mesh%dr(:,1) = [mesh%dl(1), 0.0_dp]
    mesh%dr(:,2) = [mesh%dl(1), 0.0_dp]
    mesh%dr(:,3) = [0.0_dp, mesh%dl(2)]
    mesh%dr(:,4) = [0.0_dp, mesh%dl(2)]
    mesh%NumCellFaces = 4
    allocate(mesh%mCellCenter(1:mesh%NumAllCells, 1:3))
    allocate(mesh%CellToCell(4, mesh%NumAllCells))
    allocate(mesh%mCellFacePeriodic(4, mesh%NumAllCells))
    allocate(mesh%CellMDIndex(2, mesh%NumAllCells))
    ! ----------------------
    ! Fill the cell information for the interior cells.
    !$omp parallel do private(iCell, yCell, xCell) collapse(2)
    do yCell = 1, yNumCells
      do xCell = 1, xNumCells
        iCell = cellToIndex(xCell, yCell)
        mesh%CellMDIndex(:,iCell) = [xCell, yCell]
        mesh%mCellCenter(iCell, :) = [xDelta*(xCell+0.5_dp)&
                                    ,yDelta*(yCell+0.5_dp), 0.0_dp]
        !---
        mesh%CellToCell(:,iCell) = &
          & [cellToIndex(xCell+1, yCell), &
          &  cellToIndex(xCell-1, yCell), &
          &  cellToIndex(xCell, yCell+1), &
          &  cellToIndex(xCell, yCell-1)]
        !---
        mesh%mCellFacePeriodic(:,iCell) = .false.
        if (xCell==1) mesh%mCellFacePeriodic(2, iCell) = .true.
        if (yCell==1) mesh%mCellFacePeriodic(4, iCell) = .true.
        if (xCell==xNumCells) mesh%mCellFacePeriodic(1, iCell) = .true.
        if (yCell==yNumCells) mesh%mCellFacePeriodic(3, iCell) = .true.
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
        mesh%mCellCenter(iCell, :) = [xDelta*(xCell+0.5_dp)&
                              ,yDelta*(yCell+0.5_dp), 0.0_dp]
        mesh%CellToCell(iCell, :) &
          = [ cellToIndex(1, yCell), 0, 0, 0 ]
      end do
      do xCell = xNumCells+1, xNumCells+xNumLayers 
        iCell = cellToIndex(xCell, yCell)
        mesh%mCellCenter(iCell, :) = [xDelta*(xCell+0.5_dp)&
                              , yDelta*(yCell+0.5_dp), 0.0_dp]
        mesh%CellToCell(iCell, :) &
          = [ 0, cellToIndex(xNumCells, yCell), 0, 0 ]
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(iCell, yCell, xCell)
    do xCell = 1, xNumCells
      do yCell = 1-yNumLayers, 0
        iCell = cellToIndex(xCell, yCell)
        mesh%mCellCenter(iCell, :) = [xDelta*(xCell+0.5_dp)&
                              ,yDelta*(yCell+0.5_dp), 0.0_dp]
        mesh%CellToCell(iCell, :) &
          = [ 0, 0, cellToIndex(xCell, 1), 0 ]
      end do
      do yCell = yNumCells+1, yNumCells+yNumLayers
        iCell = cellToIndex(xCell, yCell)
        mesh%mCellCenter(iCell, :) = [xDelta*(xCell+0.5_dp)&
                              ,yDelta*(yCell+0.5_dp), 0.0_dp]
        mesh%CellToCell(iCell, :) &
          = [ 0, 0, 0, cellToIndex(xCell, yNumCells) ]
      end do
    end do
    !$omp end parallel do
  end block
end subroutine tMesh_InitRect

end module StormRuler_Mesh
