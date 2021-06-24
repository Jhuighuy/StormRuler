!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module StormRuler_Mesh

use StormRuler_Helpers
use StormRuler_IO

use, intrinsic :: iso_fortran_env, only: error_unit

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Semi-structured 2D mesh.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
type :: Mesh2D
  ! ----------------------
  ! Mesh dimension.
  ! ----------------------
  integer :: Dim

  ! ----------------------
  ! Number of faces (edges in 2D) per each cell.
  ! ----------------------
  integer :: NumCellFaces

  ! ----------------------
  ! Number of interior cells.
  ! This value should be used for field operations.
  ! ----------------------
  integer :: NumCells
  ! ----------------------
  ! Number of boundary cells.
  ! Cells may be counted multiple times..
  ! ----------------------
  integer :: NumBCCells
  ! ----------------------
  ! Total number of cells (including interior and ghost cells).
  ! This value should be used for field allocation.
  ! ----------------------
  integer :: NumAllCells
  ! ----------------------
  ! Cell connectivity table.
  ! Shape is [-NumFaces/2,+NumFaces/2]×[1,NumAllCells].
  ! ----------------------
  integer, allocatable :: CellToCell(:,:), CellToCell1(:,:)

  ! ----------------------
  ! Multidimensional index bounds table.
  ! Shape is [1,Dim].
  ! ----------------------
  integer, allocatable :: MDIndexBounds(:)
  ! ----------------------
  ! Cell multidimensional index table.
  ! Shape is [1,Dim]×[1,NumCells].
  ! ----------------------
  integer, allocatable :: CellMDIndex(:,:)

  ! ----------------------
  ! Number of boundary condition marks.
  ! ----------------------
  integer :: NumBCMs
  ! ----------------------
  ! BC mark to boundary cell index (in CSR format).
  ! Shape is [1,NumBCMs].
  ! ----------------------
  integer, allocatable :: BCMs(:)
  ! ----------------------
  ! Shape is [1,NumBCCells].
  ! ----------------------
  integer, allocatable :: BCMToCell(:)
  ! ----------------------
  ! Shape is [1,NumBCCells].
  ! ----------------------
  integer, allocatable :: BCMToCellFace(:)

  ! ----------------------
  ! Temporal step value.
  ! ----------------------
  real(dp) :: dt
  ! ----------------------
  ! Distance between centers of the adjacent cells per face.
  ! Shape is [-NumFaces/2,+NumFaces/2].
  ! ----------------------
  real(dp), allocatable :: dl(:)
  ! ----------------------
  ! Difference between centers of the adjacent cells per face.
  ! Shape is [1,Dim]×[-NumFaces/2,+NumFaces/2].
  ! ----------------------
  real(dp), allocatable :: dr(:,:)
  ! ----------------------
  ! Cell center coordinates.
  ! Shape is [1,Dim]×[1,NumCells].
  ! ----------------------
  real(dp), allocatable :: CellCenter(:,:)

contains
  procedure :: InitRect => Mesh2D_InitRect
  procedure :: InitFromImage2D => Mesh2D_InitFromImage2D
  procedure :: InitFromImage3D => Mesh2D_InitFromImage3D
  ! ----------------------
  procedure :: PrintCellToCell_Neato => Mesh2D_PrintCellToCell_Neato
end type Mesh2D
!! -----------------------------------------------------------------  

private Mesh2D_InitRect


!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

contains

#$let STENCIL = { &
  & 2:{ +1:[+1,0],-1:[-1,0],   &
  &     +2:[0,+1],-2:[0,-1] }, &
  & 3:{ +1:[+1,0,0],-1:[-1,0,0], &
  &     +2:[0,+1,0],-2:[0,-1,0], &
  &     +3:[0,0,+1],-3:[0,0,-1] } }

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize a mesh from image.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do dim = 2, 3
subroutine Mesh2D_InitFromImage${dim}$D(mesh,image,fluidColor,colorToBCM,numBCLayers)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(inout) :: mesh
  integer, intent(in) :: image(@{0:}@),fluidColor,colorToBCM(:),numBCLayers
  ! >>>>>>>>>>>>>>>>>>>>>>
#$define iCellMD @{iCellMD$$}@
  integer :: iCell,iBCCell,iGCell,$iCellMD,iCellFace,iBCM
  integer, allocatable :: cache(@:)
  allocate(cache,mold=image); cache(@:) = 0
  ! ----------------------
  ! Process image data in order to:
  ! 1. Generate initial "Cell MD Index" ⟺ "Cell (Plain Index)" table 
  !    for interior cells and compute number of interior and boundary cells.
  ! 2. Generate and pre-fill "BCM" ⟺ "Num BC Cells per BCM" CSR-style table.
  ! ----------------------
  mesh%Dim = $dim
  mesh%NumCellFaces = ${2*dim}$
  mesh%MDIndexBounds = shape(image)-2
  mesh%NumBCMs = size(colorToBCM,dim=1)
  allocate(mesh%BCMs(mesh%NumBCMs+1)); mesh%BCMs(:) = 0
  ! ----------------------
  iCell = 0; mesh%BCMs(:) = 0
#$do rank = dim, 1, -1
  do iCellMD$rank = 1, mesh%MDIndexBounds($rank) 
#$end do
    if (image($iCellMD) == fluidColor) then
      iCell = iCell + 1
      cache($iCellMD) = iCell
#$for iCellFace,iCellFaceMD in STENCIL[dim]
#$define iiCellMD @{iCellMD$$+${iCellFaceMD[$$-1]}$}@
      if (image($iiCellMD) /= fluidColor) then
        iBCM = Find(colorToBCM,image($iiCellMD))
        if (iBCM == 0) then
          write(error_unit,*) 'INVALID IMAGE COLOR', IntToPixel(image($iiCellMD))
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
  associate(nac=>mesh%NumAllCells,ncf=>mesh%NumCellFaces)
    allocate(mesh%CellMDIndex(mesh%Dim,nac))
    allocate(mesh%CellToCell((-ncf/2):(+ncf/2),nac)); mesh%CellToCell(:,:) = 0
  end associate
  ! ----------------------
  iBCCell = mesh%NumCells + 1
#$do rank = dim, 1, -1
  do iCellMD$rank = 1, mesh%MDIndexBounds($rank) 
#$end do
    if (image($iCellMD) == fluidColor) then
      iCell = cache($iCellMD)
      mesh%CellToCell(0,iCell) = iCell
      mesh%CellMDIndex(:,iCell) = [$iCellMD]
#$for iCellFace,iCellFaceMD in STENCIL[dim]
#$define iiCellMD @{iCellMD$$+${iCellFaceMD[$$-1]}$}@
      if (image($iiCellMD) == fluidColor) then
        mesh%CellToCell($iCellFace,iCell) = cache($iiCellMD)
      else
        iBCM = Find(colorToBCM,image($iiCellMD))
        mesh%CellToCell($iCellFace,iCell) = -iBCM
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
  ! 5. Allocate and fill the "BCM" ⟺ "BC Cell" 
  !    and "BCM" ⟺ "BC Cell Face" CSR-style tables and fill it. 
  ! ----------------------
  associate(nbcc=>mesh%NumBCCells)
    allocate(mesh%BCMToCell(nbcc))
    allocate(mesh%BCMToCellFace(nbcc))
  end associate
  ! ----------------------
  iBCCell = mesh%NumCells + 1
  do iCell = 1, mesh%NumCells
    do iCellFace = -mesh%NumCellFaces/2, +mesh%NumCellFaces/2
      iBCM = -mesh%CellToCell(+iCellFace,iCell) 
      if (iBCM > 0) then
        ! Classic CSR insertion procedure here.
        associate(ptr=>mesh%BCMs(iBCM))
          mesh%BCMToCell(ptr+1) = iBCCell
          mesh%BCMToCellFace(ptr+1) = iCellFace; ptr = ptr + 1
        end associate
        ! Generate connectivity for BC cell and ghost cells.
        mesh%CellToCell(+iCellFace,iCell) = iBCCell
        do iGCell = iBCCell, iBCCell+numBCLayers-1
          mesh%CellToCell(0,iGCell) = iGCell
          mesh%CellToCell(+iCellFace,iGCell) = iGCell+1
          mesh%CellToCell(-iCellFace,iGCell) = iGCell-1
          ! Fill MD index information by linear extrapolation.
          associate(iCellMD=>mesh%CellMDIndex(:,iCell), &
            &      iiCellMD=>mesh%CellMDIndex(:,mesh%CellToCell(-iCellFace,iCell)))
            mesh%CellMDIndex(:,iGCell) = &
              & iCellMD + (iGCell-iBCCell+1)*(iCellMD-iiCellMD)
          end associate
        end do
        mesh%CellToCell(-iCellFace,iBCCell) = iCell
        mesh%CellToCell(+iCellFace,iBCCell+numBCLayers-1) = 0
        ! Increment BC cell index.
        iBCCell = iBCCell+numBCLayers
      end if
    end do
  end do
  ! ----------------------
  ! Classic CSR pointer correction procedure.
  mesh%BCMs(:) = eoshift(mesh%BCMs(:),shift=-1) + 1
  print *, 'ERROR PRONE CHECK:', iBCCell-1
  ! ----------------------
end subroutine Mesh2D_InitFromImage${dim}$D
#$end do

#$del STENCIL

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh connectivity in '.dot' format.
!! Output may be visualized with: 'neato -n -Tpng c2c.dot > c2c.png' 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mesh2D_PrintCellToCell_Neato(mesh,file)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(inout) :: mesh
  character(len=*), intent(in) :: file
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: unit
  integer :: iCell,iiCell,iCellFace
  integer :: iBCM,iBCMPtr,iBCCell,iiBCCell,iBCCellFace
  ! ----------------------
  integer, parameter :: DPI = 72
  character(len=10), parameter :: SHAPES(2) = &
    & [character(len=10) :: 'box', 'diamond']
  character(len=10), parameter :: PALETTE(3) = &
    & [character(len=10) :: 'red', 'green', 'blue']
  ! ----------------------
  if (mesh%Dim /= 2) then
    error stop 'Only 2D meshes can be printed to Neato'
  end if
  ! ----------------------
  open(newunit=unit,file=file,status='replace')
  write(unit,'(A)') '# Generated by StromRuler/Mesh2Neato'
  write(unit,'(A)') 'digraph CellToCell {'
  ! ----------------------
  ! Dump cell MD indices of the cells as graph nodes.
  ! Difference node shapes are use for different face 
  ! directions of the BC/ghost cells faces.
  ! ----------------------
100 format('  ',A,'[label=',A,' ','pos=',A,']')
200 format('  ',A,'[label=',A,' pos=',A,' shape=',A,']')
  do iCell = 1, mesh%NumCells
    associate(iCellPos=>DPI*mesh%CellMDIndex(:,iCell))
      write(unit,100) 'C'//I2S(iCell), &
        & '"C"', '"'//I2S(iCellPos(1))//','//I2S(iCellPos(2))//'!"'
    end associate
  end do
  do iBCM = 1, mesh%NumBCMs
    do iBCMPtr = mesh%BCMs(iBCM), mesh%BCMs(iBCM+1)-1
      iBCCell = mesh%BCMToCell(iBCMPtr)
      iBCCellFace = mesh%BCMToCellFace(iBCMPtr)
      do while(iBCCell /= 0)
        associate(iBCCellPos=>DPI*mesh%CellMDIndex(:,iBCCell))
          write(unit,200) 'C'//I2S(iBCCell), '"B'//I2S(iBCM)//'"', &
            & '"'//I2S(iBCCellPos(1))//','//I2S(iBCCellPos(2))//'!"', &
            & trim(SHAPES(abs(iBCCellFace)))
        end associate
        iBCCell = mesh%CellToCell(iBCCellFace,iBCCell)
      end do
    end do
  end do
  ! ----------------------
  ! Dump cell to cell connectivity as graph edges.
  ! BC/ghost cells are colored by the BC mark.
  ! ----------------------
300 format('  ',A,'->',A)
400 format('  ',A,'->',A,' [color=',A,']')
  do iCell = 1, mesh%NumCells
    do iCellFace = -mesh%NumCellFaces/2, +mesh%NumCellFaces/2
      if (iCellFace /= 0) then
        iiCell = mesh%CellToCell(iCellFace,iCell)
        write(unit,300) 'C'//I2S(iCell), 'C'//I2S(iiCell)
      end if
    end do
  end do
  do iBCM = 1, mesh%NumBCMs
    do iBCMPtr = mesh%BCMs(iBCM), mesh%BCMs(iBCM+1)-1
      iBCCell = mesh%BCMToCell(iBCMPtr)
      iBCCellFace = mesh%BCMToCellFace(iBCMPtr)
      do
        associate(color=>PALETTE(iBCM))
          iiBCCell = mesh%CellToCell(-iBCCellFace,iBCCell) 
          write(unit,400) 'C'//I2S(iBCCell), 'C'//I2S(iiBCCell), trim(color)
          iiBCCell = mesh%CellToCell(+iBCCellFace,iBCCell)
          if (iiBCCell == 0) exit
          write(unit,400) 'C'//I2S(iBCCell), 'C'//I2S(iiBCCell), trim(color)
        end associate
        iBCCell = iiBCCell
      end do
    end do
  end do
  ! ----------------------
  ! Close file and exit.
  write(unit,'(A)') '}'
  close(unit)
end subroutine Mesh2D_PrintCellToCell_Neato

!! -----------------------------------------------------------------  
!! Initialize a 2D rectangular mesh.
subroutine Mesh2D_InitRect(mesh,xDelta,xNumCells,xPeriodic &
                               ,yDelta,yNumCells,yPeriodic,numLayers)
  class(Mesh2D), intent(inout) :: mesh
  real(dp), intent(in) :: xDelta,yDelta
  integer, intent(in) :: xNumCells,yNumCells
  logical, intent(in) :: xPeriodic,yPeriodic
  integer, intent(in), optional :: numLayers
  integer :: xNumLayers,yNumLayers
  integer, allocatable :: cellToIndex(:,:)
  ! ----------------------
  ! Set up amount of layers.
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
    integer :: xCell,yCell
    allocate(cellToIndex(1-xNumLayers:xNumCells+xNumLayers &
                       , 1-yNumLayers:yNumCells+yNumLayers))
    ! ----------------------
    ! Generate indices for the interior cells.
    associate(iCell => mesh%NumCells)
      iCell = 0
      do yCell = 1, yNumCells
        do xCell = 1, xNumCells
          iCell = iCell + 1
          cellToIndex(xCell,yCell) = iCell
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
          cellToIndex(0,yCell) = cellToIndex(xNumCells,yCell)
          cellToIndex(xNumCells+1,yCell) = cellToIndex(1,yCell)
        end do
        !$omp end parallel do
      else
        do yCell = 1, yNumCells
          do xCell = 1-xNumLayers, 0
            iCell = iCell + 1
            cellToIndex(xCell,yCell) = iCell
          end do
          do xCell = xNumCells+1, xNumCells+xNumLayers
            iCell = iCell + 1
            cellToIndex(xCell,yCell) = iCell
          end do
        end do
      end if
      if (yPeriodic) then
        yNumLayers = 0
        !$omp parallel do private(xCell)
        do xCell = 1, xNumCells
          cellToIndex(xCell,0) = cellToIndex(xCell,yNumCells)
          cellToIndex(xCell,yNumCells+1) = cellToIndex(xCell,1)
        end do
        !$omp end parallel do
      else
        do xCell = 1, xNumCells
          do yCell = 1-yNumLayers, 0
            iCell = iCell + 1
            cellToIndex(xCell,yCell) = iCell
          end do
          do yCell = yNumCells+1, yNumCells+yNumLayers
            iCell = iCell + 1
            cellToIndex(xCell,yCell) = iCell
          end do
        end do
      end if
    end associate
  end block
  ! ----------------------
  ! Fill the cell geometry and connectivity.
  block
    integer :: iCell, xCell,yCell
    allocate(mesh%dl(1:4))
    mesh%dl(:) = [xDelta,xDelta,yDelta,yDelta]
    allocate(mesh%dr(1:2,1:4))
    mesh%dr(:,1) = [mesh%dl(1),0.0_dp]
    mesh%dr(:,2) = [0.0_dp,mesh%dl(2)]
    mesh%dr(:,3) = [mesh%dl(1),0.0_dp]
    mesh%dr(:,4) = [0.0_dp,mesh%dl(2)]
    mesh%NumCellFaces = 4
    allocate(mesh%CellCenter(1:mesh%NumAllCells,1:3))
    allocate(mesh%CellToCell(1:mesh%NumAllCells,1:4))
    allocate(mesh%CellToCell1(-2:2,mesh%NumAllCells))
    associate(cellCenter => mesh%CellCenter &
            , cellToCell => mesh%CellToCell)
      ! ----------------------
      ! Fill the cell information for the interior cells.
      !$omp parallel do private(iCell,yCell,xCell) collapse(2)
      do yCell = 1, yNumCells
        do xCell = 1, xNumCells
          iCell = cellToIndex(xCell,yCell)
          cellCenter(iCell,:) = [xDelta*(xCell+0.5_dp)&
                                ,yDelta*(yCell+0.5_dp),0.0_dp]
          cellToCell(iCell,:) &
            = [ cellToIndex(xCell+1,yCell), cellToIndex(xCell-1,yCell) &
              , cellToIndex(xCell,yCell+1), cellToIndex(xCell,yCell-1) ]
          !---
          mesh%CellToCell1(:,iCell) = &
            & [cellToIndex(xCell,yCell-1), &
            &  cellToIndex(xCell-1,yCell), &
            &  cellToIndex(xCell , yCell), &
            &  cellToIndex(xCell+1,yCell), &
            &  cellToIndex(xCell,yCell+1)]
          end do
      end do
      !$omp end parallel do
      ! ----------------------
      ! Fill the cell information for the ghost cells.
      ! (Remember, the the ghost cells have a single connection to the interior cell.)
      !$omp parallel do private(iCell,yCell,xCell)
      do yCell = 1, yNumCells
        do xCell = 1-xNumLayers, 0
          iCell = cellToIndex(xCell,yCell)
          cellCenter(iCell,:) = [xDelta*(xCell+0.5_dp)&
                                ,yDelta*(yCell+0.5_dp),0.0_dp]
          cellToCell(iCell,:) &
            = [ cellToIndex(1,yCell), 0, 0, 0 ]
        end do
        do xCell = xNumCells+1, xNumCells+xNumLayers
          iCell = cellToIndex(xCell,yCell)
          cellCenter(iCell,:) = [xDelta*(xCell+0.5_dp)&
                               , yDelta*(yCell+0.5_dp),0.0_dp]
          cellToCell(iCell,:) &
            = [ 0, cellToIndex(xNumCells,yCell), 0, 0 ]
        end do
      end do
      !$omp end parallel do
      !$omp parallel do private(iCell,yCell,xCell)
      do xCell = 1, xNumCells
        do yCell = 1-yNumLayers, 0
          iCell = cellToIndex(xCell,yCell)
          cellCenter(iCell,:) = [xDelta*(xCell+0.5_dp)&
                                ,yDelta*(yCell+0.5_dp),0.0_dp]
          cellToCell(iCell,:) &
            = [ 0, 0, cellToIndex(xCell,1), 0 ]
        end do
        do yCell = yNumCells+1, yNumCells+yNumLayers
          iCell = cellToIndex(xCell,yCell)
          cellCenter(iCell,:) = [xDelta*(xCell+0.5_dp)&
                                ,yDelta*(yCell+0.5_dp),0.0_dp]
          cellToCell(iCell,:) &
            = [ 0, 0, 0, cellToIndex(xCell,yNumCells) ]
        end do
      end do
      !$omp end parallel do
    end associate
  end block
end subroutine Mesh2D_InitRect

end module StormRuler_Mesh
