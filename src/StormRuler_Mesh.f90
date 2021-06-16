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

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Semi-structured 2D mesh.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
type :: Mesh2D
  ! ----------------------
  ! Number of interior cells.
  ! This value should be used for field operations.
  ! ----------------------
  integer :: NumCells
  ! ----------------------
  ! Total number of cells (including interior and ghost cells).
  ! This value should be used for field allocation.
  ! ----------------------
  integer :: NumAllCells

  ! ----------------------
  ! Number of faces (edges in 2D) per each cell.
  ! ----------------------
  integer :: NumCellFaces

  ! ----------------------
  ! Mesh dimension.
  ! ----------------------
  integer :: Dimensions
  ! ----------------------
  ! Multidimensional index bounds table.
  ! Shape is [1,Dimensions].
  ! ----------------------
  integer, allocatable :: MDIndexBounds(:)
  ! ----------------------
  ! Cell multidimensional index table.
  ! Shape is [1,Dimensions]×[1,NumCells].
  ! ----------------------
  integer, allocatable :: CellMDIndex(:,:)

  ! ----------------------
  ! Cell connectivity table.
  ! Shape is [-NumCellFaces/2,+NumCellFaces/2]×[1,NumAllCells].
  ! ----------------------
  integer, allocatable :: CellToCell(:,:), CellToCell1(:,:)

  ! ----------------------
  ! Number of boundary condition marks.
  ! ----------------------
  integer :: NumBCMs
  ! ----------------------
  ! BC mark to boundary cell index (in CSR format).
  ! ----------------------
  integer, allocatable :: BCMs(:)
  integer, allocatable :: BCMToCell(:)
  integer, allocatable :: BCMToCellFace(:)

  ! ----------------------
  ! Temporal step value.
  ! ----------------------
  real(dp) :: dt
  ! ----------------------
  ! Distance between centers of the adjacent cells per face.
  ! Shape is [-NumCellFaces/2,+NumCellFaces/2].
  ! ----------------------
  real(dp), allocatable :: dl(:)
  ! ----------------------
  ! Difference between centers of the adjacent cells per face.
  ! Shape is [1,Dimensions]×[-NumCellFaces/2,+NumCellFaces/2].
  ! ----------------------
  real(dp), allocatable :: dr(:,:)
  ! ----------------------
  ! Cell center coordinates.
  ! Shape is [1,Dimensions]×[1,NumCells].
  ! ----------------------
  real(dp), allocatable :: CellCenter(:,:)

contains
  procedure :: InitRect => Mesh2D_InitRect
end type Mesh2D
!! -----------------------------------------------------------------  

private Mesh2D_InitRect

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
!subroutine Mesh2D_Init(mesh, &
!                       xNumCells,yNumCells)
!  ! <<<<<<<<<<<<<<<<<<<<<<
!  class(Mesh2D), intent(inout) :: mesh
!  ! >>>>>>>>>>>>>>>>>>>>>>
!end subroutine Mesh2D_Init
!! -----------------------------------------------------------------  

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
