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
module StormRuler_Matrix

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Helpers, only: I2S, R2S, PixelToInt
use StormRuler_Mesh, only: tMesh
use StormRuler_KrylovSolvers, only: @{tMatVecFunc$$@|@0, NUM_RANKS}@

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! ELLPACK sparse square matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tSparseMatrix
  ! ----------------------
  ! Number of matrix rows/columns.
  ! ----------------------
  integer(ip) :: NumRows
  ! ----------------------
  ! Number of row entries.
  ! ----------------------
  integer(ip) :: NumRowColumns
  ! ----------------------
  ! Column numbers per row.
  ! Shape is [1, NumRowColumns]×[1, NumRows].
  ! ----------------------
  integer(ip), allocatable :: RowColumnIndices(:,:)
  ! ----------------------
  ! Column values per row.
  ! Shape is [1, NumRowColumns]×[1, NumRows].
  ! ----------------------
  real(dp), allocatable :: RowColumnValues(:,:)
end type tSparseMatrix

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Matrix_SaveTo_Image2D(mat, image)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tSparseMatrix), intent(in) :: mat
  integer(ip), intent(out), allocatable :: image(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: iRow, iRowColumn

  ! ----------------------
  ! Allocate an image and fill it.
  ! ----------------------
  allocate( image(mat%NumRows, mat%NumRows) )
  image(:,:) = PixelToInt([255, 255, 255])
  do iRow = 1, mat%NumRows
    do iRowColumn = 1, mat%NumRowColumns
      associate(iColumn => mat%RowColumnIndices(iRowColumn, iRow))
        image(iRow, iColumn) = PixelToInt([0, 0, 0])
      end associate
    end do
  end do
end subroutine Matrix_SaveTo_Image2D

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Matrix_SaveTo_MatrixMarket(mat, file)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tSparseMatrix), intent(in) :: mat
  character(len=*), intent(in) :: file
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: unit
  integer(ip) :: iRow, iRowColumn
  
  open(newunit=unit, file=file, status='replace')

  write(unit, "(A)") &
    & '%%MatrixMarket matrix coordinate real general'
  write(unit, "(A, ' ', A, ' ', A)") &
    & I2S(mat%NumRows), I2S(mat%NumRows), I2S(mat%NumRows*mat%NumRowColumns)
  do iRow = 1_ip, mat%NumRows
    do iRowColumn = 1_ip, mat%NumRowColumns
      associate(iColumn => mat%RowColumnIndices(iRowColumn, iRow), &
        &   columnValue => mat%RowColumnValues(iRowColumn, iRow))
        write(unit, "(A, ' ', A, ' ', A)") &
          & I2S(iRow), I2S(iColumn), R2S(columnValue)
      end associate
    end do
  end do

  close(unit)
end subroutine Matrix_SaveTo_MatrixMarket

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Reconstruct a matrix using a matrix-vector product function,
!! using a very basic approach.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, 0*NUM_RANKS
subroutine MakeMatrix_Basic$rank(mesh, mat, stencilHalfWidth, symmetric, MatVec, env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  class(tSparseMatrix), intent(out) :: mat
  integer(ip), intent(in) :: stencilHalfWidth
  logical, intent(in) :: symmetric
  procedure(tMatVecFunc$rank) :: MatVec
  class(*), intent(in) :: env
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: iCell, iCellFace, iiCell
  integer(ip) :: iRow, iRowColumn, iColumn

  real(dp), allocatable :: u(@:,:), Au(@:,:)
  allocate( u(mesh%NumAllCells), Au(mesh%NumAllCells) )

  ! ----------------------
  ! Build a matrix portrait based on the stencil width.
  ! ----------------------
  mat%NumRows = mesh%NumCells
  mat%NumRowColumns = 1 + stencilHalfWidth*mesh%NumCellFaces
  allocate( mat%RowColumnIndices(mat%NumRowColumns, mat%NumRows) )
  mat%RowColumnIndices(:,:) = 0
  ! ----------------------
  !#omp parallel do schedule(static) &
  !#omp & default(private) shared(mesh, mat, stencilHalfWidth) if(.false.)
  do iCell = 1_ip, mesh%NumCells
    mat%RowColumnIndices(1_ip, iCell) = iCell
    do iCellFace = 1_ip, mesh%NumCellFaces
      iiCell = iCell
      do iRowColumn = &
        & 2_ip + stencilHalfWidth*(iCellFace-1), &
        & 2_ip + stencilHalfWidth*iCellFace
        iiCell = mesh%CellToCell(iCellFace, iiCell)
        mat%RowColumnIndices(iRowColumn, iCell) = iiCell
      end do
    end do
  end do
  !#omp end parallel do

  ! ----------------------
  ! Reconstruct matrix values.
  ! ----------------------
  allocate( mat%RowColumnValues(mat%NumRowColumns, mat%NumRows) )
  mat%RowColumnValues(:,:) = 0.0_dp
  do iCell = 1_ip, mesh%NumCells
    ! ----------------------
    ! Reconstruct row using matrix-vector product.
    ! ----------------------
    u(:) = 0.0_dp; u(iCell) = 1.0_dp
    call MatVec(mesh, Au, u, env)

    ! ----------------------
    ! Fill matrix row.
    ! ----------------------
    do iRow = 1, mat%NumRows
      do iRowColumn = 1, mat%NumRowColumns
        associate(iColumn => mat%RowColumnIndices(iRowColumn, iRow))
          if (iColumn == iCell) &
            mat%RowColumnValues(iRowColumn, iRow) = Au(iRow)
        end associate
      end do
    end do
  end do
end subroutine MakeMatrix_Basic$rank
#$end do

end module StormRuler_Matrix
