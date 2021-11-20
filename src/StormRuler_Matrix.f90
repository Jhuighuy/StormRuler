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

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: Fill
use StormRuler_BLAS, only: tMatVecFunc

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Marix column coloring.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tColumnColoring
  ! ----------------------
  ! Number of the colors.
  ! ----------------------
  integer(ip) :: NumColors

  ! ----------------------
  ! Shape is [1, Cells].
  ! ----------------------
  integer(ip), allocatable :: ColumnColor(:)

  ! ----------------------
  ! Shape is [1, NumColors + 1].
  ! ----------------------
  integer(ip), allocatable :: ColorColumns(:)
  ! ----------------------
  ! Shape is [1, Cells].
  ! ----------------------
  integer(ip), allocatable :: ColorColumnIndices(:)
end type tColumnColoring

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compressed sparse column (CSC) matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tColumnMatrix
  ! ----------------------
  ! Shape is [1, NumCells + 1].
  ! ----------------------
  integer(ip), allocatable :: ColumnPtr(:)
  ! ----------------------
  ! Column nonzeroes indices.
  ! Shape is [1, ColumnPtrs(NumCells + 1)].
  ! ----------------------
  integer(ip), allocatable :: RowIndices(:)
  ! ----------------------
  ! Column nonzeroes values.
  ! Shape is [1, ColumnPtrs(NumCells + 1)].
  ! TODO: Should be [1, NumVars]×[1, NumVars]×Shape in block case.
  ! ----------------------
  real(dp), allocatable :: RowEntries(:)
end type tColumnMatrix

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

subroutine InitBandedColumnMatrix(matrix, mesh, power)
  class(tColumnMatrix), intent(inout) :: matrix
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: power

  integer(ip) :: cell, cellFace, column, row
  integer(ip) :: bandwidth, halfBandwidth

  ! ----------------------
  ! Compute bandwidth.
  ! ----------------------
  halfBandwidth = 0
  do cell = 1, mesh%NumCells
    do cellFace = 1, mesh%NumCellFaces
      associate(cellCell => mesh%CellToCell(cellFace, cell))

        if (cellCell <= mesh%NumCells) then
          halfBandwidth = max(halfBandwidth, abs(cell - cellCell))
        end if

      end associate
    end do
  end do

  halfBandwidth = halfBandwidth*power
  bandwidth = 2*halfBandwidth + 1

  ! ----------------------
  ! Fill the columns entries.
  ! ----------------------
  allocate(matrix%ColumnPtr(mesh%NumCells + 1))
  !! TODO: size is overestimated here.
  !! (nc + 1)*nc/2 - ...
  associate(size => mesh%NumCells*bandwidth)
    allocate(matrix%RowIndices(size))
    allocate(matrix%RowEntries(size))
  end associate

  matrix%ColumnPtr(1) = 1
  do column = 1, mesh%NumCells

    associate(columnPtr => matrix%ColumnPtr(column + 1))
      columnPtr = matrix%ColumnPtr(column)
      do row = max(column - halfBandwidth, 1), &
             & min(column + halfBandwidth, mesh%NumCells)
        matrix%RowIndices(columnPtr) = row
        matrix%RowEntries(columnPtr) = 0.0_dp
        columnPtr = columnPtr + 1
      end do
    end associate

  end do

end subroutine InitBandedColumnMatrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Generate column coloring using the banded algorithm.
!!
!! Estimated number of colors is ≈ 𝘣𝘸(𝓐ᵖ) = 𝑝⋅(𝘣𝘸(𝓐) - 1) + 1,
!! where 𝘣𝘸(𝓐) is the bandwidth of 𝓐.
!! 
!! Banded algorithm is simple and fast, but suboptimal in 
!! the most practical cases (for 𝑝 = 1 it generates ≈2 times
!! more colors than more sophisticated coloring algorithms).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ColorColumns_Banded(coloring, mesh, matrix)
  class(tColumnColoring), intent(inout) :: coloring
  class(tMesh), intent(in) :: mesh
  class(tColumnMatrix), intent(in) :: matrix

  integer(ip) :: cell, color
  integer(ip) :: bandwidth

  ! ----------------------
  ! Compute bandwidth.
  ! ----------------------
  bandwidth = 401
  coloring%NumColors = bandwidth

  ! ----------------------
  ! Fill the column to color table.
  ! ----------------------
  allocate(coloring%ColumnColor(mesh%NumCells))
  do cell = 1, mesh%NumCells

    color = mod(cell - 1, coloring%NumColors) + 1 
    coloring%ColumnColor(cell) = color

  end do

  ! ----------------------
  ! Fill the color to column table.
  ! ----------------------
  allocate(coloring%ColorColumns(coloring%NumColors + 1))
  allocate(coloring%ColorColumnIndices(mesh%NumCells))

  coloring%ColorColumns(1) = 1
  do color = 1, coloring%NumColors

    associate(colorColumn => coloring%ColorColumns(color + 1))
      colorColumn = coloring%ColorColumns(color)
      do cell = color, mesh%NumCells, bandwidth
        coloring%ColorColumnIndices(colorColumn) = cell 
        colorColumn = colorColumn + 1
      end do
    end associate

  end do

  print *, coloring%ColorColumns

end subroutine ColorColumns_Banded

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

subroutine ReconstructMatrix(matrix, coloring, mesh, MatVec, mold)
  class(tColumnMatrix), intent(inout) :: matrix
  class(tColumnColoring), intent(in) :: coloring
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: mold
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: cell, row, column, color
  real(dp), pointer :: tDat(:,:)
  type(tArray) :: Q, t, s

  call AllocArray(t, mold=mold)
  call AllocArray(Q, shape=[mold%mShape, coloring%NumColors])

  ! ----------------------
  ! Reconstruct the compressed columns.
  ! ----------------------
  do color = 1, coloring%NumColors
    call Fill(mesh, t, 0.0_dp)
    do column = coloring%ColorColumns(color), coloring%ColorColumns(color + 1) - 1
      call t%Get(tDat)
      tDat(:,coloring%ColorColumnIndices(column)) = 1.0_dp
    end do
    s = Q%Slice(color); call MatVec(mesh, s, t)
  end do

  ! ----------------------
  ! Decompress the columns.
  ! ----------------------
  do cell = 1, mesh%NumCells
    color = coloring%ColumnColor(cell)
    t = Q%Slice(color); call t%Get(tDat)
    do column = matrix%ColumnPtr(cell), matrix%ColumnPtr(cell + 1) - 1

      matrix%RowEntries(column) = tDat(1,matrix%RowIndices(column))

    end do
  end do

end subroutine ReconstructMatrix

subroutine ColumnMatVec(mesh, matrix, Ax, x)
  class(tColumnMatrix), intent(inout) :: matrix
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(inout) :: Ax, x

  integer(ip) :: cell, row, column, color

  real(dp) :: a
  real(dp), pointer :: y(:), z(:)

  call Ax%Get(y)
  call x%Get(z)

  y(:) = 0.0_dp
  do cell = 1, mesh%NumCells
    do column = matrix%ColumnPtr(cell), matrix%ColumnPtr(cell + 1) - 1
      row = matrix%RowIndices(column)
      a = matrix%RowEntries(column)

      y(row) = y(row) + a*z(cell)

    end do
  end do

end subroutine ColumnMatVec

end module StormRuler_Matrix
