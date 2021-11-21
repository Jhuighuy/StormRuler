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
!! Compressed sparse row (CSR) block matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tRowMatrix

  ! ----------------------
  ! Row nonzeroes pointers.
  ! Shape is [1, NumCells + 1].
  ! ----------------------
  integer(ip), allocatable :: RowPtrs(:)

  ! ----------------------
  ! Row nonzeroes column indices.
  ! Shape is [1, RowPtrs(NumCells + 1) - 1].
  ! ----------------------
  integer(ip), allocatable :: ColumnIndices(:)

  ! ----------------------
  ! Row nonzeroes column coefficients.
  ! Shape is [1, NumVars]×[1, NumVars]×[1, RowPtrs(NumCells + 1) - 1].
  ! ----------------------
  real(dp), allocatable :: ColumnCoeffs(:,:,:)

end type tRowMatrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compressed sparse column (CSC) block matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tColumnMatrix

  ! ----------------------
  ! Column nonzeroes pointers.
  ! Shape is [1, NumCells + 1].
  ! ----------------------
  integer(ip), allocatable :: ColumnPtrs(:)

  ! ----------------------
  ! Column nonzeroes row indices.
  ! Shape is [1, ColumnPtrs(NumCells + 1) - 1].
  ! ----------------------
  integer(ip), allocatable :: RowIndices(:)

  ! ----------------------
  ! Column nonzeroes row coefficients.
  ! Shape is [1, NumVars]×[1, NumVars]×[1, RowPtrs(NumCells + 1) - 1].
  ! ----------------------
  real(dp), allocatable :: RowCoeffs(:,:,:)

end type tColumnMatrix

interface SparseMatVec
  module procedure RowMatVec
  module procedure ColumnMatVec
end interface SparseMatVec

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize compressed sparse column matrix 
!! with the specified power of mesh bandwidth.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitBandedColumnMatrix(matrix, mesh, power)
  class(tColumnMatrix), intent(inout) :: matrix
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: power

  integer(ip) :: cell, cellFace, row, column
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

  halfBandwidth = power*halfBandwidth
  bandwidth = 2*halfBandwidth + 1

  ! ----------------------
  ! Fill the columns pointers and row entries.
  ! ----------------------
  allocate(matrix%ColumnPtrs(mesh%NumCells + 1))
  !! TODO: size is overestimated here.
  !! (nc + 1)*nc/2 - ...
  associate(size => mesh%NumCells*bandwidth)
    allocate(matrix%RowIndices(size))
    allocate(matrix%RowCoeffs(1,1,size))
  end associate

  matrix%ColumnPtrs(1) = 1
  do column = 1, mesh%NumCells

    associate(columnPtr => matrix%ColumnPtrs(column + 1))
      columnPtr = matrix%ColumnPtrs(column)
      do row = max(column - halfBandwidth, 1), &
             & min(column + halfBandwidth, mesh%NumCells)
        matrix%RowIndices(columnPtr) = row
        matrix%RowCoeffs(:,:,columnPtr) = 0.0_dp
        columnPtr = columnPtr + 1
      end do
    end associate

  end do

end subroutine InitBandedColumnMatrix

subroutine ColumnToRowMatrix(mesh, rowMat, colMat)
  class(tMesh), intent(inout) :: mesh
  class(tRowMatrix), intent(inout) :: rowMat
  class(tColumnMatrix), intent(in) :: colMat

  integer(ip) :: row, column, columnPtr, index, numEntries

  allocate(rowMat%RowPtrs, mold=colMat%ColumnPtrs)
  allocate(rowMat%ColumnIndices, mold=colMat%RowIndices)
  allocate(rowMat%ColumnCoeffs, mold=colMat%RowCoeffs)

  numEntries = colMat%ColumnPtrs(size(colMat%ColumnPtrs)) - 1
  
  ! ----------------------
  ! Fill row pointers.
  ! ----------------------
  rowMat%RowPtrs(:) = 0
  do index = 1, numEntries

    associate(rowPtr => rowMat%RowPtrs(colMat%RowIndices(index) + 1))
      rowPtr = rowPtr + 1
    end associate
  
  end do

  rowMat%RowPtrs(1) = 1
  do row = 1, mesh%NumCells
  
    rowMat%RowPtrs(row + 1) = rowMat%RowPtrs(row) + rowMat%RowPtrs(row + 1)
  
  end do

  ! ----------------------
  ! Fill column entries.
  ! ----------------------
  do column = 1, mesh%NumCells
    do columnPtr = colMat%ColumnPtrs(column), colMat%ColumnPtrs(column + 1) - 1

      associate(row => colMat%RowIndices(columnPtr), &
           & coeffs => colMat%RowCoeffs(:,:,columnPtr))
        associate(rowPtr => rowMat%RowPtrs(row))
          rowMat%ColumnIndices(rowPtr) = column
          rowMat%ColumnCoeffs(:,:,rowPtr) = coeffs
          rowPtr = rowPtr + 1
        end associate
      end associate

    end do
  end do

  rowMat%RowPtrs = eoshift(rowMat%RowPtrs, shift=-1)
  rowMat%RowPtrs(1) = 1

end subroutine ColumnToRowMatrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compressed sparse column matrix-vector product. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine RowMatVec(mesh, mat, yArr, xArr)
  class(tMesh), intent(inout) :: mesh
  class(tRowMatrix), intent(in) :: mat
  class(tArray), intent(inout) :: xArr, yArr

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(RowMatVec_Kernel)

contains
  subroutine RowMatVec_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowPtr

    y(:,row) = 0.0_dp

    do rowPtr = mat%RowPtrs(row), mat%RowPtrs(row + 1) - 1

      associate(column => mat%ColumnIndices(rowPtr), &
              & coeffs => mat%ColumnCoeffs(:,:,rowPtr))
        y(:,row) = y(:,row) + matmul(coeffs, x(:,column))
      end associate

    end do

  end subroutine RowMatVec_Kernel
end subroutine RowMatVec

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compressed sparse column matrix-vector product. 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ColumnMatVec(mesh, mat, yArr, xArr)
  class(tMesh), intent(inout) :: mesh
  class(tColumnMatrix), intent(in) :: mat
  class(tArray), intent(inout) :: xArr, yArr

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)
  
  call Fill(mesh, yArr, 0.0_dp)
  call mesh%RunCellKernel(ColumnMatVec_Kernel)

contains
  subroutine ColumnMatVec_Kernel(column)
    integer(ip), intent(in) :: column

    integer(ip) :: columnPtr

    do columnPtr = mat%ColumnPtrs(column), mat%ColumnPtrs(column + 1) - 1

      associate(row => mat%RowIndices(columnPtr), &
           & coeffs => mat%RowCoeffs(:,:,columnPtr))
        y(:,row) = y(:,row) + matmul(coeffs, x(:,column))
      end associate

    end do

  end subroutine ColumnMatVec_Kernel
end subroutine ColumnMatVec

end module StormRuler_Matrix
