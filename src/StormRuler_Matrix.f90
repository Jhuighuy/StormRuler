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
use StormRuler_Helpers!, only: 

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
type :: tMatrix

  ! ----------------------
  ! Row nonzeroes addresses.
  ! Shape is [1, NumCells + 1].
  ! ----------------------
  integer(ip), allocatable :: RowAddrs(:)

  ! ----------------------
  ! Row nonzeroes column indices.
  ! Shape is [1, RowAddrs(NumCells + 1) - 1].
  ! ----------------------
  integer(ip), allocatable :: ColIndices(:)

  ! ----------------------
  ! Row nonzeroes column coefficients.
  ! Shape is [1, NumVars]×[1, NumVars]×[1, RowAddrs(NumCells + 1) - 1].
  ! ----------------------
  real(dp), allocatable :: ColCoeffs(:,:,:)

end type tMatrix

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize a basic sparse matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitMatrix(mesh, mat, power)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(inout) :: mat
  integer(ip), intent(in) :: power

  type :: tMatrixRow
    integer(ip), allocatable :: ColIndices(:)
  end type tMatrixRow

  integer(ip) :: k, row, col, colCol, colAddr, colColAddr
  type(tMatrixRow), allocatable :: matRows(:)

  ! ----------------------
  ! Generate matrix rows.
  ! ----------------------
  allocate(matRows(mesh%NumCells))
  do row = 1, mesh%NumCells

    col = row
    matRows(row)%ColIndices = [col]

    do k = 1, power
      do colAddr = 1, size(matRows(row)%ColIndices)
        col = matRows(row)%ColIndices(colAddr)

        do colColAddr = 1, mesh%NumCellFaces
          colCol = mesh%CellToCell(colColAddr,col)

          if ((colCol <= mesh%NumCells).and.&
              & (IndexOf(colCol, matRows(row)%ColIndices) == 0)) then
            matRows(row)%ColIndices = [matRows(row)%ColIndices, colCol]
          end if
          
        end do
      end do
    end do

    call BubbleSort(matRows(row)%ColIndices)

  end do

  ! ----------------------
  ! Compress the rows.
  ! ----------------------
  allocate(mat%RowAddrs(mesh%NumCells + 1))
  mat%RowAddrs(1) = 1
  do row = 1, mesh%NumCells
    mat%RowAddrs(row + 1) = mat%RowAddrs(row) + size(matRows(row)%ColIndices)
  end do

  associate(nnz => mat%RowAddrs(mesh%NumCells + 1) - 1)
    allocate(mat%ColIndices(nnz), mat%ColCoeffs(1, 1, nnz))
  end associate
  do row = 1, mesh%NumCells
    associate(first => mat%RowAddrs(row), last => mat%RowAddrs(row + 1) - 1)

      mat%ColIndices(first:last) = matRows(row)%ColIndices(:)
      mat%ColCoeffs(:,:,first:last) = 0.0_dp

    end associate
  end do

end subroutine InitMatrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Sparse matrix-vector product: 𝒚 ← 𝓐𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine MatrixVector(mesh, mat, yArr, xArr)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tArray), intent(inout) :: xArr, yArr

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(MatrixVector_Kernel)

contains
  subroutine MatrixVector_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr

    y(:,row) = 0.0_dp

    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      associate(col => mat%ColIndices(rowAddr), &
           & coeffs => mat%ColCoeffs(:,:,rowAddr))

        y(:,row) = y(:,row) + matmul(coeffs, x(:,col))

      end associate
    end do

  end subroutine MatrixVector_Kernel
end subroutine MatrixVector

end module StormRuler_Matrix
