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

include 'mkl_spblas.fi'

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

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Sparse matrix-vector product: 𝒚 ← 𝓐𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine MatrixVector(mesh, mat, yArr, xArr)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(MatrixVector_Kernel)

contains
  subroutine MatrixVector_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col

    y(:,row) = 0.0_dp

    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      col = mat%ColIndices(rowAddr)
      y(:,row) = y(:,row) + matmul(mat%ColCoeffs(:,:,rowAddr), x(:,col))
    end do

  end subroutine MatrixVector_Kernel
end subroutine MatrixVector

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Sparse partial matrix-vector product: 𝒚 ← 𝓓𝒙, 𝒚 ← 𝓛𝒙 or 𝒚 ← 𝓤𝒙,
!! where 𝓓 is the (block-)diagonal of 𝓐, 𝓛 and 𝓤 are lower and upper 
!! strict (block-)triangular parts of 𝓐, 𝓐 = 𝓛 + 𝓓 + 𝓤.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine PartialMatrixVector(mesh, part, mat, yArr, xArr)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  character, intent(in) :: part

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  select case(part)
    case('d', 'D')
      call mesh%RunCellKernel(DiagMatrixVector_Kernel)
    case('l', 'L')
      call mesh%RunCellKernel(LowerMatrixVector_Kernel)
    case('u', 'U')
      call mesh%RunCellKernel(UpperMatrixVector_Kernel)
  end select

contains
subroutine DiagMatrixVector_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col

    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      col = mat%ColIndices(rowAddr)
      if (row == col) then
        y(:,row) = matmul(mat%ColCoeffs(:,:,rowAddr), x(:,col))
      end if
    end do

  end subroutine DiagMatrixVector_Kernel
  subroutine LowerMatrixVector_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col

    y(:,row) = 0.0_dp

    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      col = mat%ColIndices(rowAddr)
      if (row < col) then
        y(:,row) = y(:,row) + matmul(mat%ColCoeffs(:,:,rowAddr), x(:,col))
      end if
    end do

  end subroutine LowerMatrixVector_Kernel
  subroutine UpperMatrixVector_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col

    y(:,row) = 0.0_dp

    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      col = mat%ColIndices(rowAddr)
      if (row > col) then
        y(:,row) = y(:,row) + matmul(mat%ColCoeffs(:,:,rowAddr), x(:,col))
      end if
    end do

  end subroutine UpperMatrixVector_Kernel
end subroutine PartialMatrixVector

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve equation 𝓓𝒙 = 𝒚, where 𝓓 is the (block-)diagonal of 𝓐.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SolveDiag(mesh, mat, yArr, xArr)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tArray), intent(inout) :: xArr, yArr

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(SolveDiag_Kernel)

contains
  subroutine SolveDiag_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col

    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      col = mat%ColIndices(rowAddr)
      if (row == col) then
        y(:,row) = DenseSolve(mat%ColCoeffs(:,:,rowAddr), x(:,col))
      end if
    end do

  end subroutine SolveDiag_Kernel
end subroutine SolveDiag

subroutine SolveTrianular(mesh, part, mat, yArr, xArr)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  character, intent(in) :: part

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mkl_dcsrtrsv(part, 'N', 'N', mesh%NumCells, &
      & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, x, y)

end subroutine SolveTrianular

end module StormRuler_Matrix
