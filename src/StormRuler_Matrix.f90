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
use StormRuler_Parameters, only: gUseMKL

use StormRuler_Helpers, only: IndexOf, BubbleSort, &
  & InverseCompressMapping, DenseSolve

use StormRuler_Mesh, only: tMesh, tKernelFunc
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: Fill
use StormRuler_BLAS, only: tMatVecFunc

#$if HAS_MKL
use StormRuler_MKL, only: mkl_dcsrgemv, mkl_dbsrgemv, &
  & mkl_dcsrtrsv, mkl_dbsrtrsv
#$end if

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Block compressed sparse row (CSR) block matrix.
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

contains

  ! ----------------------
  ! Get block size of a matrix.
  ! ----------------------
  procedure :: BlockSize => GetMatrixBlockSize

  ! ----------------------
  ! Check if matrix is a block matrix.
  ! ----------------------
  procedure :: IsBlock => IsBlockMatrix

end type tMatrix

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Parallelization information for the triangular solvers.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tParallelTriangularContext

  ! ----------------------
  ! Number of levels in the DAG.
  ! ----------------------
  integer(ip) :: NumLevels

  integer(ip), allocatable :: LevelAddrs(:)

  integer(ip), allocatable :: LevelRowIndices(:)

end type tParallelTriangularContext

interface SolveTriangular
  module procedure SolveTriangular
  module procedure ParallelSolveTriangular
end interface SolveTriangular

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
!! Get block size of a matrix.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
pure integer(ip) function GetMatrixBlockSize(mat) result(blockSize)
  class(tMatrix), intent(in) :: mat

  blockSize = size(mat%ColCoeffs(:,1,1))

end function GetMatrixBlockSize

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Check if matrix is block.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
pure logical function IsBlockMatrix(mat)
  class(tMatrix), intent(in) :: mat

  IsBlockMatrix = mat%BlockSize() /= 1

end function IsBlockMatrix

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

#$if HAS_MKL
  if (gUseMKL) then
    if (mat%IsBlock()) then
      call mkl_dbsrgemv('N', mesh%NumCells, mat%BlockSize(), &
        & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, x, y)
    else
      call mkl_dcsrgemv('N', mesh%NumCells, &
        & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, x, y)
    end if
    return
  end if
#$end if

  !! TODO: optimization in the non-block case.
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

  !! TODO: Does the MKL implementation for these operation exist?

  !! TODO: optimization in the non-block case.
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
      if (col == row) then
        y(:,row) = matmul(mat%ColCoeffs(:,:,rowAddr), x(:,col))
        return
      end if
    end do

  end subroutine DiagMatrixVector_Kernel
  subroutine LowerMatrixVector_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col

    y(:,row) = 0.0_dp

    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      col = mat%ColIndices(rowAddr)
      if (col < row) then
        y(:,row) = y(:,row) + matmul(mat%ColCoeffs(:,:,rowAddr), x(:,col))
      else
        return
      end if
    end do

  end subroutine LowerMatrixVector_Kernel
  subroutine UpperMatrixVector_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col

    y(:,row) = 0.0_dp

    do rowAddr = mat%RowAddrs(row + 1) - 1, mat%RowAddrs(row), -1
      col = mat%ColIndices(rowAddr)
      if (col > row) then
        y(:,row) = y(:,row) + matmul(mat%ColCoeffs(:,:,rowAddr), x(:,col))
      else
        return
      end if
    end do

  end subroutine UpperMatrixVector_Kernel
end subroutine PartialMatrixVector

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! (Block-)diagonal solver helper.
!! ----------------------------------------------------------------- !!
subroutine SolveDiagHelper(rowAddrs, colIndices, colCoeffs, y, b, row)
  integer(ip), intent(in) :: rowAddrs(*), colIndices(*)
  real(dp), intent(in) :: colCoeffs(*), b(*)
  real(dp), intent(inout) :: y(*)
  integer(ip), intent(in) :: row

  integer(ip) :: rowAddr, col

  do rowAddr = rowAddrs(row), rowAddrs(row + 1) - 1
    col = colIndices(rowAddr)
    if (row == col) then
      y(row) = b(col)/colCoeffs(rowAddr)
      return
    end if
  end do

end subroutine SolveDiagHelper
subroutine SolveBlockDiagHelper(size, rowAddrs, colIndices, colCoeffs, y, b, row)
  integer(ip), intent(in) :: size, rowAddrs(*), colIndices(*)
  real(dp), intent(in) :: colCoeffs(size,size,*), b(size,*)
  real(dp), intent(inout) :: y(size,*)
  integer(ip), intent(in) :: row

  integer(ip) :: rowAddr, col

  do rowAddr = rowAddrs(row), rowAddrs(row + 1) - 1
    col = colIndices(rowAddr)
    if (row == col) then
      y(:,row) = DenseSolve(colCoeffs(:,:,rowAddr), b(:,col))
      return
    end if
  end do

end subroutine SolveBlockDiagHelper

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve equation 𝓓𝒚 = 𝒃, where 𝓓 is the (block-)diagonal of 𝓐.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SolveDiag(mesh, mat, yArr, bArr)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tArray), intent(inout) :: bArr, yArr

  integer(ip) :: size
  real(dp), pointer :: b(:,:), y(:,:)

  size = mat%BlockSize()
  call bArr%Get(b); call yArr%Get(y)

  !! TODO: Does the MKL implementation for this operation exist?
  if (size == 1) then
    call mesh%RunCellKernel(SolveDiag_Kernel)
  else
    call mesh%RunCellKernel(SolveBlockDiag_Kernel)
  end if

contains
  subroutine SolveDiag_Kernel(row)
    integer(ip), intent(in) :: row

    call SolveDiagHelper( &
      & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)

  end subroutine SolveDiag_Kernel
  subroutine SolveBlockDiag_Kernel(row)
    integer(ip), intent(in) :: row

    call SolveBlockDiagHelper(size, &
      & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)

  end subroutine SolveBlockDiag_Kernel
end subroutine SolveDiag

!! ----------------------------------------------------------------- !!
!! (Block-)lower triangular solver helper.
!! ----------------------------------------------------------------- !!
subroutine SolveLowerTriangHelper(rowAddrs, colIndices, colCoeffs, y, b, row)
  integer(ip), intent(in) :: rowAddrs(*), colIndices(*)
  real(dp), intent(in) :: colCoeffs(*), b(*)
  real(dp), intent(inout) :: y(*)
  integer(ip), intent(in) :: row

  integer(ip) :: rowAddr, col

  y(row) = 0.0

  do rowAddr = rowAddrs(row), rowAddrs(row + 1) - 1
    col = colIndices(rowAddr)
    if (col < row) then
      y(row) = y(row) + colCoeffs(rowAddr)*y(col)
    else if (col == row) then
      y(row) = (b(row) - y(row))/colCoeffs(rowAddr)
      return
    end if
  end do

end subroutine SolveLowerTriangHelper
subroutine SolveBlockLowerTriangHelper(size, rowAddrs, colIndices, colCoeffs, y, b, row)
  integer(ip), intent(in) :: size, rowAddrs(*), colIndices(*)
  real(dp), intent(in) :: colCoeffs(size,size,*), b(size,*)
  real(dp), intent(inout) :: y(size,*)
  integer(ip), intent(in) :: row

  integer(ip) :: rowAddr, col

  y(:,row) = 0.0

  do rowAddr = rowAddrs(row), rowAddrs(row + 1) - 1
    col = colIndices(rowAddr)
    if (col < row) then
      y(:,row) = y(:,row) + matmul(colCoeffs(:,:,rowAddr), y(:,col))
    else if (col == row) then
      y(:,row) = DenseSolve(colCoeffs(:,:,rowAddr), b(:,row) - y(:,row))
      return
    end if
  end do

end subroutine SolveBlockLowerTriangHelper

!! ----------------------------------------------------------------- !!
!! (Block-)upper triangular solver helper.
!! ----------------------------------------------------------------- !!
subroutine SolveUpperTriangHelper(rowAddrs, colIndices, colCoeffs, y, b, row)
  integer(ip), intent(in) :: rowAddrs(*), colIndices(*)
  real(dp), intent(in) :: colCoeffs(*), b(*)
  real(dp), intent(inout) :: y(*)
  integer(ip), intent(in) :: row

  integer(ip) :: rowAddr, col

  y(row) = 0.0

  do rowAddr = rowAddrs(row + 1) - 1, rowAddrs(row), -1
    col = colIndices(rowAddr)
    if (col > row) then
      y(row) = y(row) + colCoeffs(rowAddr)*y(col)
    else if (col == row) then
      y(row) = (b(row) - y(row))/colCoeffs(rowAddr)
      return
    end if
  end do

end subroutine SolveUpperTriangHelper
subroutine SolveBlockUpperTriangHelper(size, rowAddrs, colIndices, colCoeffs, y, b, row)
  integer(ip), intent(in) :: size, rowAddrs(*), colIndices(*)
  real(dp), intent(in) :: colCoeffs(size,size,*), b(size,*)
  real(dp), intent(inout) :: y(size,*)
  integer(ip), intent(in) :: row

  integer(ip) :: rowAddr, col

  y(:,row) = 0.0

  do rowAddr = rowAddrs(row + 1) - 1, rowAddrs(row), -1
    col = colIndices(rowAddr)
    if (col > row) then
      y(:,row) = y(:,row) + matmul(colCoeffs(:,:,rowAddr), y(:,col))
    else if (col == row) then
      y(:,row) = DenseSolve(colCoeffs(:,:,rowAddr), b(:,row) - y(:,row))
      return
    end if
  end do

end subroutine SolveBlockUpperTriangHelper

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve equation (𝓛 + 𝓓)𝒚 = 𝒃, or (𝓓 + 𝓤)𝒚 = 𝒃,
!! where 𝓓 is the (block-)diagonal of 𝓐, 𝓛 and 𝓤 are lower and upper 
!! strict (block-)triangular parts of 𝓐, 𝓐 = 𝓛 + 𝓓 + 𝓤.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SolveTriangular(mesh, part, mat, yArr, bArr)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tArray), intent(in) :: bArr
  class(tArray), intent(inout) :: yArr
  character, intent(in) :: part

  integer(ip) :: size, row
  real(dp), pointer :: b(:,:), y(:,:)

  size = mat%BlockSize()
  call bArr%Get(b); call yArr%Get(y)

#$if HAS_MKL
  if (gUseMKL) then
    select case(part)
      case('l', 'L')
        if (size == 1) then
          call mkl_dcsrtrsv('L', 'N', 'N', mesh%NumCells, &
            & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, b, y)
        else
          call mkl_dbsrtrsv('L', 'N', 'N', mesh%NumCells, size, &
            & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, b, y)
        end if
      case('u', 'U')
        if (size == 1) then
          call mkl_dcsrtrsv('U', 'N', 'N', mesh%NumCells, &
            & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, b, y)
        else
          call mkl_dbsrtrsv('U', 'N', 'N', mesh%NumCells, size, &
            & mat%ColCoeffs, mat%RowAddrs, mat%ColIndices, b, y)
        end if
    end select
    return
  end if
#$end if

  select case(part)
    case('l', 'L')
      if (size == 1) then
        do row = 1, mesh%NumCells
          call SolveLowerTriangular_SeqKernel(row)
        end do
      else
        do row = 1, mesh%NumCells
          call SolveBlockLowerTriangular_SeqKernel(row)
        end do
      end if
    case('u', 'U')
      if (size == 1) then
        do row = mesh%NumCells, 1, -1
          call SolveUpperTriangular_SeqKernel(row)
        end do
      else
        do row = mesh%NumCells, 1, -1
          call SolveBlockUpperTriangular_SeqKernel(row)
        end do
      end if
  end select

contains
  subroutine SolveLowerTriangular_SeqKernel(row)
    integer(ip), intent(in) :: row

    call SolveLowerTriangHelper( &
      & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)

  end subroutine SolveLowerTriangular_SeqKernel
  subroutine SolveBlockLowerTriangular_SeqKernel(row)
    integer(ip), intent(in) :: row

    call SolveBlockLowerTriangHelper(size, &
      & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)

  end subroutine SolveBlockLowerTriangular_SeqKernel
  subroutine SolveUpperTriangular_SeqKernel(row)
    integer(ip), intent(in) :: row

    call SolveUpperTriangHelper( &
      & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)

  end subroutine SolveUpperTriangular_SeqKernel
  subroutine SolveBlockUpperTriangular_SeqKernel(row)
    integer(ip), intent(in) :: row

    call SolveBlockUpperTriangHelper(size, &
      & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)

  end subroutine SolveBlockUpperTriangular_SeqKernel
end subroutine SolveTriangular

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize the parallel triangular context.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitParallelContext(mesh, part, mat, ctx)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tParallelTriangularContext), intent(inout) :: ctx
  character, intent(in) :: part

  integer(ip) :: row, width
  integer(ip), allocatable :: rowLevels(:)

  ctx%NumLevels = 0
  allocate(rowLevels(mesh%NumCells))
  rowLevels(:) = 0
  select case(part)
    case('l', 'L')
      do row = 1, mesh%NumCells
        call InitParallelLowerTriangularContext_SeqKernel(row)
      end do
    case('u', 'U')
      do row = mesh%NumCells, 1, -1
        call InitParallelUpperTriangularSolver_SeqKernel(row)
      end do
  end select

  call InverseCompressMapping(ctx%LevelAddrs, &
    & ctx%LevelRowIndices, ctx%NumLevels, rowLevels, width)
  
  print *, 'num levels: ', ctx%NumLevels, 'width: ', width

contains
  subroutine InitParallelLowerTriangularContext_SeqKernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col, level

    rowLevels(row) = 0
    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      col = mat%ColIndices(rowAddr)
      if (col < row) then
        rowLevels(row) = max(rowLevels(row), rowLevels(col))
      else if (col == row) then
        rowLevels(row) = rowLevels(row) + 1
        ctx%NumLevels = max(ctx%NumLevels, rowLevels(row))
        return
      end if
    end do


  end subroutine InitParallelLowerTriangularContext_SeqKernel
  subroutine InitParallelUpperTriangularSolver_SeqKernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr, col, level

    rowLevels(row) = 0
    do rowAddr = mat%RowAddrs(row + 1) - 1, mat%RowAddrs(row), -1
      col = mat%ColIndices(rowAddr)
      if (col > row) then
        rowLevels(row) = max(rowLevels(row), rowLevels(col))
      else if (col == row) then
        rowLevels(row) = rowLevels(row) + 1
        ctx%NumLevels = max(ctx%NumLevels, rowLevels(row))
        return
      end if
    end do

  end subroutine InitParallelUpperTriangularSolver_SeqKernel
end subroutine InitParallelContext

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve equation (𝓛 + 𝓓)𝒚 = 𝒃, or (𝓓 + 𝓤)𝒚 = 𝒃,
!! where 𝓓 is the (block-)diagonal of 𝓐, 𝓛 and 𝓤 are lower and upper 
!! strict (block-)triangular parts of 𝓐, 𝓐 = 𝓛 + 𝓓 + 𝓤.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ParallelSolveTriangular(mesh, part, mat, ctx, yArr, bArr)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tArray), intent(in) :: bArr
  class(tArray), intent(inout) :: yArr
  class(tParallelTriangularContext), intent(inout) :: ctx
  character, intent(in) :: part

  integer(ip) :: size
  real(dp), pointer :: b(:,:), y(:,:)

  !call SolveTriangular(mesh, part, mat, yArr, bArr)
  !return

  size = mat%BlockSize()
  call bArr%Get(b); call yArr%Get(y)

  select case(part)
    case('l', 'L')
      if (size == 1) then
        call ParallelSolveLowerTriangular()
      else
        call ParallelSolveBlockLowerTriangular()
      end if
    case('u', 'U')
      if (size == 1) then
        call ParallelSolveUpperTriangular()
      else
        call ParallelSolveBlockUpperTriangular()
      end if
  end select

contains
  subroutine ParallelSolveLowerTriangular()
    integer(ip) :: row, level, levelAddr

    !$omp parallel
    do level = 1, ctx%NumLevels
      !$omp do private(row,levelAddr) schedule(static)
      do levelAddr = ctx%LevelAddrs(level), ctx%LevelAddrs(level + 1) - 1

        row = ctx%LevelRowIndices(levelAddr)
        call SolveLowerTriangHelper( &
          & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)

      end do
      !$omp end do
    end do
    !$omp end parallel

  end subroutine ParallelSolveLowerTriangular
  subroutine ParallelSolveBlockLowerTriangular()
    integer(ip) :: row, level, levelAddr

    !$omp parallel
    do level = 1, ctx%NumLevels
      !$omp do private(row,levelAddr) schedule(static)
      do levelAddr = ctx%LevelAddrs(level), ctx%LevelAddrs(level + 1) - 1
        
        row = ctx%LevelRowIndices(levelAddr)
        call SolveBlockLowerTriangHelper(size, &
          & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)
  
      end do
      !$omp end do
    end do
    !$omp end parallel

  end subroutine ParallelSolveBlockLowerTriangular
  subroutine ParallelSolveUpperTriangular()
    integer(ip) :: row, level, levelAddr

    !$omp parallel
    do level = 1, ctx%NumLevels
      !$omp do private(row,levelAddr) schedule(static)
      do levelAddr = ctx%LevelAddrs(level), ctx%LevelAddrs(level + 1) - 1
  
        row = ctx%LevelRowIndices(levelAddr)
        call SolveUpperTriangHelper( &
          & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)
  
      end do
      !$omp end do
    end do
    !$omp end parallel

  end subroutine ParallelSolveUpperTriangular
  subroutine ParallelSolveBlockUpperTriangular()

    integer(ip) :: row, level, levelAddr

    !$omp parallel
    do level = 1, ctx%NumLevels
      !$omp do private(row,levelAddr) schedule(static)
      do levelAddr = ctx%LevelAddrs(level), ctx%LevelAddrs(level + 1) - 1

        row = ctx%LevelRowIndices(levelAddr)
        call SolveBlockUpperTriangHelper(size, &
          & mat%RowAddrs, mat%ColIndices, mat%ColCoeffs, y, b, row)
  
      end do
      !$omp end do
    end do
    !$omp end parallel

  end subroutine ParallelSolveBlockUpperTriangular
end subroutine ParallelSolveTriangular

end module StormRuler_Matrix
