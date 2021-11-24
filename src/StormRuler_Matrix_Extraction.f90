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
module StormRuler_Matrix_Extraction

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: tMatVecFunc, Fill
use StormRuler_Matrix!, only: tColumnMatrix

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Marix column labeling.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tMatrixLabeling

  ! ----------------------
  ! Number of the colors.
  ! ----------------------
  integer(ip) :: NumLabels = 0

  ! ----------------------
  ! Shape is [1, Cells].
  ! ----------------------
  integer(ip), allocatable :: ColLabels(:)

end type tMatrixLabeling

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Generate column labeling using the banded algorithm.
!!
!! Estimated number of label is equal to the bandwidth of 𝓐.
!! 
!! Banded algorithm is simple and fast, but generated the 
!! the suboptimal labeling in the most practical cases.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LabelColumns_Banded(mesh, mat, labeling)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tMatrixLabeling), intent(inout) :: labeling

  integer(ip) :: row, rowAddr, col

  ! ----------------------
  ! Compute number of colors as matrix bandwidth.
  ! ----------------------
  labeling%NumLabels = 0
  do row = 1, mesh%NumCells
    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1

      col = mat%ColIndices(rowAddr)
      labeling%NumLabels = max(labeling%NumLabels, abs(row - col))

    end do
  end do
  labeling%NumLabels = 2*labeling%NumLabels + 1

  ! ----------------------
  ! Label columns.
  ! ----------------------
  allocate(labeling%ColLabels(mesh%NumCells))
  do col = 1, mesh%NumCells

    ! ----------------------
    ! Columns with distance more than bandwidth 
    ! are guaranteed not to intersect.
    ! ----------------------
    labeling%ColLabels(col) = mod(col - 1, labeling%NumLabels) + 1

  end do

end subroutine LabelColumns_Banded

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Generate column labeling using the pattern-based greedy algorithm.
!!
!! Estimated number of labels is approximately equal to
!! the maximum amount of the non-zero entries per row of 𝓐.
!!
!! Labeling construction may be slow.
!!
!! References:
!! [1] Curtis, A. R. and Powell M. J. D., and Reid J. K. 
!!     “On the Estimation of Sparse Jacobian Matrices.“ 
!!     IMA Journal of Applied Mathematics 13 (1974): 117–119.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LabelColumns_Patterned(mesh, mat, labeling)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tMatrixLabeling), intent(inout) :: labeling

  type :: tMatrixColumn
    integer(ip), allocatable :: RowIndices(:)
    type(tMatrixColumn), pointer :: NextCol => null()
  end type tMatrixColumn

  integer(ip) :: row, rowAddr, col, colCol
  type(tMatrixColumn), allocatable, target :: matCols(:)
  type(tMatrixColumn), pointer :: matCol

  ! ----------------------
  ! Generate matrix columns.
  ! ----------------------
  allocate(matCols(mesh%NumCells))
  do row = 1, mesh%NumCells
    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      col = mat%ColIndices(rowAddr)
      !! TODO: refactor with `SortedInsert`.
      if (allocated(matCols(col)%RowIndices)) then
        matCols(col)%RowIndices = [matCols(col)%RowIndices, row]
      else
        matCols(col)%RowIndices = [row]
      end if
    end do
  end do
  do col = 1, mesh%NumCells
    call BubbleSort(matCols(col)%RowIndices)
  end do

  ! ----------------------
  ! Label columns with greedy algorithm.
  ! ----------------------
  allocate(labeling%ColLabels(mesh%NumCells))
  labeling%NumLabels = 0; labeling%ColLabels(:) = 0
  do col = 1, mesh%NumCells
    if (labeling%ColLabels(col) /= 0) cycle

    ! ----------------------
    ! Assign the label.
    ! ----------------------
    labeling%ColLabels(col) = labeling%NumLabels + 1
    labeling%NumLabels = labeling%ColLabels(col)

    ! ----------------------
    ! Try to assign the same label to the following column.
    ! ----------------------
    do colCol = col + 1, mesh%NumCells
      if (labeling%ColLabels(colCol) /= 0) cycle
      labeling%ColLabels(colCol) = labeling%ColLabels(col)
      
      ! ----------------------
      ! Check for intersection with the previously assigned columns.
      ! ----------------------
      matCol => matCols(col)
      do while (associated(matCol))
        if (ArraysIntersect(matCol%RowIndices, matCols(colCol)%RowIndices)) then
          labeling%ColLabels(colCol) = 0; exit
        end if
        matCol => matCol%NextCol
      end do

      ! ----------------------
      ! If no intersection detected, link the current column with the first one.
      ! ----------------------
      if (.not.associated(matCol)) then
        matCols(colCol)%NextCol => matCols(col)%NextCol
        matCols(col)%NextCol => matCols(colCol)
      end if

    end do
  end do

end subroutine LabelColumns_Patterned

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Extract matrix from the matrix-vector product function
!! with the specified column labeling.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ExtractMatrix(mesh, mat, labeling, MatVec, mold)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(inout) :: mat
  class(tMatrixLabeling), intent(in) :: labeling
  class(tArray), intent(in) :: mold
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: label

  real(dp), pointer :: t(:,:), Q(:,:,:,:)
  type(tArray) :: QArr, tArr, sArr

  call AllocArray(tArr, mold=mold)
  call AllocArray(QArr, shape=[mold%mShape, labeling%NumLabels])
  call tArr%Get(t); call QArr%Get(q)

  ! ----------------------
  ! Prepare the probes and extract the merged columns.
  ! ----------------------
  do label = 1, labeling%NumLabels
    call mesh%RunCellKernel(PrepareProbe_Kernel)
    sArr = QArr%Slice(label); call MatVec(mesh, sArr, tArr)
  end do

  ! ----------------------
  ! Compress the rows to the CSR matrix.
  ! ----------------------
  call mesh%RunCellKernel(CompressRows_Kernel)

contains
  subroutine PrepareProbe_Kernel(col)
    integer(ip), intent(in) :: col

    t(:,col) = merge(1.0_dp, 0.0_dp, labeling%ColLabels(col) == label)

  end subroutine PrepareProbe_Kernel
  subroutine CompressRows_Kernel(row)
    integer(ip), intent(in) :: row

    integer(ip) :: rowAddr

    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      associate(label => labeling%ColLabels(mat%ColIndices(rowAddr)))
        mat%ColCoeffs(:,:,rowAddr) = Q(:,:,row,label)
      end associate
    end do

  end subroutine CompressRows_Kernel
end subroutine ExtractMatrix

end module StormRuler_Matrix_Extraction
