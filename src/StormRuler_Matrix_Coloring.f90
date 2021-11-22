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
module StormRuler_Matrix_Coloring

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
  integer(ip) :: NumLabels

  ! ----------------------
  ! Shape is [1, Cells].
  ! ----------------------
  integer(ip), allocatable :: ColLbls(:)

  ! ----------------------
  ! Label addresses.
  ! Shape is [1, NumLabels + 1].
  ! ----------------------
  integer(ip), allocatable :: LabelAddrs(:)
  ! ----------------------
  ! Label column indices.
  ! Shape is [1, Cells].
  ! ----------------------
  integer(ip), allocatable :: LabelColIndices(:)

end type tMatrixLabeling

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Generate column labeling using the banded algorithm.
!!
!! Estimated number of label is ≈ 𝘣𝘸(𝓐ᵖ) = 𝑝⋅(𝘣𝘸(𝓐) - 1) + 1,
!! where 𝘣𝘸(𝓐) is the bandwidth of 𝓐.
!! 
!! Banded algorithm is simple and fast, but suboptimal in 
!! the most practical cases (for 𝑝 = 1 it generates ≈2 times
!! more labels than more sophisticated labeling algorithms).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LabelColumns_Banded(mesh, mat, labeling)
  class(tMesh), intent(in) :: mesh
  class(tMatrix), intent(in) :: mat
  class(tMatrixLabeling), intent(inout) :: labeling

  integer(ip) :: row, rowAddr, col, lbl

  ! ----------------------
  ! Compute number of colors as matrix bandwidth.
  ! ----------------------
  labeling%NumLabels = 0
  do row = 1, mesh%NumCells
    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1
      associate(col => mat%ColIndices(rowAddr))

        labeling%NumLabels = max(labeling%NumLabels, abs(row - col))

      end associate
    end do
  end do
  labeling%NumLabels = 2*labeling%NumLabels + 1

  ! ----------------------
  ! Fill the column to label table.
  ! ----------------------
  allocate(labeling%ColLbls(mesh%NumCells))
  do col = 1, mesh%NumCells

    lbl = mod(col - 1, labeling%NumLabels) + 1 
    labeling%ColLbls(col) = lbl
  
  end do

  ! ----------------------
  ! Fill the label to column table.
  ! ----------------------
  allocate(labeling%LabelAddrs(labeling%NumLabels + 1))
  allocate(labeling%LabelColIndices(mesh%NumCells))
  labeling%LabelAddrs(1) = 1
  do lbl = 1, labeling%NumLabels
    associate(lblAddr => labeling%LabelAddrs(lbl + 1))
      
      lblAddr = labeling%LabelAddrs(lbl)
      do col = lbl, mesh%NumCells, labeling%NumLabels
        labeling%LabelColIndices(lblAddr) = col 
        lblAddr = lblAddr + 1
      end do

    end associate
  end do

end subroutine LabelColumns_Banded

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

subroutine ExtractMatrix(mesh, mat, labeling, MatVec, mold)
  class(tMesh), intent(inout) :: mesh
  class(tMatrix), intent(inout) :: mat
  class(tMatrixLabeling), intent(in) :: labeling
  class(tArray), intent(in) :: mold
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: row, rowAddr, col, lbl, lblAddr

  real(dp), pointer :: t(:,:), Q(:,:,:,:)
  type(tArray) :: QArr, tArr, sArr

  call AllocArray(tArr, mold=mold)
  call AllocArray(QArr, shape=[mold%mShape, labeling%NumLabels])
  call QArr%Get(q)

  ! ----------------------
  ! Reconstruct the compressed columns.
  ! ----------------------
  do lbl = 1, labeling%NumLabels

    call Fill(mesh, tArr, 0.0_dp)

    do lblAddr = labeling%LabelAddrs(lbl), labeling%LabelAddrs(lbl + 1) - 1
      col = labeling%LabelColIndices(lblAddr)
      call tArr%Get(t); t(:,col) = 1.0_dp
    end do
  
    sArr = QArr%Slice(lbl); call MatVec(mesh, sArr, tArr)
  
  end do

  ! ----------------------
  ! Decompress the columns to the CSR matrix.
  ! ----------------------
  do row = 1, mesh%NumCells
    do rowAddr = mat%RowAddrs(row), mat%RowAddrs(row + 1) - 1

      lbl = labeling%ColLbls(mat%ColIndices(rowAddr))
      mat%ColCoeffs(:,:,rowAddr) = Q(:,:,row,lbl)

    end do
  end do

end subroutine ExtractMatrix

end module StormRuler_Matrix_Coloring
