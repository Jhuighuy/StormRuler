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
use StormRuler_Matrix, only: tColumnMatrix

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
  integer(ip), allocatable :: ColorPtrs(:)
  ! ----------------------
  ! Shape is [1, Cells].
  ! ----------------------
  integer(ip), allocatable :: ColorColumnIndices(:)

end type tColumnColoring

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

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

  integer(ip) :: column, color
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
  do column = 1, mesh%NumCells

    color = mod(column - 1, coloring%NumColors) + 1 
    coloring%ColumnColor(column) = color

  end do

  ! ----------------------
  ! Fill the color to column table.
  ! ----------------------
  allocate(coloring%ColorPtrs(coloring%NumColors + 1))
  allocate(coloring%ColorColumnIndices(mesh%NumCells))

  coloring%ColorPtrs(1) = 1
  do color = 1, coloring%NumColors
    
    associate(colorPtr => coloring%ColorPtrs(color + 1))
      colorPtr = coloring%ColorPtrs(color)
      do column = color, mesh%NumCells, bandwidth
        coloring%ColorColumnIndices(colorPtr) = column 
        colorPtr = colorPtr + 1
      end do
    end associate

  end do

end subroutine ColorColumns_Banded

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

subroutine ReconstructMatrix(matrix, coloring, mesh, MatVec, mold)
  class(tColumnMatrix), intent(inout) :: matrix
  class(tColumnColoring), intent(in) :: coloring
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: mold
  procedure(tMatVecFunc) :: MatVec

  integer(ip) :: row, column, columnPtr, color, colorPtr
  real(dp), pointer :: tDat(:,:)
  type(tArray) :: Q, t, s

  call AllocArray(t, mold=mold)
  call AllocArray(Q, shape=[mold%mShape, coloring%NumColors])

  ! ----------------------
  ! Reconstruct the compressed columns.
  ! ----------------------
  do color = 1, coloring%NumColors

    call Fill(mesh, t, 0.0_dp)
    do colorPtr = coloring%ColorPtrs(color), coloring%ColorPtrs(color + 1) - 1
      column = coloring%ColorColumnIndices(colorPtr)
      call t%Get(tDat); tDat(:,column) = 1.0_dp
    end do
  
    s = Q%Slice(color); call MatVec(mesh, s, t)
  
  end do

  ! ----------------------
  ! Decompress the columns.
  ! ----------------------
  do column = 1, mesh%NumCells

    color = coloring%ColumnColor(column)
    t = Q%Slice(color); call t%Get(tDat)
    do columnPtr = matrix%ColumnPtrs(column), matrix%ColumnPtrs(column + 1) - 1
      matrix%RowCoeffs(:,:,columnPtr) = tDat(1,matrix%RowIndices(columnPtr))
    end do

  end do

end subroutine ReconstructMatrix

end module StormRuler_Matrix_Coloring
