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
module StormRuler_FDM_Convection

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8, gCylCoords
use StormRuler_Helpers, only: Flip

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray

use StormRuler_FDM_Operators, only: gTruncErrorOrder
use StormRuler_FDM_Base_Flux

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface FDM_Convection_Central
  module procedure FDM_Convection_Central
end interface FDM_Convection_Central

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate convection: 𝒗 ← 𝒗 - 𝜆∇⋅(𝒂⃗𝒖).
!! Shape of 𝒖, 𝒗 is [1, NumVars]×[1, NumAllCells],
!! shape of 𝒂 is [1, NumDims]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Convection_Central(mesh, vArr, lambda, uArr, aArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: uArr, aArr
  class(tArray), intent(inout) :: vArr
  real(dp), intent(in) :: lambda

  real(dp), pointer :: u(:,:), v(:,:), aVec(:,:)

  call uArr%Get(u); call vArr%Get(v); call aArr%Get(aVec)

  call mesh%RunCellKernel(FDM_Convection_Central_Kernel)

contains
  subroutine FDM_Convection_Central_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    integer(ip) :: dim
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell

    ! ----------------------
    ! For each direction do:
    ! ----------------------
    do dim = 1, mesh%NumDims
      iCellFace = 2*(dim - 1) + 1

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
              & lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (gTruncErrorOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (gTruncErrorOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (gTruncErrorOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate convection increment.
      ! ----------------------
      associate(dlInv => lambda/mesh%dl(iCellFace))
        select case(gTruncErrorOrder)
          case(1:2)
            v(:,iCell) = v(:,iCell) + &
              & dlInv*( FR_C2(aVec(dim,lCell), u(:,lCell), &
              &               aVec(dim,iCell), u(:,iCell)) &
              &       - FR_C2(aVec(dim,iCell), u(:,iCell), &
              &               aVec(dim,rCell), u(:,rCell)) )
          ! ----------------------
          case(3:4)
            v(:,iCell) = v(:,iCell) + &
              & dlInv*( FR_C4(aVec(dim, rCell), u(:, rCell), &
              &               aVec(dim, iCell), u(:, iCell), &
              &               aVec(dim, lCell), u(:, lCell), &
              &               aVec(dim,llCell), u(:,llCell)) &
              &       - FR_C4(aVec(dim,rrCell), u(:,rrCell), &
              &               aVec(dim, rCell), u(:, rCell), &
              &               aVec(dim, iCell), u(:, iCell), &
              &               aVec(dim, lCell), u(:, lCell)) )
          ! ----------------------
          case(5:6)
            v(:,iCell) = v(:,iCell) + &
              & dlInv*( FR_C6(aVec(dim, rrCell), u(:, rrCell), &
              &               aVec(dim,  rCell), u(:,  rCell), &
              &               aVec(dim,  iCell), u(:,  iCell), &
              &               aVec(dim,  lCell), u(:,  lCell), &
              &               aVec(dim, llCell), u(:, llCell), &
              &               aVec(dim,lllCell), u(:,lllCell)) &
              &       - FR_C6(aVec(dim,rrrCell), u(:,rrrCell), &
              &               aVec(dim, rrCell), u(:, rrCell), &
              &               aVec(dim,  rCell), u(:,  rCell), &
              &               aVec(dim,  iCell), u(:,  iCell), &
              &               aVec(dim,  lCell), u(:,  lCell), &
              &               aVec(dim, llCell), u(:, llCell)) )
          ! ----------------------
          case(7:8)
            v(:,iCell) = v(:,iCell) + &
              & dlInv*( FR_C8(aVec(dim, rrrCell), u(:, rrrCell), &
              &               aVec(dim,  rrCell), u(:,  rrCell), &
              &               aVec(dim,   rCell), u(:,   rCell), &
              &               aVec(dim,   iCell), u(:,   iCell), &
              &               aVec(dim,   lCell), u(:,   lCell), &
              &               aVec(dim,  llCell), u(:,  llCell), &
              &               aVec(dim, lllCell), u(:, lllCell), &
              &               aVec(dim,llllCell), u(:,llllCell)) &
              &       - FR_C8(aVec(dim,rrrrCell), u(:,rrrrCell), &
              &               aVec(dim, rrrCell), u(:, rrrCell), &
              &               aVec(dim,  rrCell), u(:,  rrCell), &
              &               aVec(dim,   rCell), u(:,   rCell), &
              &               aVec(dim,   iCell), u(:,   iCell), &
              &               aVec(dim,   lCell), u(:,   lCell), &
              &               aVec(dim,  llCell), u(:,  llCell), &
              &               aVec(dim, lllCell), u(:, lllCell)) )
        end select
      end associate
    end do

  end subroutine FDM_Convection_Central_Kernel
end subroutine FDM_Convection_Central

end module StormRuler_FDM_Convection
