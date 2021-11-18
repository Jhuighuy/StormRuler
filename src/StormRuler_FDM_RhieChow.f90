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
module StormRuler_FDM_RhieChow

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8, gCylCoords
use StormRuler_Helpers, only: Flip

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray

use StormRuler_FDM_Operators, only: gTruncErrorOrder
use StormRuler_FDM_Base, only: WFD4_C2

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface FDM_RhieChow_Correction
  module procedure FDM_RhieChow_Correction
end interface FDM_RhieChow_Correction

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate Rhie-Chow type correction: 𝒗 ← 𝒗 - 𝜆𝜏𝓡𝓒(𝒑,𝛒),
!! where 𝓡𝓒(𝒑,𝛒) = ∂(1/𝛒⋅∂³𝒑/∂𝑥³)/∂𝑥 + … + ∂(1/𝛒⋅∂³𝒑/∂𝑧)/∂𝑧.
!!
!! Rhie-Chow correction is used to eleminate the checkerboard 
!! pressure phenomenon that may lead to the pressure-velocity 
!! decoupling in the incompressible simulations.
!!
!! shape of 𝒗,𝒑,𝛒 is [1,NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_RhieChow_Correction(mesh, vArr, lambda, tau, pArr, rhoAny)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: pArr
  class(tArray), intent(in) :: rhoAny
  class(tArray), intent(inout) :: vArr
  real(dp), intent(in) :: lambda, tau

  real(dp), pointer :: v(:,:), p(:), rho(:)

  call vArr%Get(v); call pArr%Get(p); call rhoAny%Get(rho) 

  call mesh%RunCellKernel(FDM_RhieChow_Correction_Kernel)

contains
  subroutine FDM_RhieChow_Correction_Kernel(iCell)
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
        &       lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        rrCell = mesh%CellToCell(rCellFace, rCell)
        llCell = mesh%CellToCell(lCellFace, lCell)
        if (gTruncErrorOrder >= 3) then
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
      ! Compute Rhie-Chow correction increment.
      ! ----------------------
      associate(dlInv => tau*lambda/mesh%dl(iCellFace))
        select case(gTruncErrorOrder)
          case(1:4)
            v(:,iCell) = v(:,iCell) - &
              &      ( dlInv*WFD4_C2(p(llCell), &
              &   1.0_dp/rho(lCell), p( lCell), &
              &   1.0_dp/rho(iCell), p( iCell), &
              &   1.0_dp/rho(rCell), p( rCell), &
              &                      p(rrCell)) )
        end select
      end associate
    end do

  end subroutine FDM_RhieChow_Correction_Kernel
end subroutine FDM_RhieChow_Correction

end module StormRuler_FDM_RhieChow
