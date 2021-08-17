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
module StormRuler_FDM_BCs

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Helpers, only: Flip, @{tSMapFunc$$@|@0, NUM_RANKS}@
use StormRuler_Mesh, only: tMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface FDM_ApplyBCs
#$do rank = 0, NUM_RANKS
  module procedure FDM_ApplyBCs$rank
#$end do
end interface FDM_ApplyBCs

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the third-order boundary conditions: αu + β∂u/∂n̅ = γ + f(̅x).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FDM_ApplyBCs$rank(mesh, iBCM, u, alpha, beta, gamma, f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: iBCM
  real(dp), intent(in) :: alpha, beta, gamma
  real(dp), intent(inout) :: u(@:,:)
  procedure(tSMapFunc$rank), optional :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer(ip) :: iBCMPtr
  ! ----------------------
  associate(bcmFirst => mesh%BCMs(iBCM), &
    &        bcmLast => mesh%BCMs(iBCM+1)-1, &
    &      bcmToCell => mesh%BCMToCell, &
    &  bcmToCellFace => mesh%BCMToCellFace, &
    &     cellToCell => mesh%CellToCell, &
    &    cellMDIndex => mesh%CellMDIndex, &
    &        mLambda => (0.5_dp*alpha - beta/mesh%dl), &
    &     pLambdaInv => 1.0_dp/(0.5_dp*alpha + beta/mesh%dl))
    ! ----------------------
    ! For each BC cell with the specific mark do:
    ! ----------------------
    !#omp parallel do schedule(static) &
    !#omp & default(none) private(iBCMPtr) shared(u, alpha, beta, gamma)
    do iBCMPtr = bcmFirst, bcmLast; block
      integer(ip) :: iCell, iBCCell, iBCCellFace, iGCell
      iBCCell = bcmToCell(iBCMPtr)
      iBCCellFace = bcmToCellFace(iBCMPtr)
      iCell = cellToCell(Flip(iBCCellFace), iBCCell)
      ! ----------------------
      ! Compute the FDM-approximate (second order) boundary conditions.
      ! ----------------------
      if (present(f)) then
        associate(x => 0.5_dp*( cellMDIndex(:,iCell) + &
          &                   cellMDIndex(:,iBCCell) ))
          u(@:,iBCCell) = pLambdaInv(iBCCellFace) * &
            & (gamma + f(x, u(@:,iCell)) - mLambda(iBCCellFace)*u(@:,iCell))
        end associate
      else
        u(@:,iBCCell) = pLambdaInv(iBCCellFace) * &
          & (gamma - mLambda(iBCCellFace)*u(@:,iCell))
      end if
      ! ----------------------
      ! Propagate the boundary condition towards the ghost cells.
      ! ----------------------
      iGCell = cellToCell(iBCCellFace, iBCCell)
      do while(iGCell /= 0)
        u(@:,iGCell) = u(@:,iBCCell)
        iGCell = cellToCell(iBCCellFace, iGCell)
      end do
    end block; end do
    !#omp end parallel do    
  end associate
end subroutine FDM_ApplyBCs$rank
#$end do

end module StormRuler_FDM_BCs
