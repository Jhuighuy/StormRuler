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

#$use 'StormRuler_Parameters.f90'
  
use StormRuler_Parameters, only: dp
use StormRuler_Helpers, only: Flip, ^{MathSpatialFunc$$^|^0,NUM_RANKS}^
use StormRuler_Mesh, only: Mesh2D

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface FDM_ApplyBCs
#$do rank = 0, NUM_RANKS-1
  module procedure FDM_ApplyBCs$rank
#$end do
end interface FDM_ApplyBCs

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the third-order boundary conditions: αu + β∂u/∂n̅ = γ + f(̅x).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine FDM_ApplyBCs$rank(mesh,iBCM,u,alpha,beta,gamma,f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  integer, intent(in) :: iBCM
  real(dp), intent(in) :: alpha,beta,gamma
  real(dp), intent(inout) :: u(^:,:)
  procedure(MathSpatialFunc$rank), optional :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iBCMPtr
  ! ----------------------
  associate(bcms=>mesh%BCMs, &
    &     bcmToCell=>mesh%BCMToCell, &
    & bcmToCellFace=>mesh%BCMToCellFace, &
    &    cellToCell=>mesh%CellToCell, &
    &   cellMDIndex=>mesh%CellMDIndex, &
    &       mLambda=>(0.5_dp*alpha - beta/(1+0*mesh%dl)), &
    &    pLambdaInv=>1.0_dp/(0.5_dp*alpha + beta/(1+0*mesh%dl)))
    ! ----------------------
    !$omp parallel do default(none) shared(u,alpha,beta,gamma)
    do iBCMPtr = bcms(iBCM), bcms(iBCM+1)-1; block
      integer :: iCell,iGCell,iBCCell,iBCCellFace
      iBCCell = bcmToCell(iBCMPtr)
      iBCCellFace = bcmToCellFace(iBCMPtr)
      iCell = cellToCell(Flip(iBCCellFace),iBCCell)
      ! ----------------------
      ! Compute the FDM-approximate (second order) boundary conditions.
      ! ----------------------
      if (present(f)) then
        associate(x=>0.5_dp*( cellMDIndex(:,iCell) + &
          &                   cellMDIndex(:,iBCCell) ))
          u(^:,iBCCell) = pLambdaInv(iBCCellFace) * &
            & (gamma + f(x,u) - mLambda(iBCCellFace)*u(^:,iCell))
        end associate
      else
        u(^:,iBCCell) = pLambdaInv(iBCCellFace) * &
          & (gamma - mLambda(iBCCellFace)*u(^:,iCell))
      end if
      ! ----------------------
      ! Propagate the boundary condition towards the ghost cells.
      ! ----------------------
      iGCell = cellToCell(iBCCellFace,iBCCell)
      do while(iGCell /= 0)
        u(^:,iGCell) = u(^:,iBCCell)
        iGCell = cellToCell(iBCCellFace,iGCell)
      end do
    end block; end do
    !$omp end parallel do    
  end associate
end subroutine FDM_ApplyBCs$rank
#$end do

end module StormRuler_FDM_BCs
