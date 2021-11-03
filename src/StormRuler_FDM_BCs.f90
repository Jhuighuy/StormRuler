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
use StormRuler_Helpers, only: Flip
#$for type_, _ in SCALAR_TYPES
!use StormRuler_Helpers, only: tSMapFunc$type_
#$end for
use StormRuler_Mesh, only: tMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface FDM_ApplyBCs
#$do rank = 0, NUM_RANKS-3
  module procedure FDM_ApplyBCs$rank
#$end do
end interface FDM_ApplyBCs

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the generalized third-order 
!! boundary conditions: 𝛼𝒖 + 𝛽∂𝒖/∂𝑛 = 𝛾 + 𝑓(𝑟).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-3
subroutine FDM_ApplyBCs$rank(mesh, iBCM, u, alpha, beta, gamma)!, f)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: iBCM
  real(dp), intent(in) :: alpha, beta, gamma
  real(dp), intent(inout) :: u(@:,:)
  !procedure(tSMapFuncR$rank), optional :: f
  
  integer(ip) :: iBCMPtr

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
    do iBCMPtr = bcmFirst, bcmLast; block
      integer(ip) :: iCell, iBCCell, iBCCellFace, iGCell
      iBCCell = bcmToCell(iBCMPtr)
      iBCCellFace = bcmToCellFace(iBCMPtr)
      iCell = cellToCell(Flip(iBCCellFace), iBCCell)

      ! ----------------------
      ! Compute the FDM-approximate (second order) boundary conditions.
      ! ----------------------
      !if (present(f)) then
      !  associate(x => 0.5_dp*( cellMDIndex(:,iCell) + &
      !    &                   cellMDIndex(:,iBCCell) ))
      !    u(@:,iBCCell) = pLambdaInv(iBCCellFace) * &
      !      & (gamma + f(x, u(@:,iCell)) - mLambda(iBCCellFace)*u(@:,iCell))
      !  end associate
      !else
        u(@:,iBCCell) = pLambdaInv(iBCCellFace) * &
          & (gamma - mLambda(iBCCellFace)*u(@:,iCell))
      !end if

      ! ----------------------
      ! Propagate the boundary condition towards the ghost cells.
      ! ----------------------
      iGCell = cellToCell(iBCCellFace, iBCCell)
      do while(iGCell /= 0)
        u(@:,iGCell) = u(@:,iBCCell)
        iGCell = cellToCell(iBCCellFace, iGCell)
      end do
    end block; end do

  end associate
end subroutine FDM_ApplyBCs$rank
#$end do

subroutine FDM_ApplyBCs_SlipWall(mesh, iBCM, v)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: iBCM
  real(dp), intent(inout) :: v(:,:)

  integer(ip) :: dim
  integer(ip) :: iBCMPtr
  integer(ip) :: iCell, iBCCell, iBCCellFace, iGCell

  do iBCMPtr = mesh%BCMs(iBCM), mesh%BCMs(iBCM+1)-1

    iBCCell = mesh%BCMToCell(iBCMPtr)
    iBCCellFace = mesh%BCMToCellFace(iBCMPtr)
    iCell = mesh%CellToCell(Flip(iBCCellFace), iBCCell)

    dim = (iBCCellFace - 1)/2 + 1

    v(:,iBCCell) = v(:,iCell)  
    v(dim,iBCCell) = -v(dim,iCell)  

  end do

end subroutine FDM_ApplyBCs_SlipWall

subroutine FDM_ApplyBCs_InOutLet(mesh, iBCM, v)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: iBCM
  real(dp), intent(inout) :: v(:,:)

  integer(ip) :: dim
  integer(ip) :: iBCMPtr
  integer(ip) :: iCell, iBCCell, iBCCellFace, iGCell

  real(dp) :: R, RR, RRR

  R = 0.0_dp

  do iBCMPtr = mesh%BCMs(iBCM), mesh%BCMs(iBCM+1)-1

    iBCCell = mesh%BCMToCell(iBCMPtr)
    iBCCellFace = mesh%BCMToCellFace(iBCMPtr)
    iCell = mesh%CellToCell(Flip(iBCCellFace), iBCCell)

    dim = (iBCCellFace - 1)/2 + 1

    R = max(R, mesh%CellCenter(1, iCell))

  end do

  R = R + 0.5_dp*mesh%dl(1)

  do iBCMPtr = mesh%BCMs(iBCM), mesh%BCMs(iBCM+1)-1

    iBCCell = mesh%BCMToCell(iBCMPtr)
    iBCCellFace = mesh%BCMToCellFace(iBCMPtr)
    iCell = mesh%CellToCell(Flip(iBCCellFace), iBCCell)

    RR = mesh%CellCenter(1, iCell)

    RRR = merge(0.1_dp, 10.0_dp, iBCM == 2)

    v(:,iBCCell) = 0.0_dp
    v(2,iBCCell) = -RRR*( R**2 - RR**2 )

  end do

end subroutine FDM_ApplyBCs_InOutLet

end module StormRuler_FDM_BCs
