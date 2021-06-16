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

module StormRuler_FDM_Operators

use StormRuler_Helpers
use StormRuler_Arithmetics
use StormRuler_Mesh
#$use 'StormRuler_Parameters.f90'

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

integer, parameter :: ACCURACY_ORDER = 4

private :: FD1_C2,FD1_C4,FD1_C6,FD1_C8

interface FDM_Gradient_Central
#$do rank = 0, NUM_RANKS-1
  module procedure FDM_Gradient_Central$rank
#$end do
end interface FDM_Gradient_Central

interface FDM_Divergence_Central
#$do rank = 0, NUM_RANKS-1
  module procedure FDM_Divergence_Central$rank
#$end do
end interface FDM_Divergence_Central

private :: FD1_F1,FD1_F2,FD1_F3,FD1_F4
private :: FD1_F5,FD1_F6,FD1_F7,FD1_F8

interface FDM_Gradient_Forward
#$do rank = 0, NUM_RANKS-1
  module procedure FDM_Gradient_Forward$rank
#$end do
end interface FDM_Gradient_Forward

interface FDM_Divergence_Backward
#$do rank = 0, NUM_RANKS-1
  module procedure FDM_Divergence_Backward$rank
#$end do
end interface FDM_Divergence_Backward

interface FDM_Convection_Central
#$do rank = 0, NUM_RANKS-1
  module procedure FDM_Convection_Central$rank
#$end do
end interface FDM_Convection_Central

private :: FD2_C2,FD2_C4,FD2_C6,FD2_C8

interface FDM_Laplacian_Central
#$do rank = 0, NUM_RANKS
  module procedure FDM_Laplacian_Central$rank
#$end do
end interface FDM_Laplacian_Central

interface FDM_LaplacianF_Central
#$do rank = 0, NUM_RANKS
  module procedure FDM_LaplacianF_Central$rank
#$end do
end interface FDM_LaplacianF_Central

interface FDM_Bilaplacian_Central
#$do rank = 0, NUM_RANKS
  module procedure FDM_Bilaplacian_Central$rank
#$end do
end interface FDM_Bilaplacian_Central

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_C2(u_l,u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_r,u_l
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C2 = 0.5_dp*(u_r - u_l)
end function FD1_C2

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_C4(u_ll,u_l,u_r,u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll,u_l,u_r,u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C4 = &
    & (-01.0_dp/12.0_dp)*u_rr + &
    & (+02.0_dp/03.0_dp)*u_r  + &
    & (-02.0_dp/03.0_dp)*u_l  + &
    & (+01.0_dp/12.0_dp)*u_ll
end function FD1_C4

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_C6(u_lll,u_ll,u_l,u_r,u_rr,u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll,u_ll,u_l,u_r,u_rr,u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C6 = &
    & (+01.0_dp/60.0_dp)*u_rrr + &
    & (-03.0_dp/20.0_dp)*u_rr  + &
    & (+03.0_dp/04.0_dp)*u_r   + &
    & (-03.0_dp/04.0_dp)*u_l   + &
    & (+03.0_dp/20.0_dp)*u_ll  + &
    & (-01.0_dp/60.0_dp)*u_lll
end function FD1_C6

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_C8(u_llll,u_lll,u_ll,u_l,u_r,u_rr,u_rrr,u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_llll,u_lll,u_ll,u_l,u_r,u_rr,u_rrr,u_rrrr
   ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C8 = &
    & (-001.0_dp/280.0_dp)*u_rrrr + &
    & (+004.0_dp/105.0_dp)*u_rrr  + &
    & (-001.0_dp/005.0_dp)*u_rr   + &
    & (+004.0_dp/005.0_dp)*u_r    + &
    & (-004.0_dp/005.0_dp)*u_l    + &
    & (+001.0_dp/005.0_dp)*u_ll   + &
    & (-004.0_dp/105.0_dp)*u_lll  + &
    & (+001.0_dp/280.0_dp)*u_llll
end function FD1_C8

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The central FDM-approximate gradient: v̅ ← v̅ - λ∇ₕu.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Gradient_Central$rank(mesh,vBar,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:),lambda
  real(dp), intent(inout) :: vBar(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell,iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell1, &
      &        dnInv=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do default(none) shared(u,vBar)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        rCell = cellToCell(+iCellFace,iCell)
        lCell = cellToCell(-iCellFace,iCell)
        if (ACCURACY_ORDER>=3) then
          rrCell = cellToCell(+iCellFace,rCell)
          llCell = cellToCell(-iCellFace,lCell)
          if (ACCURACY_ORDER>=5) then
            rrrCell = cellToCell(+iCellFace,rrCell)
            lllCell = cellToCell(-iCellFace,llCell)
            if (ACCURACY_ORDER>=7) then
              rrrrCell = cellToCell(+iCellFace,rrrCell)
              llllCell = cellToCell(-iCellFace,lllCell)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate gradient increment.
        ! ----------------------
        select case(ACCURACY_ORDER)
          ! ----------------------
          case(1:2)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &       ( dnInv(:,iCellFace).outer. &
              &               FD1_C2(u(@:,rCell), &
              &                      u(@:,lCell)) )
          ! ----------------------
          case(3:4)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &       ( dnInv(:,iCellFace).outer. &
              &              FD1_C4(u(@:,llCell), &
              &                     u(@:, lCell), &
              &                     u(@:, rCell), &
              &                     u(@:,rrCell)) )
          ! ----------------------
          case(5:6)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &       ( dnInv(:,iCellFace).outer. &
              &             FD1_C6(u(@:,lllCell), &
              &                    u(@:, llCell), &
              &                    u(@:,  lCell), &
              &                    u(@:,  rCell), &
              &                    u(@:, rrCell), &
              &                    u(@:,rrrCell)) )
          ! ----------------------
          case(7:8)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &       ( dnInv(:,iCellFace).outer. &
              &            FD1_C8(u(@:,llllCell), &
              &                   u(@:, lllCell), &
              &                   u(@:,  llCell), &
              &                   u(@:,   lCell), &
              &                   u(@:,   rCell), &
              &                   u(@:,  rrCell), &
              &                   u(@:, rrrCell), &
              &                   u(@:,rrrrCell)) )
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient_Central$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The central FDM-approximate divergence: v ← v - λ∇ₕ⋅u̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Divergence_Central$rank(mesh,v,lambda,uBar)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: uBar(:,@:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell,iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell1, &
      &        dnInv=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do default(none) shared(uBar,v)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        rCell = cellToCell(+iCellFace,iCell)
        lCell = cellToCell(-iCellFace,iCell)
        if (ACCURACY_ORDER>=3) then
          rrCell = cellToCell(+iCellFace,rCell)
          llCell = cellToCell(-iCellFace,lCell)
          if (ACCURACY_ORDER>=5) then
            rrrCell = cellToCell(+iCellFace,rrCell)
            lllCell = cellToCell(-iCellFace,llCell)
            if (ACCURACY_ORDER>=7) then
              rrrrCell = cellToCell(+iCellFace,rrrCell)
              llllCell = cellToCell(-iCellFace,lllCell)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate divergence increment.
        ! ----------------------
        select case(ACCURACY_ORDER)
          ! ----------------------
          case(1:2)
            v(@:,iCell) = v(@:,iCell) - &
              & ( dnInv(:,iCellFace).inner. &
              &    FD1_C2(uBar(:,@:,lCell), &
              &           uBar(:,@:,rCell)) )
          ! ----------------------
          case(3:4)
            v(@:,iCell) = v(@:,iCell) - &
              & ( dnInv(:,iCellFace).inner. &
              &   FD1_C4(uBar(:,@:,llCell), &
              &          uBar(:,@:, lCell), &
              &          uBar(:,@:, rCell), &
              &          uBar(:,@:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(@:,iCell) = v(@:,iCell) - &
              & ( dnInv(:,iCellFace).inner. &
              &   FD1_C6(uBar(:,@:,lllCell), &
              &          uBar(:,@:, llCell), &
              &          uBar(:,@:,  lCell), &
              &          uBar(:,@:,  rCell), &
              &          uBar(:,@:, rrCell), &
              &          uBar(:,@:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(@:,iCell) = v(@:,iCell) - &
              & ( dnInv(:,iCellFace).inner. &
              &   FD1_C8(uBar(:,@:,llllCell), &
              &          uBar(:,@:, lllCell), &
              &          uBar(:,@:,  llCell), &
              &          uBar(:,@:,   lCell), &
              &          uBar(:,@:,   rCell), &
              &          uBar(:,@:,  rrCell), &
              &          uBar(:,@:, rrrCell), &
              &          uBar(:,@:,rrrrCell)) )
          end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Divergence_Central$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! First order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F1(u,u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u,u_r
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F1 = u_r - u
end function FD1_F1

!! ----------------------------------------------------------------- !!
!! Second order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F2(u,u_r,u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u,u_r,u_rr
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F2 = &
    & (-1.5_dp)*u   + &
    & (+2.0_dp)*u_r + &
    & (-0.5_dp)*u_rr
end function FD1_F2

!! ----------------------------------------------------------------- !!
!! Third order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F3(u_l,u,u_r,u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l,u,u_r,u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F3 = &
    & (-1.0_dp/3.0_dp)*u_l + &
    & (-1.0_dp/2.0_dp)*u   + &
    &                  u_r + &
    & (-1.0_dp/6.0_dp)*u_rr
end function FD1_F3

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F4(u_l,u,u_r,u_rr,u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l,u,u_r,u_rr,u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F4 = &
    & (-01.0_dp/04.0_dp)*u_l  + &
    & (-05.0_dp/06.0_dp)*u    + &
    & (+03.0_dp/02.0_dp)*u_r  + &
    & (-01.0_dp/02.0_dp)*u_rr + &
    & (+01.0_dp/12.0_dp)*u_rrr
end function FD1_F4

!! ----------------------------------------------------------------- !!
!! Fifth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F5(u_ll,u_l,u,u_r,u_rr,u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll,u_l,u,u_r,u_rr,u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F5 = &
    & (+01.0_dp/20.0_dp)*u_ll + &
    & (-01.0_dp/02.0_dp)*u_l  + &
    & (-01.0_dp/03.0_dp)*u    + &
    &                    u_r  + &
    & (-01.0_dp/04.0_dp)*u_rr + &
    & (+01.0_dp/30.0_dp)*u_rrr
end function FD1_F5

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F6(u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F6 = &
    & (+01.0_dp/30.0_dp)*u_ll  + &
    & (-02.0_dp/05.0_dp)*u_l   + &
    & (-07.0_dp/12.0_dp)*u     + &
    & (+04.0_dp/03.0_dp)*u_r   + &
    & (-01.0_dp/02.0_dp)*u_rr  + &
    & (+02.0_dp/15.0_dp)*u_rrr + &
    & (-01.0_dp/60.0_dp)*u_rrrr
end function FD1_F6

!! ----------------------------------------------------------------- !!
!! Seventh order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F7(u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F7 = &
    & (-001.0_dp/105.0_dp)*u_lll + &
    & (+001.0_dp/010.0_dp)*u_ll  + &
    & (-003.0_dp/005.0_dp)*u_l   + &
    & (-001.0_dp/004.0_dp)*u     + &
    &                      u_r   + &
    & (-003.0_dp/010.0_dp)*u_rr  + &
    & (+001.0_dp/015.0_dp)*u_rrr + &
    & (-001.0_dp/140.0_dp)*u_rrrr
end function FD1_F7

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F8(u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr,u_rrrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr,u_rrrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F8 = &
    & (-001.0_dp/168.0_dp)*u_lll  + &
    & (+001.0_dp/014.0_dp)*u_ll   + &
    & (-001.0_dp/002.0_dp)*u_l    + &
    & (-009.0_dp/020.0_dp)*u      + &
    & (+005.0_dp/004.0_dp)*u_r    + &
    & (-001.0_dp/002.0_dp)*u_rr   + &
    & (+001.0_dp/006.0_dp)*u_rrr  + &
    & (-001.0_dp/028.0_dp)*u_rrrr + &
    & (+001.0_dp/280.0_dp)*u_rrrrr
end function FD1_F8

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The forward FDM-approximate gradient: v̅ ← v̅ - λ∇ₕu.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Gradient_Forward$rank(mesh,vBar,lambda,u,dirAll,dirFace,dirCell)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:),lambda
  real(dp), intent(inout) :: vBar(:,@:,:)
  integer(kind=1), intent(in), optional :: dirAll,dirFace(:),dirCell(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell,iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell1, &
      &        dnInv=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do default(none) shared(u,vBar) &
    !$omp shared(dirAll,dirFace,dirCell)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell,rrrrrCell
      integer :: lCell,llCell,lllCell,llllCell,lllllCell
      integer(kind=1) :: dir
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Determine FD direction.
        ! ----------------------
        dir = merge(dirAll,1_1,present(dirAll))
        dir = merge(dirFace(iCellFace),dir,present(dirFace))
        dir = merge(dirCell(iCellFace,iCell),dir,present(dirCell))
        dir = max(min(dir,+1_1),-1_1)
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        rCell = cellToCell(dir*iCellFace,iCell)
        lCell = cellToCell(dir*iCellFace,iCell)
        if (ACCURACY_ORDER>=2) then
          rrCell = cellToCell(dir*iCellFace,rCell)
          llCell = cellToCell(dir*iCellFace,lCell)
          if (ACCURACY_ORDER>=4) then
            rrrCell = cellToCell(dir*iCellFace,rrCell)
            lllCell = cellToCell(dir*iCellFace,llCell)
            if (ACCURACY_ORDER>=6) then
              rrrrCell = cellToCell(dir*iCellFace,rrrCell)
              llllCell = cellToCell(dir*iCellFace,lllCell)
              if (ACCURACY_ORDER>=8) then
                rrrrrCell = cellToCell(dir*iCellFace,rrrrCell)
                lllllCell = cellToCell(dir*iCellFace,llllCell)
              end if
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate gradient increment.
        ! ----------------------
        select case(ACCURACY_ORDER)
          case(1)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &   dir*( dnInv(:,iCellFace).outer. &
              &               FD1_F1(u(@:,iCell), &
              &                      u(@:,rCell)) )
          ! ----------------------
          case(2)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &   dir*( dnInv(:,iCellFace).outer. &
              &              FD1_F2(u(@:, iCell), &
              &                     u(@:, rCell), &
              &                     u(@:,rrCell)) )
          ! ----------------------
          case(3)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &   dir*( dnInv(:,iCellFace).outer. &
              &              FD1_F3(u(@:, lCell), &
              &                     u(@:, iCell), &
              &                     u(@:, rCell), &
              &                     u(@:,rrCell)) )
          ! ----------------------
          case(4)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &   dir*( dnInv(:,iCellFace).outer. &
              &             FD1_F4(u(@:,  lCell), &
              &                    u(@:,  iCell), &
              &                    u(@:,  rCell), &
              &                    u(@:, rrCell), &
              &                    u(@:,rrrCell)) )
          ! ----------------------
          case(5)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &   dir*( dnInv(:,iCellFace).outer. &
              &             FD1_F5(u(@:, llCell), &
              &                    u(@:,  lCell), &
              &                    u(@:,  iCell), &
              &                    u(@:,  rCell), &
              &                    u(@:, rrCell), &
              &                    u(@:,rrrCell)) )
          ! ----------------------
          case(6)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &   dir*( dnInv(:,iCellFace).outer. &
              &            FD1_F6(u(@:,  llCell), &
              &                   u(@:,   lCell), &
              &                   u(@:,   iCell), &
              &                   u(@:,   rCell), &
              &                   u(@:,  rrCell), &
              &                   u(@:, rrrCell), &
              &                   u(@:,rrrrCell)) )
          ! ----------------------
          case(7)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &   dir*( dnInv(:,iCellFace).outer. &
              &            FD1_F7(u(@:, lllCell), &
              &                   u(@:,  llCell), &
              &                   u(@:,   lCell), &
              &                   u(@:,   iCell), &
              &                   u(@:,   rCell), &
              &                   u(@:,  rrCell), &
              &                   u(@:, rrrCell), &
              &                   u(@:,rrrrCell)) )
          ! ----------------------
          case(8)
            vBar(:,@:,iCell) = vBar(:,@:,iCell) - &
              &   dir*( dnInv(:,iCellFace).outer. &
              &           FD1_F8(u(@:,  lllCell), &
              &                  u(@:,   llCell), &
              &                  u(@:,    lCell), &
              &                  u(@:,    iCell), &
              &                  u(@:,    rCell), &
              &                  u(@:,   rrCell), &
              &                  u(@:,  rrrCell), &
              &                  u(@:, rrrrCell), &
              &                  u(@:,rrrrrCell)) )
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient_Forward$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The backward FDM-approximate divergence: v ← v - λ∇ₕ⋅u̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Divergence_Backward$rank(mesh,v,lambda,uBar, &
  &                                     dirAll,dirFace,dirCell)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: uBar(:,@:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
  integer(kind=1), intent(in), optional :: dirAll,dirFace(:),dirCell(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell,iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell1, &
      &        dnInv=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do default(none) shared(uBar,v) &
    !$omp shared(dirAll,dirFace,dirCell)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell,rrrrrCell
      integer :: lCell,llCell,lllCell,llllCell,lllllCell
      integer(kind=1) :: dir
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Determine FD direction.
        ! ----------------------
        dir = merge(dirAll,1_1,present(dirAll))
        dir = merge(dirFace(iCellFace),dir,present(dirFace))
        dir = merge(dirCell(iCellFace,iCell),dir,present(dirCell))
        dir = max(min(dir,+1_1),-1_1)
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        rCell = cellToCell(dir*iCellFace,iCell)
        lCell = cellToCell(dir*iCellFace,iCell)
        if (ACCURACY_ORDER>=2) then
          rrCell = cellToCell(dir*iCellFace,rCell)
          llCell = cellToCell(dir*iCellFace,lCell)
          if (ACCURACY_ORDER>=4) then
            rrrCell = cellToCell(dir*iCellFace,rrCell)
            lllCell = cellToCell(dir*iCellFace,llCell)
            if (ACCURACY_ORDER>=6) then
              rrrrCell = cellToCell(dir*iCellFace,rrrCell)
              llllCell = cellToCell(dir*iCellFace,lllCell)
              if (ACCURACY_ORDER>=8) then
                rrrrrCell = cellToCell(dir*iCellFace,rrrrCell)
                lllllCell = cellToCell(dir*iCellFace,llllCell)
              end if
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate divergence increment.
        ! ----------------------
        select case(ACCURACY_ORDER)
          case(1)
            v(@:,iCell) = v(@:,iCell) + &
              & dir*( dnInv(:,iCellFace).inner. &
              &        FD1_F1(uBar(:,@:,iCell), &
              &               uBar(:,@:,lCell)) )
          ! ----------------------
          case(2)
            v(@:,iCell) = v(@:,iCell) + &
              & dir*( dnInv(:,iCellFace).inner. &
              &       FD1_F2(uBar(:,@:, iCell), &
              &              uBar(:,@:, lCell), &
              &              uBar(:,@:,llCell)) )
          ! ----------------------
          case(3)
            v(@:,iCell) = v(@:,iCell) + &
              & dir*( dnInv(:,iCellFace).inner. &
              &       FD1_F3(uBar(:,@:, rCell), &
              &              uBar(:,@:, iCell), &
              &              uBar(:,@:, lCell), &
              &              uBar(:,@:,llCell)) )
          ! ----------------------
          case(4)
            v(@:,iCell) = v(@:,iCell) + &
              & dir*( dnInv(:,iCellFace).inner. &
              &       FD1_F4(uBar(:,@:,  rCell), &
              &              uBar(:,@:,  iCell), &
              &              uBar(:,@:,  lCell), &
              &              uBar(:,@:, llCell), &
              &              uBar(:,@:,lllCell)) )
          ! ----------------------
          case(5)
            v(@:,iCell) = v(@:,iCell) + &
              & dir*( dnInv(:,iCellFace).inner. &
              &       FD1_F5(uBar(:,@:, rrCell), &
              &              uBar(:,@:,  rCell), &
              &              uBar(:,@:,  iCell), &
              &              uBar(:,@:,  lCell), &
              &              uBar(:,@:, llCell), &
              &              uBar(:,@:,lllCell)) )
          ! ----------------------
          case(6)
            v(@:,iCell) = v(@:,iCell) + &
              & dir*( dnInv(:,iCellFace).inner. &
              &       FD1_F6(uBar(:,@:,  rrCell), &
              &              uBar(:,@:,   rCell), &
              &              uBar(:,@:,   iCell), &
              &              uBar(:,@:,   lCell), &
              &              uBar(:,@:,  llCell), &
              &              uBar(:,@:, lllCell), &
              &              uBar(:,@:,llllCell)) )
          ! ----------------------
          case(7)
            v(@:,iCell) = v(@:,iCell) + &
              & dir*( dnInv(:,iCellFace).inner. &
              &       FD1_F7(uBar(:,@:, rrrCell), &
              &              uBar(:,@:,  rrCell), &
              &              uBar(:,@:,   rCell), &
              &              uBar(:,@:,   iCell), &
              &              uBar(:,@:,   lCell), &
              &              uBar(:,@:,  llCell), &
              &              uBar(:,@:, lllCell), &
              &              uBar(:,@:,llllCell)) )
          ! ----------------------
          case(8)
            v(@:,iCell) = v(@:,iCell) + &
              & dir*( dnInv(:,iCellFace).inner. &
              &       FD1_F8(uBar(:,@:,  rrrCell), &
              &              uBar(:,@:,   rrCell), &
              &              uBar(:,@:,    rCell), &
              &              uBar(:,@:,    iCell), &
              &              uBar(:,@:,    lCell), &
              &              uBar(:,@:,   llCell), &
              &              uBar(:,@:,  lllCell), &
              &              uBar(:,@:, llllCell), &
              &              uBar(:,@:,lllllCell)) )
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Divergence_Backward$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The central FDM-approximate convection: v ← v - λ∇ₕ⋅uw̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Convection_Central$rank(mesh,v,lambda,u,wBar)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:),wBar(:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: fBar(:,@:,:)
  allocate(fBar(size(wBar,dim=1), &
    & @{size(u,dim=$$+1)}@,size(wBar,dim=2)))
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  ! f̅ ← uw̅,
  ! v ← v - λ∇⋅f̅.
  ! ----------------------
  call Mul_Outer(mesh,fBar,wBar,u)
  call FDM_Divergence_Central(mesh,v,lambda,fBar)
end subroutine FDM_Convection_Central$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD2_C2(u_l,u,u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l,u,u_r
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C2 = u_l - 2.0_dp*u + u_r
end function FD2_C2

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD2_C4(u_ll,u_l,u,u_r,u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll,u_l,u,u_r,u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C4 = &
    & (-1.0_dp/12.0_dp)*u_ll + &
    & (+4.0_dp/03.0_dp)*u_l  + &
    & (-5.0_dp/02.0_dp)*u    + &
    & (+4.0_dp/03.0_dp)*u_r  + &
    & (-1.0_dp/12.0_dp)*u_rr
end function FD2_C4

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD2_C6(u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C6 = &
    & (+01.0_dp/90.0_dp)*u_lll + &
    & (-03.0_dp/20.0_dp)*u_ll  + &
    & (+03.0_dp/02.0_dp)*u_l   + &
    & (-49.0_dp/18.0_dp)*u     + &
    & (+03.0_dp/02.0_dp)*u_r   + &
    & (-03.0_dp/20.0_dp)*u_rr  + &
    & (+01.0_dp/90.0_dp)*u_rrr
end function FD2_C6

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD2_C8(u_llll,u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_llll,u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C8 = &
    & (-001.0_dp/560.0_dp)*u_llll + &
    & (+008.0_dp/315.0_dp)*u_lll  + &
    & (-001.0_dp/005.0_dp)*u_ll   + &
    & (+008.0_dp/005.0_dp)*u_l    + &
    & (-205.0_dp/072.0_dp)*u      + &
    & (+008.0_dp/005.0_dp)*u_r    + &
    & (-001.0_dp/005.0_dp)*u_rr   + &
    & (+008.0_dp/315.0_dp)*u_rrr  + &
    & (-001.0_dp/560.0_dp)*u_rrrr
end function FD2_C8

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The FDM-approximate Laplacian: v ← v + λΔₕu.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
#$do rank = 0, NUM_RANKS
subroutine FDM_Laplacian_Central$rank(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell,iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell1, &
      &     dxSqrInv=>lambda/(mesh%Dx**2))
    ! ----------------------
    !$omp parallel do default(none) shared(u,v)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        rCell = cellToCell(+iCellFace,iCell)
        lCell = cellToCell(-iCellFace,iCell)
        if (ACCURACY_ORDER>=3) then
          rrCell = cellToCell(+iCellFace,rCell)
          llCell = cellToCell(-iCellFace,lCell)
          if (ACCURACY_ORDER>=5) then
            rrrCell = cellToCell(+iCellFace,rrCell)
            lllCell = cellToCell(-iCellFace,llCell)
            if (ACCURACY_ORDER>=7) then
              rrrrCell = cellToCell(+iCellFace,rrrCell)
              llllCell = cellToCell(-iCellFace,lllCell)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate Laplacian increment.
        ! ----------------------
        select case(ACCURACY_ORDER)
          case(1:2)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dxSqrInv(iCellFace) * &
              &     FD2_C2(u(@:,lCell), &
              &            u(@:,iCell), &
              &            u(@:,rCell)) )
          ! ----------------------
          case(3:4)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dxSqrInv(iCellFace) * &
              &    FD2_C4(u(@:,llCell), &
              &           u(@:, lCell), &
              &           u(@:, iCell), &
              &           u(@:, rCell), &
              &           u(@:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dxSqrInv(iCellFace) * &
              &   FD2_C6(u(@:,lllCell), &
              &          u(@:, llCell), &
              &          u(@:,  lCell), &
              &          u(@:,  iCell), &
              &          u(@:,  rCell), &
              &          u(@:, rrCell), &
              &          u(@:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dxSqrInv(iCellFace) * &
              &   FD2_C8(u(@:,llllCell), &
              &          u(@:, lllCell), &
              &          u(@:,  llCell), &
              &          u(@:,   lCell), &
              &          u(@:,   iCell), &
              &          u(@:,   rCell), &
              &          u(@:,  rrCell), &
              &          u(@:, rrrCell), &
              &          u(@:,rrrrCell)) )
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Laplacian_Central$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate nonlinear Laplacian: v ← v + λΔₕf(u).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FDM_LaplacianF_Central$rank(mesh,v,lambda,f,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
  procedure(MathFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: w(@:,:)
  allocate(w,mold=u)
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  ! w ← f(u),
  ! v ← v + λΔw.
  ! ----------------------
  call ApplyFunc(mesh,w,u,f)
  call FDM_Laplacian_Central(mesh,v,lambda,w)
end subroutine FDM_LaplacianF_Central$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate Bilaplacian: v ← v + λ(Δₕ)²u.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FDM_Bilaplacian_Central$rank(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: w(@:,:)
  allocate(w,mold=u)
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  ! w ← 0,
  ! w ← w + Δw.
  ! v ← v + λΔw.
  ! ----------------------
  call Fill(mesh,w)
  call FDM_Laplacian_Central(mesh,w,1.0_dp,u)
  call FDM_Laplacian_Central(mesh,v,lambda,w)
end subroutine FDM_Bilaplacian_Central$rank
#$end do

end module StormRuler_FDM_Operators
