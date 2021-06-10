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
module StormRuler_FiniteDifferences

use StormRuler_Helpers
use StormRuler_Arithmetics
use StormRuler_Mesh
#@use 'StormRuler_Parameters.f90'

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

integer, parameter :: ACCURACY_ORDER = 4

private :: FD1_C2,FD1_C4,FD1_C6,FD1_C8

interface FDM_Gradient_Central
#@do rank = 0, NUM_RANKS-1
  module procedure FDM_Gradient_Central$rank
#@end do
end interface FDM_Gradient_Central

interface FDM_Divergence_Central
#@do rank = 0, NUM_RANKS-1
  module procedure FDM_Divergence_Central$rank
#@end do
end interface FDM_Divergence_Central

private :: FD1_F1,FD1_F2,FD1_F3,FD1_F4
private :: FD1_F5,FD1_F6,FD1_F7,FD1_F8

interface FDM_Gradient_Forward
#@do rank = 0, NUM_RANKS-1
  module procedure FDM_Gradient_Forward$rank
#@end do
end interface FDM_Gradient_Forward

interface FDM_Divergence_Backward
#@do rank = 0, NUM_RANKS-1
  module procedure FDM_Divergence_Backward$rank
#@end do
end interface FDM_Divergence_Backward

interface FDM_Convection_Central
#@do rank = 0, NUM_RANKS-1
  module procedure FDM_Convection_Central$rank
#@end do
end interface FDM_Convection_Central

private :: FD2_C2,FD2_C4,FD2_C6,FD2_C8

interface FDM_Laplacian_Central
#@do rank = 0, NUM_RANKS
  module procedure FDM_Laplacian_Central$rank
#@end do
end interface FDM_Laplacian_Central

interface FDM_LaplacianF_Central
#@do rank = 0, NUM_RANKS
  module procedure FDM_LaplacianF_Central$rank
#@end do
end interface FDM_LaplacianF_Central

interface FDM_Bilaplacian_Central
#@do rank = 0, NUM_RANKS
  module procedure FDM_Bilaplacian_Central$rank
#@end do
end interface FDM_Bilaplacian_Central

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Second order accracy central undivided finite difference.
elemental function FD1_C2(u_r,u_l) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_r
  real(dp), intent(in) :: u_l
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = 0.5_dp*(u_r - u_l)
end function FD1_C2
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Fourth order accracy central undivided finite difference.
elemental function FD1_C4(u_r,u_rr, &
  &                       u_l,u_ll) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_r,u_rr
  real(dp), intent(in) :: u_l,u_ll
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-1.0_dp/12.0_dp)*u_rr + &
    &  (+2.0_dp/03.0_dp)*u_r  + &
    &  (-2.0_dp/03.0_dp)*u_l  + &
    &  (+1.0_dp/12.0_dp)*u_ll
end function FD1_C4
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Sixth order accracy central undivided finite difference.
elemental function FD1_C6(u_r,u_rr,u_rrr, &
  &                       u_l,u_ll,u_lll) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_r,u_rr,u_rrr
  real(dp), intent(in) :: u_l,u_ll,u_lll
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (+01.0_dp/60.0_dp)*u_rrr + &
    &  (-03.0_dp/20.0_dp)*u_rr  + &
    &  (+03.0_dp/04.0_dp)*u_r   + &
    &  (-03.0_dp/04.0_dp)*u_l   + &
    &  (+03.0_dp/20.0_dp)*u_ll  + &
    &  (-01.0_dp/60.0_dp)*u_lll
end function FD1_C6
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Eighth order accracy central undivided finite difference.
elemental function FD1_C8(u_r,u_rr,u_rrr,u_rrrr, &
  &                       u_l,u_ll,u_lll,u_llll) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_r,u_rr,u_rrr,u_rrrr
  real(dp), intent(in) :: u_l,u_ll,u_lll,u_llll
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-001.0_dp/280.0_dp)*u_rrrr + &
    &  (+004.0_dp/105.0_dp)*u_rrr  + &
    &  (-001.0_dp/005.0_dp)*u_rr   + &
    &  (+004.0_dp/005.0_dp)*u_r    + &
    &  (-004.0_dp/005.0_dp)*u_l    + &
    &  (+001.0_dp/005.0_dp)*u_ll   + &
    &  (-004.0_dp/105.0_dp)*u_lll  + &
    &  (+001.0_dp/280.0_dp)*u_llll
end function FD1_C8
!! ----------------------------------------------------------------- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The central FDM-approximate gradient: v̅ ← v̅ - λ∇ₕu.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Gradient_Central$rank(mesh,vBar,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:),lambda
  real(dp), intent(inout) :: vBar(:,@:,:)
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
      &      nOverDx=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do default(none) shared(u,vBar)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(+iCellFace,iCell)
        lCell = cellToCell(-iCellFace,iCell)
        if (ACCURACY_ORDER >= 3) then
          rrCell = cellToCell(+iCellFace,rCell)
          llCell = cellToCell(-iCellFace,lCell)
          if (ACCURACY_ORDER >= 5) then
            rrrCell = cellToCell(+iCellFace,rrCell)
            lllCell = cellToCell(-iCellFace,llCell)
            if (ACCURACY_ORDER >= 7) then
              rrrrCell = cellToCell(+iCellFace,rrrCell)
              llllCell = cellToCell(-iCellFace,lllCell)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate gradient increment.
        select case(ACCURACY_ORDER)
        case(1:2)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_C2(u(@:,rCell),u(@:,lCell)))
        case(3:4)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_C4(u(@:,rCell),u(@:,rrCell), &
            &         u(@:,lCell),u(@:,llCell)))
        case(5:6)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_C6(u(@:,rCell),u(@:,rrCell),u(@:,rrrCell), &
            &         u(@:,lCell),u(@:,llCell),u(@:,lllCell)))
        case(7:8)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_C8(u(@:,rCell),u(@:,rrCell),u(@:,rrrCell),u(@:,rrrrCell), &
            &         u(@:,lCell),u(@:,llCell),u(@:,lllCell),u(@:,llllCell)))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient_Central$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The central FDM-approximate divergence: v ← v - λ∇ₕ⋅u̅.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Divergence_Central$rank(mesh,v,lambda,uBar)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: uBar(:,@:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
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
      &      nOverDx=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do default(none) shared(uBar,v)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(+iCellFace,iCell)
        lCell = cellToCell(-iCellFace,iCell)
        if (ACCURACY_ORDER >= 3) then
          rrCell = cellToCell(+iCellFace,rCell)
          llCell = cellToCell(-iCellFace,lCell)
          if (ACCURACY_ORDER >= 5) then
            rrrCell = cellToCell(+iCellFace,rrCell)
            lllCell = cellToCell(-iCellFace,llCell)
            if (ACCURACY_ORDER >= 7) then
              rrrrCell = cellToCell(+iCellFace,rrrCell)
              llllCell = cellToCell(-iCellFace,lllCell)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate divergence increment.
        select case(ACCURACY_ORDER)
        case(1:2)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_C2(uBar(:,@:,rCell),uBar(:,@:,lCell)))
        case(3:4)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_C4(uBar(:,@:,rCell),uBar(:,@:,rrCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell)))
        case(5:6)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_C6(uBar(:,@:,rCell),uBar(:,@:,rrCell),uBar(:,@:,rrrCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell),uBar(:,@:,lllCell)))
        case(7:8)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_C8(uBar(:,@:,rCell),uBar(:,@:,rrCell),uBar(:,@:,rrrCell),uBar(:,@:,rrrrCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell),uBar(:,@:,lllCell),uBar(:,@:,llllCell)))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Divergence_Central$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! First order accracy forward undivided finite difference.
elemental function FD1_F1(u,u_r) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u,u_r
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = u_r - u
end function FD1_F1
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Second order accracy forward undivided finite difference.
elemental function FD1_F2(u,u_r,u_rr) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u,u_r,u_rr
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-1.5_dp)*u + (+2.0_dp)*u_r + (-0.5_dp)*u_rr
end function FD1_F2
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Third order accracy forward undivided finite difference.
elemental function FD1_F3(u_l,u,u_r,u_rr) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l,u,u_r,u_rr
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-1.0_dp/3.0_dp)*u_l + (-0.5_dp)*u + u_r + (-1.0_dp/6.0_dp)*u_rr
end function FD1_F3
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Fourth order accracy forward undivided finite difference.
elemental function FD1_F4(u_l,u,u_r,u_rr,u_rrr) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l,u,u_r,u_rr,u_rrr
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-0.25_dp)*u_l + (-5.0_dp/6.0_dp)*u + &
    & (+1.5_dp)*u_r + (-0.5_dp)*u_rr + (1.0_dp/12.0_dp)*u_rrr
end function FD1_F4
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Fifth order accracy forward undivided finite difference.
elemental function FD1_F5(u_ll,u_l,u,u_r,u_rr,u_rrr) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll,u_l,u,u_r,u_rr,u_rrr
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-0.25_dp)*u_l + (-5.0_dp/6.0_dp)*u + &
    & (+1.5_dp)*u_r + (-0.5_dp)*u_rr + (1.0_dp/12.0_dp)*u_rrr
end function FD1_F5
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Sixth order accracy forward undivided finite difference.
elemental function FD1_F6(u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-0.25_dp)*u_l + (-5.0_dp/6.0_dp)*u + &
    & (+1.5_dp)*u_r + (-0.5_dp)*u_rr + (1.0_dp/12.0_dp)*u_rrr
end function FD1_F6
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Seventh order accracy forward undivided finite difference.
elemental function FD1_F7(u_lll,u_ll,u_l, &
  &                       u,u_r,u_rr,u_rrr,u_rrrr) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-0.25_dp)*u_l + (-5.0_dp/6.0_dp)*u + &
    & (+1.5_dp)*u_r + (-0.5_dp)*u_rr + (1.0_dp/12.0_dp)*u_rrr
end function FD1_F7
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Eighth order accracy forward undivided finite difference.
elemental function FD1_F8(u_lll,u_ll,u_l, &
  &                       u,u_r,u_rr,u_rrr,u_rrrr,u_rrrrr) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll,u_ll,u_l,u,u_r,u_rr,u_rrr,u_rrrr,u_rrrrr
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-0.25_dp)*u_l + (-5.0_dp/6.0_dp)*u + &
    & (+1.5_dp)*u_r + (-0.5_dp)*u_rr + (1.0_dp/12.0_dp)*u_rrr
end function FD1_F8
!! ----------------------------------------------------------------- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The forward FDM-approximate gradient: v̅ ← v̅ - λ∇ₕu.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Gradient_Forward$rank(mesh,vBar,lambda,u,dir)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  integer, intent(in), optional :: dir(:)
  real(dp), intent(in) :: u(@:,:),lambda
  real(dp), intent(inout) :: vBar(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell,iCellFace
  integer :: d(1:mesh%NumCellFaces/2)
  d(:) = 1
  ! ----------------------
  ! Fast exit in case λ=0.
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell1, &
      &      nOverDx=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do default(none) shared(d,u,vBar)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell,rrrrrCell
      integer :: lCell,llCell,lllCell,llllCell,lllllCell
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(+iCellFace,iCell)
        lCell = cellToCell(-iCellFace,iCell)
        if (ACCURACY_ORDER >= 2) then
          rrCell = cellToCell(+iCellFace,rCell)
          llCell = cellToCell(-iCellFace,lCell)
          if (ACCURACY_ORDER >= 4) then
            rrrCell = cellToCell(+iCellFace,rrCell)
            lllCell = cellToCell(-iCellFace,llCell)
            if (ACCURACY_ORDER >= 6) then
              rrrrCell = cellToCell(+iCellFace,rrrCell)
              llllCell = cellToCell(-iCellFace,lllCell)
              if (ACCURACY_ORDER >= 8) then
                rrrrrCell = cellToCell(+iCellFace,rrrrCell)
                lllllCell = cellToCell(-iCellFace,llllCell)
              end if
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate gradient increment.
        select case(d(iCellFace)*ACCURACY_ORDER)
        case(+1)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F1(u(@:,iCell),u(@:,rCell)))
        case(-1)
          vBar(:,@:,iCell) += &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F1(u(@:,iCell),u(@:,lCell)))
        case(+2)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F2(u(@:,iCell), &
            &         u(@:,rCell),u(@:,rrCell)))
        case(-2)
          vBar(:,@:,iCell) += &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F2(u(@:,iCell), &
            &         u(@:,lCell),u(@:,llCell)))
        case(+3)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F3(u(@:,lCell),u(@:,iCell), &
            &         u(@:,rCell),u(@:,rrCell)))
        case(-3)
          vBar(:,@:,iCell) += &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F3(u(@:,rCell),u(@:,iCell), &
            &         u(@:,lCell),u(@:,llCell)))
        case(+4)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F4(u(@:,lCell),u(@:,iCell), &
            &         u(@:,rCell),u(@:,rrCell),u(@:,rrrCell)))
        case(-4)
          vBar(:,@:,iCell) += &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F4(u(@:,rCell),u(@:,iCell), &
            &         u(@:,lCell),u(@:,llCell),u(@:,lllCell)))
        case(+5)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F5(u(@:,llCell),u(@:,lCell),u(@:,iCell), &
            &         u(@:,rCell),u(@:,rrCell),u(@:,rrrCell)))
        case(-5)
          vBar(:,@:,iCell) += &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F5(u(@:,rrCell),u(@:,rCell),u(@:,iCell), &
            &         u(@:,lCell),u(@:,llCell),u(@:,lllCell)))
        case(+6)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F6(u(@:,llCell),u(@:,lCell),u(@:,iCell), &
            &         u(@:,rCell),u(@:,rrCell),u(@:,rrrCell),u(@:,rrrrCell)))
        case(-6)
          vBar(:,@:,iCell) += &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F6(u(@:,rrCell),u(@:,rCell),u(@:,iCell), &
            &         u(@:,lCell),u(@:,llCell),u(@:,lllCell),u(@:,llllCell)))
        case(+7)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F7(u(@:,lllCell),u(@:,llCell),u(@:,lCell),u(@:,iCell), &
            &         u(@:,rCell),u(@:,rrCell),u(@:,rrrCell),u(@:,rrrrCell)))
        case(-7)
          vBar(:,@:,iCell) += &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F7(u(@:,rrrCell),u(@:,rrCell),u(@:,rCell),u(@:,iCell), &
            &         u(@:,lCell),u(@:,llCell),u(@:,lllCell),u(@:,llllCell)))
        case(+8)
          vBar(:,@:,iCell) -= &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F8(u(@:,lllCell),u(@:,llCell),u(@:,lCell),u(@:,iCell), &
            &         u(@:,rCell),u(@:,rrCell),u(@:,rrrCell),u(@:,rrrrCell),u(@:,rrrrrCell)))
        case(-8)
          vBar(:,@:,iCell) += &
            & (nOverDx(:,iCellFace).outer.&
            &  FD1_F8(u(@:,rrrCell),u(@:,rrCell),u(@:,rCell),u(@:,iCell), &
            &         u(@:,lCell),u(@:,llCell),u(@:,lllCell),u(@:,llllCell),u(@:,lllllCell)))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient_Forward$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The backward FDM-approximate divergence: v ← v - λ∇ₕ⋅u̅.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Divergence_Backward$rank(mesh,v,lambda,uBar)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: uBar(:,@:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell,iCellFace
  integer :: d(1:mesh%NumCellFaces/2)
  d(:) = 1
  ! ----------------------
  ! Fast exit in case λ=0.
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell1, &
      &      nOverDx=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do default(none) shared(d,uBar,v)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell,rrrrrCell
      integer :: lCell,llCell,lllCell,llllCell,lllllCell
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(+iCellFace,iCell)
        lCell = cellToCell(-iCellFace,iCell)
        if (ACCURACY_ORDER >= 2) then
          rrCell = cellToCell(+iCellFace,rCell)
          llCell = cellToCell(-iCellFace,lCell)
          if (ACCURACY_ORDER >= 4) then
            rrrCell = cellToCell(+iCellFace,rrCell)
            lllCell = cellToCell(-iCellFace,llCell)
            if (ACCURACY_ORDER >= 6) then
              rrrrCell = cellToCell(+iCellFace,rrrCell)
              llllCell = cellToCell(-iCellFace,lllCell)
              if (ACCURACY_ORDER >= 8) then
                rrrrrCell = cellToCell(+iCellFace,rrrrCell)
                lllllCell = cellToCell(-iCellFace,llllCell)
              end if
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate divergence increment.
        select case(d(iCellFace)*ACCURACY_ORDER)
        case(-1)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F1(uBar(:,@:,iCell),uBar(:,@:,rCell)))
        case(+1)
          v(@:,iCell) += &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F1(uBar(:,@:,iCell),uBar(:,@:,lCell)))
        case(-2)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F2(uBar(:,@:,iCell), &
            &         uBar(:,@:,rCell),uBar(:,@:,rrCell)))
        case(+2)
          v(@:,iCell) += &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F2(uBar(:,@:,iCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell)))
        case(-3)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F3(uBar(:,@:,lCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,rCell),uBar(:,@:,rrCell)))
        case(+3)
          v(@:,iCell) += &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F3(uBar(:,@:,rCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell)))
        case(-4)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F4(uBar(:,@:,lCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,rCell),uBar(:,@:,rrCell),uBar(:,@:,rrrCell)))
        case(+4)
          v(@:,iCell) += &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F4(uBar(:,@:,rCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell),uBar(:,@:,lllCell)))
        case(-5)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F5(uBar(:,@:,llCell),uBar(:,@:,lCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,rCell),uBar(:,@:,rrCell),uBar(:,@:,rrrCell)))
        case(+5)
          v(@:,iCell) += &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F5(uBar(:,@:,rrCell),uBar(:,@:,rCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell),uBar(:,@:,lllCell)))
        case(-6)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F6(uBar(:,@:,llCell),uBar(:,@:,lCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,rCell),uBar(:,@:,rrCell),uBar(:,@:,rrrCell),uBar(:,@:,rrrrCell)))
        case(+6)
          v(@:,iCell) += &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F6(uBar(:,@:,rrCell),uBar(:,@:,rCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell),uBar(:,@:,lllCell),uBar(:,@:,llllCell)))
        case(-7)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F7(uBar(:,@:,lllCell),uBar(:,@:,llCell),uBar(:,@:,lCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,rCell),uBar(:,@:,rrCell),uBar(:,@:,rrrCell),uBar(:,@:,rrrrCell)))
        case(+7)
          v(@:,iCell) += &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F7(uBar(:,@:,rrrCell),uBar(:,@:,rrCell),uBar(:,@:,rCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell),uBar(:,@:,lllCell),uBar(:,@:,llllCell)))
        case(-8)
          v(@:,iCell) -= &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F8(uBar(:,@:,lllCell),uBar(:,@:,llCell),uBar(:,@:,lCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,rCell),uBar(:,@:,rrCell),uBar(:,@:,rrrCell),uBar(:,@:,rrrrCell), &
            &                                                               uBar(:,@:,rrrrrCell)))
        case(+8)
          v(@:,iCell) += &
            & (nOverDx(:,iCellFace).inner.&
            &  FD1_F8(uBar(:,@:,rrrCell),uBar(:,@:,rrCell),uBar(:,@:,rCell),uBar(:,@:,iCell), &
            &         uBar(:,@:,lCell),uBar(:,@:,llCell),uBar(:,@:,lllCell),uBar(:,@:,llllCell), &
            &                                                               uBar(:,@:,lllllCell)))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Divergence_Backward$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The central FDM-approximate convection: v ← v - λ∇ₕ⋅uw̅.
#@do rank = 0, NUM_RANKS-1
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
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  ! f̅ ← uw̅,
  ! v ← v - λ∇⋅f̅.
  call Mul_Outer(mesh,fBar,wBar,u)
  call FDM_Divergence_Central(mesh,v,lambda,fBar)
end subroutine FDM_Convection_Central$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Second order accracy central undivided second finite difference.
elemental function FD2_C2(u,u_r,u_l) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u
  real(dp), intent(in) :: u_r
  real(dp), intent(in) :: u_l
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = u_r - 2.0_dp*u + u_l
end function FD2_C2
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Fourth order accracy central undivided second finite difference.
elemental function FD2_C4(u, &
  &                       u_r,u_rr, &
  &                       u_l,u_ll) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u
  real(dp), intent(in) :: u_r,u_rr
  real(dp), intent(in) :: u_l,u_ll
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-1.0_dp/12.0_dp)*u_rr + &
    &  (+4.0_dp/03.0_dp)*u_r  + &
    &  (-5.0_dp/02.0_dp)*u    + &
    &  (+4.0_dp/03.0_dp)*u_l  + &
    &  (-1.0_dp/12.0_dp)*u_ll
end function FD2_C4
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Sixth order accracy central undivided second finite difference.
elemental function FD2_C6(u, &
  &                       u_r,u_rr,u_rrr, &
  &                       u_l,u_ll,u_lll) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u
  real(dp), intent(in) :: u_r,u_rr,u_rrr
  real(dp), intent(in) :: u_l,u_ll,u_lll
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (+01.0_dp/90.0_dp)*u_rrr + &
    &  (-03.0_dp/20.0_dp)*u_rr  + &
    &  (+03.0_dp/02.0_dp)*u_r   + &
    &  (-49.0_dp/18.0_dp)*u     + &
    &  (+03.0_dp/02.0_dp)*u_l   + &
    &  (-03.0_dp/20.0_dp)*u_ll  + &
    &  (+01.0_dp/90.0_dp)*u_lll
end function FD2_C6
!! ----------------------------------------------------------------- !!

!! ----------------------------------------------------------------- !!
!! Eighth order accracy central undivided second finite difference.
elemental function FD2_C8(u, &
  &                       u_r,u_rr,u_rrr,u_rrrr, &
  &                       u_l,u_ll,u_lll,u_llll) result(du)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u
  real(dp), intent(in) :: u_r,u_rr,u_rrr,u_rrrr
  real(dp), intent(in) :: u_l,u_ll,u_lll,u_llll
  real(dp) :: du
  ! >>>>>>>>>>>>>>>>>>>>>>
  du = (-001.0_dp/560.0_dp)*u_rrrr + &
    &  (+008.0_dp/315.0_dp)*u_rrr  + &
    &  (-001.0_dp/005.0_dp)*u_rr   + &
    &  (+008.0_dp/005.0_dp)*u_r    + &
    &  (-205.0_dp/072.0_dp)*u      + &
    &  (+008.0_dp/005.0_dp)*u_l    + &
    &  (-001.0_dp/005.0_dp)*u_ll   + &
    &  (+008.0_dp/315.0_dp)*u_lll  + &
    &  (-001.0_dp/560.0_dp)*u_llll
end function FD2_C8
!! ----------------------------------------------------------------- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The FDM-approximate Laplacian: v ← v + λΔₕu.
#@do rank = 0, NUM_RANKS
subroutine FDM_Laplacian_Central$rank(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:),lambda
  real(dp), intent(inout) :: v(@:,:)
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
      &     invDxSqr=>lambda/mesh%Dx**2)
    ! ----------------------
    !$omp parallel do default(none) shared(u,v)
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces/2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(+iCellFace,iCell)
        lCell = cellToCell(-iCellFace,iCell)
        if (ACCURACY_ORDER >= 3) then
          rrCell = cellToCell(+iCellFace,rCell)
          llCell = cellToCell(-iCellFace,lCell)
          if (ACCURACY_ORDER >= 5) then
            rrrCell = cellToCell(+iCellFace,rrCell)
            lllCell = cellToCell(-iCellFace,llCell)
            if (ACCURACY_ORDER >= 7) then
              rrrrCell = cellToCell(+iCellFace,rrrCell)
              llllCell = cellToCell(-iCellFace,lllCell)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate Laplacian increment.
        select case(ACCURACY_ORDER)
        case(1:2)
          v(@:,iCell) += &
            & invDxSqr(iCellFace)* &
            & FD2_C2(u(@:,iCell),  &
            &        u(@:,rCell),u(@:,lCell))
        case(3:4)
          v(@:,iCell) += &
            & invDxSqr(iCellFace)* &
            & FD2_C4(u(@:,iCell),  &
            &        u(@:,rCell),u(@:,rrCell), &
            &        u(@:,lCell),u(@:,llCell))
        case(5:6)
          v(@:,iCell) += &
            & invDxSqr(iCellFace)* &
            & FD2_C6(u(@:,iCell),  &
            &        u(@:,rCell),u(@:,rrCell),u(@:,rrrCell), &
            &        u(@:,lCell),u(@:,llCell),u(@:,lllCell))
        case(7:8)
          v(@:,iCell) += &
            & invDxSqr(iCellFace)* &
            & FD2_C8(u(@:,iCell),  &
            &        u(@:,rCell),u(@:,rrCell),u(@:,rrrCell),u(@:,rrrrCell), &
            &        u(@:,lCell),u(@:,llCell),u(@:,lllCell),u(@:,llllCell))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Laplacian_Central$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate nonlinear Laplacian: v ← v + λΔₕf(u).
#@do rank = 0, NUM_RANKS
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
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  ! w ← f(u),
  ! v ← v + λΔw.
  call ApplyFunc(mesh,w,u,f)
  call FDM_Laplacian_Central(mesh,v,lambda,w)
end subroutine FDM_LaplacianF_Central$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate Bilaplacian: v ← v + λ(Δₕ)²u.
#@do rank = 0, NUM_RANKS
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
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  ! w ← 0,
  ! w ← w + Δw.
  ! v ← v + λΔw.
  call Fill(mesh,w)
  call FDM_Laplacian_Central(mesh,w,1.0_dp,u)
  call FDM_Laplacian_Central(mesh,v,lambda,w)
end subroutine FDM_Bilaplacian_Central$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

end module StormRuler_FiniteDifferences
