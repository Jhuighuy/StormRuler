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

#$use 'StormRuler_Parameters.f90'

use StormRuler_Parameters, only: dp
use StormRuler_Helpers, only: Flip, SafeInverse, &
  & @{tMapFunc$$@|@0, NUM_RANKS}@, operator(.inner.), operator(.outer.)
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: Fill, Mul_Outer, FuncProd

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

integer, parameter :: gAccuracyOrder = 2

private :: FD1_C2, FD1_C4, FD1_C6, FD1_C8

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

private :: FD1_F1, FD1_F2, FD1_F3, FD1_F4
private :: FD1_F5, FD1_F6, FD1_F7, FD1_F8

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

private :: FD2_C2, FD2_C4, FD2_C6, FD2_C8

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

private :: WFD2_C2, WFD2_C4, WFD2_C6, WFD2_C8

interface FDM_DivWGrad_Central
#$do rank = 0, NUM_RANKS
  module procedure FDM_DivWGrad_Central$rank
#$end do
end interface FDM_DivWGrad_Central

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_C2(u_l, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_r, u_l
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C2 = 0.5_dp*(u_r - u_l)
end function FD1_C2

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_C4(u_ll, u_l, u_r, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll, u_l, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C4 = &
    & ( (-01.0_dp/12.0_dp)*u_rr + &
    &   (+02.0_dp/03.0_dp)*u_r  + &
    &   (-02.0_dp/03.0_dp)*u_l  + &
    &   (+01.0_dp/12.0_dp)*u_ll )
end function FD1_C4

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_C6(u_lll, u_ll, u_l, u_r, u_rr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll, u_ll, u_l, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C6 = &
    & ( (+01.0_dp/60.0_dp)*u_rrr + &
    &   (-03.0_dp/20.0_dp)*u_rr  + &
    &   (+03.0_dp/04.0_dp)*u_r   + &
    &   (-03.0_dp/04.0_dp)*u_l   + &
    &   (+03.0_dp/20.0_dp)*u_ll  + &
    &   (-01.0_dp/60.0_dp)*u_lll )
end function FD1_C6

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_C8(u_llll, u_lll, u_ll, u_l, u_r, u_rr, u_rrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_llll, u_lll, u_ll, u_l, u_r, u_rr, u_rrr, u_rrrr
   ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C8 = &
    & ( (-001.0_dp/280.0_dp)*u_rrrr + &
    &   (+004.0_dp/105.0_dp)*u_rrr  + &
    &   (-001.0_dp/005.0_dp)*u_rr   + &
    &   (+004.0_dp/005.0_dp)*u_r    + &
    &   (-004.0_dp/005.0_dp)*u_l    + &
    &   (+001.0_dp/005.0_dp)*u_ll   + &
    &   (-004.0_dp/105.0_dp)*u_lll  + &
    &   (+001.0_dp/280.0_dp)*u_llll )
end function FD1_C8

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate gradient: v̅ ← v̅ - λ∇u.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Gradient_Central$rank(mesh, v_bar, lambda, u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), pointer :: u(@:,:), v_bar(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells => mesh%NumCells, &
    &   numCellFaces => mesh%NumCellFaces, &
    &     cellToCell => mesh%CellToCell, &
    &         dr_inv => lambda*SafeInverse(mesh%dr))
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    !$omp parallel do schedule(static) &
    !$omp & default(none) private(iCell, iCellFace) shared(u, v_bar)
    do iCell = 1, numCells; block
      integer :: rCell, rrCell, rrrCell, rrrrCell
      integer :: lCell, llCell, lllCell, llllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        associate(rCellFace => iCellFace, &
          &       lCellFace => Flip(iCellFace))
          rCell = cellToCell(rCellFace, iCell)
          lCell = cellToCell(lCellFace, iCell)
          if (gAccuracyOrder >= 3) then
            rrCell = cellToCell(rCellFace, rCell)
            llCell = cellToCell(lCellFace, lCell)
            if (gAccuracyOrder >= 5) then
              rrrCell = cellToCell(rCellFace, rrCell)
              lllCell = cellToCell(lCellFace, llCell)
              if (gAccuracyOrder >= 7) then
                rrrrCell = cellToCell(rCellFace, rrrCell)
                llllCell = cellToCell(lCellFace, lllCell)
              end if
            end if
          end if
        end associate
        ! ----------------------
        ! Compute FDM-approximate gradient increment.
        ! ----------------------
        select case(gAccuracyOrder)
          case(1:2)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &          ( dr_inv(:,iCellFace).outer. &
              &                     FD1_C2(u(@:,lCell), &
              &                            u(@:,rCell)) )
          ! ----------------------
          case(3:4)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &          ( dr_inv(:,iCellFace).outer. &
              &                    FD1_C4(u(@:,llCell), &
              &                           u(@:, lCell), &
              &                           u(@:, rCell), &
              &                           u(@:,rrCell)) )
          ! ----------------------
          case(5:6)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &          ( dr_inv(:,iCellFace).outer. &
              &                   FD1_C6(u(@:,lllCell), &
              &                          u(@:, llCell), &
              &                          u(@:,  lCell), &
              &                          u(@:,  rCell), &
              &                          u(@:, rrCell), &
              &                          u(@:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &          ( dr_inv(:,iCellFace).outer. &
              &                   FD1_C8(u(@:,llllCell), &
              &                          u(@:, lllCell), &
              &                          u(@:,  llCell), &
              &                          u(@:,   lCell), &
              &                          u(@:,   rCell), &
              &                          u(@:,  rrCell), &
              &                          u(@:, rrrCell), &
              &                          u(@:,rrrrCell)) )
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient_Central$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate divergence: v ← v - λ∇⋅u̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Divergence_Central$rank(mesh, v, lambda, u_bar)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), pointer :: u_bar(:,@:,:), v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells => mesh%NumCells, &
    &   numCellFaces => mesh%NumCellFaces, &
    &     cellToCell => mesh%CellToCell, &
    &         dr_inv => lambda*SafeInverse(mesh%dr))
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    !$omp parallel do schedule(static) &
    !$omp & default(none) private(iCell, iCellFace) shared(u_bar, v)
    do iCell = 1, numCells; block
      integer :: rCell, rrCell, rrrCell, rrrrCell
      integer :: lCell, llCell, lllCell, llllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        associate(rCellFace => iCellFace, &
            &     lCellFace => Flip(iCellFace))
          rCell = cellToCell(rCellFace, iCell)
          lCell = cellToCell(lCellFace, iCell)
          if (gAccuracyOrder >= 3) then
            rrCell = cellToCell(rCellFace, rCell)
            llCell = cellToCell(lCellFace, lCell)
            if (gAccuracyOrder >= 5) then
              rrrCell = cellToCell(rCellFace, rrCell)
              lllCell = cellToCell(lCellFace, llCell)
              if (gAccuracyOrder >= 7) then
                rrrrCell = cellToCell(rCellFace, rrrCell)
                llllCell = cellToCell(lCellFace, lllCell)
              end if
            end if
          end if
        end associate
        ! ----------------------
        ! Compute FDM-approximate divergence increment.
        ! ----------------------
        select case(gAccuracyOrder)
          case(1:2)
            v(@:,iCell) = v(@:,iCell) - &
              & ( dr_inv(:,iCellFace).inner. &
              &      FD1_C2(u_bar(:,@:,lCell), &
              &             u_bar(:,@:,rCell)) )
          ! ----------------------
          case(3:4)
            v(@:,iCell) = v(@:,iCell) - &
              & ( dr_inv(:,iCellFace).inner. &
              &     FD1_C4(u_bar(:,@:,llCell), &
              &            u_bar(:,@:, lCell), &
              &            u_bar(:,@:, rCell), &
              &            u_bar(:,@:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(@:,iCell) = v(@:,iCell) - &
              & ( dr_inv(:,iCellFace).inner. &
              &    FD1_C6(u_bar(:,@:,lllCell), &
              &           u_bar(:,@:, llCell), &
              &           u_bar(:,@:,  lCell), &
              &           u_bar(:,@:,  rCell), &
              &           u_bar(:,@:, rrCell), &
              &           u_bar(:,@:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(@:,iCell) = v(@:,iCell) - &
              & ( dr_inv(:,iCellFace).inner. &
              &   FD1_C8(u_bar(:,@:,llllCell), &
              &          u_bar(:,@:, lllCell), &
              &          u_bar(:,@:,  llCell), &
              &          u_bar(:,@:,   lCell), &
              &          u_bar(:,@:,   rCell), &
              &          u_bar(:,@:,  rrCell), &
              &          u_bar(:,@:, rrrCell), &
              &          u_bar(:,@:,rrrrCell)) )
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
real(dp) elemental function FD1_F1(u, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u, u_r
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F1 = u_r - u
end function FD1_F1

!! ----------------------------------------------------------------- !!
!! Second order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F2(u, u_r, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F2 = &
    & ( (-1.5_dp)*u   + &
    &   (+2.0_dp)*u_r + &
    &   (-0.5_dp)*u_rr )
end function FD1_F2

!! ----------------------------------------------------------------- !!
!! Third order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F3(u_l, u, u_r, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l, u, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F3 = &
    & ( (-1.0_dp/3.0_dp)*u_l + &
    &   (-1.0_dp/2.0_dp)*u   + &
    &                    u_r + &
    &   (-1.0_dp/6.0_dp)*u_rr )
end function FD1_F3

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F4(u_l, u, u_r, u_rr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l, u, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F4 = &
    & ( (-01.0_dp/04.0_dp)*u_l  + &
    &   (-05.0_dp/06.0_dp)*u    + &
    &   (+03.0_dp/02.0_dp)*u_r  + &
    &   (-01.0_dp/02.0_dp)*u_rr + &
    &   (+01.0_dp/12.0_dp)*u_rrr )
end function FD1_F4

!! ----------------------------------------------------------------- !!
!! Fifth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F5(u_ll, u_l, u, u_r, u_rr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll, u_l, u, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F5 = &
    & ( (+01.0_dp/20.0_dp)*u_ll + &
    &   (-01.0_dp/02.0_dp)*u_l  + &
    &   (-01.0_dp/03.0_dp)*u    + &
    &                      u_r  + &
    &   (-01.0_dp/04.0_dp)*u_rr + &
    &   (+01.0_dp/30.0_dp)*u_rrr )
end function FD1_F5

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F6(u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F6 = &
    & ( (+01.0_dp/30.0_dp)*u_ll  + &
    &   (-02.0_dp/05.0_dp)*u_l   + &
    &   (-07.0_dp/12.0_dp)*u     + &
    &   (+04.0_dp/03.0_dp)*u_r   + &
    &   (-01.0_dp/02.0_dp)*u_rr  + &
    &   (+02.0_dp/15.0_dp)*u_rrr + &
    &   (-01.0_dp/60.0_dp)*u_rrrr )
end function FD1_F6

!! ----------------------------------------------------------------- !!
!! Seventh order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F7(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F7 = &
    & ( (-001.0_dp/105.0_dp)*u_lll + &
    &   (+001.0_dp/010.0_dp)*u_ll  + &
    &   (-003.0_dp/005.0_dp)*u_l   + &
    &   (-001.0_dp/004.0_dp)*u     + &
    &                        u_r   + &
    &   (-003.0_dp/010.0_dp)*u_rr  + &
    &   (+001.0_dp/015.0_dp)*u_rrr + &
    &   (-001.0_dp/140.0_dp)*u_rrrr )
end function FD1_F7

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD1_F8(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr, u_rrrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr, u_rrrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F8 = &
    & ( (-001.0_dp/168.0_dp)*u_lll  + &
    &   (+001.0_dp/014.0_dp)*u_ll   + &
    &   (-001.0_dp/002.0_dp)*u_l    + &
    &   (-009.0_dp/020.0_dp)*u      + &
    &   (+005.0_dp/004.0_dp)*u_r    + &
    &   (-001.0_dp/002.0_dp)*u_rr   + &
    &   (+001.0_dp/006.0_dp)*u_rrr  + &
    &   (-001.0_dp/028.0_dp)*u_rrrr + &
    &   (+001.0_dp/280.0_dp)*u_rrrrr )
end function FD1_F8

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The forward FDM-approximate gradient: v̅ ← v̅ - λ∇u.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Gradient_Forward$rank(mesh, v_bar, lambda, u, dirAll, dirFace, dirCellFace)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), pointer :: u(@:,:), v_bar(:,@:,:)
  integer(1), intent(in), optional :: dirAll, dirFace(:), dirCellFace(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  associate(numCells => mesh%NumCells, &
    &   numCellFaces => mesh%NumCellFaces, &
    &     cellToCell => mesh%CellToCell, &
    &         dr_inv => lambda*SafeInverse(mesh%dr))
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    !$omp parallel do schedule(static) &
    !$omp & default(none) private(iCell, iCellFace) &
    !$omp & shared(u, v_bar, dirAll, dirFace, dirCellFace)
    do iCell = 1, numCells; block
      integer(1) :: dir
      integer :: rCell, rrCell, rrrCell, rrrrCell, rrrrrCell
      integer :: lCell, llCell, lllCell, llllCell, lllllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Determine FD direction (default is forward).
        ! ----------------------
        dir = 1_1
        if (present(dirAll)) dir = dirAll
        if (present(dirFace)) dir = dirFace(iCellFace)
        if (present(dirCellFace)) dir = dirCellFace(iCellFace, iCell)
        ! ----------------------
        ! Find indices of the adjacent cells using the FD direction.
        ! ----------------------
        associate(inc => (1_1-dir)/2_1)
          associate(rCellFace => iCellFace+inc, &
            &       lCellFace => Flip(iCellFace+inc))
            rCell = cellToCell(rCellFace, iCell)
            lCell = cellToCell(lCellFace, iCell)
            if (gAccuracyOrder >= 2) then
              rrCell = cellToCell(rCellFace, rCell)
              llCell = cellToCell(lCellFace, lCell)
              if (gAccuracyOrder >= 4) then
                rrrCell = cellToCell(rCellFace, rrCell)
                lllCell = cellToCell(lCellFace, llCell)
                if (gAccuracyOrder >= 6) then
                  rrrrCell = cellToCell(rCellFace, rrrCell)
                  llllCell = cellToCell(lCellFace, lllCell)
                  if (gAccuracyOrder >= 8) then
                    rrrrrCell = cellToCell(rCellFace, rrrrCell)
                    lllllCell = cellToCell(lCellFace, llllCell)
                  end if
                end if
              end if
            end if
          end associate
        end associate
        ! ----------------------
        ! Compute FDM-approximate gradient increment.
        ! ----------------------
        select case(gAccuracyOrder)
          case(1)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &      dir*( dr_inv(:,iCellFace).outer. &
              &                     FD1_F1(u(@:,iCell), &
              &                            u(@:,rCell)) )
          ! ----------------------
          case(2)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &      dir*( dr_inv(:,iCellFace).outer. &
              &                    FD1_F2(u(@:, iCell), &
              &                           u(@:, rCell), &
              &                           u(@:,rrCell)) )
          ! ----------------------
          case(3)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &      dir*( dr_inv(:,iCellFace).outer. &
              &                    FD1_F3(u(@:, lCell), &
              &                           u(@:, iCell), &
              &                           u(@:, rCell), &
              &                           u(@:,rrCell)) )
          ! ----------------------
          case(4)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &      dir*( dr_inv(:,iCellFace).outer. &
              &                   FD1_F4(u(@:,  lCell), &
              &                          u(@:,  iCell), &
              &                          u(@:,  rCell), &
              &                          u(@:, rrCell), &
              &                          u(@:,rrrCell)) )
          ! ----------------------
          case(5)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &      dir*( dr_inv(:,iCellFace).outer. &
              &                   FD1_F5(u(@:, llCell), &
              &                          u(@:,  lCell), &
              &                          u(@:,  iCell), &
              &                          u(@:,  rCell), &
              &                          u(@:, rrCell), &
              &                          u(@:,rrrCell)) )
          ! ----------------------
          case(6)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &      dir*( dr_inv(:,iCellFace).outer. &
              &                  FD1_F6(u(@:,  llCell), &
              &                         u(@:,   lCell), &
              &                         u(@:,   iCell), &
              &                         u(@:,   rCell), &
              &                         u(@:,  rrCell), &
              &                         u(@:, rrrCell), &
              &                         u(@:,rrrrCell)) )
          ! ----------------------
          case(7)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &      dir*( dr_inv(:,iCellFace).outer. &
              &                  FD1_F7(u(@:, lllCell), &
              &                         u(@:,  llCell), &
              &                         u(@:,   lCell), &
              &                         u(@:,   iCell), &
              &                         u(@:,   rCell), &
              &                         u(@:,  rrCell), &
              &                         u(@:, rrrCell), &
              &                         u(@:,rrrrCell)) )
          ! ----------------------
          case(8)
            v_bar(:,@:,iCell) = v_bar(:,@:,iCell) - &
              &      dir*( dr_inv(:,iCellFace).outer. &
              &                 FD1_F8(u(@:,  lllCell), &
              &                        u(@:,   llCell), &
              &                        u(@:,    lCell), &
              &                        u(@:,    iCell), &
              &                        u(@:,    rCell), &
              &                        u(@:,   rrCell), &
              &                        u(@:,  rrrCell), &
              &                        u(@:, rrrrCell), &
              &                        u(@:,rrrrrCell)) )
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient_Forward$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The backward FDM-approximate divergence: v ← v - λ∇⋅u̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Divergence_Backward$rank(mesh, v, lambda, u_bar, dirAll, dirFace, dirCellFace)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), pointer :: u_bar(:,@:,:), v(@:,:)
  integer(1), intent(in), optional :: dirAll, dirFace(:), dirCellFace(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells => mesh%NumCells, &
    &   numCellFaces => mesh%NumCellFaces, &
    &     cellToCell => mesh%CellToCell, &
    &         dr_inv => lambda*SafeInverse(mesh%dr))
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    !$omp parallel do schedule(static) &
    !$omp & default(none) private(iCell, iCellFace) &
    !$omp & shared(u_bar, v, dirAll, dirFace, dirCellFace)
    do iCell = 1, numCells; block
      integer(1) :: dir
      integer :: rCell, rrCell, rrrCell, rrrrCell, rrrrrCell
      integer :: lCell, llCell, lllCell, llllCell, lllllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Determine FD direction (default is backward).
        ! ----------------------
        dir = -1_1
        if (present(dirAll)) dir = dirAll
        if (present(dirFace)) dir = dirFace(iCellFace)
        if (present(dirCellFace)) dir = dirCellFace(iCellFace, iCell)
        ! ----------------------
        ! Find indices of the adjacent cells using the FD direction.
        ! ----------------------
        associate(inc => (1_1-dir)/2_1)
          associate(rCellFace => iCellFace+inc, &
            &       lCellFace => Flip(iCellFace+inc))
            rCell = cellToCell(rCellFace, iCell)
            lCell = cellToCell(lCellFace, iCell)
            if (gAccuracyOrder >= 2) then
              rrCell = cellToCell(rCellFace, rCell)
              llCell = cellToCell(lCellFace, lCell)
              if (gAccuracyOrder >= 4) then
                rrrCell = cellToCell(rCellFace, rrCell)
                lllCell = cellToCell(lCellFace, llCell)
                if (gAccuracyOrder >= 6) then
                  rrrrCell = cellToCell(rCellFace, rrrCell)
                  llllCell = cellToCell(lCellFace, lllCell)
                  if (gAccuracyOrder >= 8) then
                    rrrrrCell = cellToCell(rCellFace, rrrrCell)
                    lllllCell = cellToCell(lCellFace, llllCell)
                  end if
                end if
              end if
            end if
          end associate
        end associate
        ! ----------------------
        ! Compute FDM-approximate divergence increment.
        ! ----------------------
        select case(gAccuracyOrder)
          case(1)
            v(@:,iCell) = v(@:,iCell) - &
              & dir*( dr_inv(:,iCellFace).inner. &
              &          FD1_F1(u_bar(:,@:,rCell), &
              &                 u_bar(:,@:,iCell)) )
          ! ----------------------
          case(2)
            v(@:,iCell) = v(@:,iCell) - &
              & dir*( dr_inv(:,iCellFace).inner. &
              &         FD1_F2(u_bar(:,@:, iCell), &
              &                u_bar(:,@:, rCell), &
              &                u_bar(:,@:,rrCell)) )
          ! ----------------------
          case(3)
            v(@:,iCell) = v(@:,iCell) - &
              & dir*( dr_inv(:,iCellFace).inner. &
              &         FD1_F3(u_bar(:,@:, lCell), &
              &                u_bar(:,@:, iCell), &
              &                u_bar(:,@:, rCell), &
              &                u_bar(:,@:,rrCell)) )
          ! ----------------------
          case(4)
            v(@:,iCell) = v(@:,iCell) - &
              & dir*( dr_inv(:,iCellFace).inner. &
              &         FD1_F4(u_bar(:,@:,  lCell), &
              &                u_bar(:,@:,  iCell), &
              &                u_bar(:,@:,  rCell), &
              &                u_bar(:,@:, rrCell), &
              &                u_bar(:,@:,rrrCell)) )
          ! ----------------------
          case(5)
            v(@:,iCell) = v(@:,iCell) - &
              & dir*( dr_inv(:,iCellFace).inner. &
              &         FD1_F5(u_bar(:,@:, llCell), &
              &                u_bar(:,@:,  lCell), &
              &                u_bar(:,@:,  iCell), &
              &                u_bar(:,@:,  rCell), &
              &                u_bar(:,@:, rrCell), &
              &                u_bar(:,@:,rrrCell)) )
          ! ----------------------
          case(6)
            v(@:,iCell) = v(@:,iCell) - &
              & dir*( dr_inv(:,iCellFace).inner. &
              &       FD1_F6(u_bar(:,@:,  llCell), &
              &              u_bar(:,@:,   lCell), &
              &              u_bar(:,@:,   iCell), &
              &              u_bar(:,@:,   rCell), &
              &              u_bar(:,@:,  rrCell), &
              &              u_bar(:,@:, rrrCell), &
              &              u_bar(:,@:,rrrrCell)) )
          ! ----------------------
          case(7)
            v(@:,iCell) = v(@:,iCell) - &
              & dir*( dr_inv(:,iCellFace).inner. &
              &       FD1_F7(u_bar(:,@:, lllCell), &
              &              u_bar(:,@:,  llCell), &
              &              u_bar(:,@:,   lCell), &
              &              u_bar(:,@:,   iCell), &
              &              u_bar(:,@:,   rCell), &
              &              u_bar(:,@:,  rrCell), &
              &              u_bar(:,@:, rrrCell), &
              &              u_bar(:,@:,rrrrCell)) )
          ! ----------------------
          case(8)
            v(@:,iCell) = v(@:,iCell) - &
              & dir*( dr_inv(:,iCellFace).inner. &
              &       FD1_F8(u_bar(:,@:,  lllCell), &
              &              u_bar(:,@:,   llCell), &
              &              u_bar(:,@:,    lCell), &
              &              u_bar(:,@:,    iCell), &
              &              u_bar(:,@:,    rCell), &
              &              u_bar(:,@:,   rrCell), &
              &              u_bar(:,@:,  rrrCell), &
              &              u_bar(:,@:, rrrrCell), &
              &              u_bar(:,@:,rrrrrCell)) )
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
!! The central FDM-approximate convection: v ← v - λ∇⋅uw̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine FDM_Convection_Central$rank(mesh, v, lambda, u, w_bar)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda 
  real(dp), intent(in), pointer :: u(@:,:), v(@:,:), w_bar(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable, target :: f_bar(:,@:,:)
  allocate(f_bar(size(w_bar, dim=1), @{size(u, dim=$$)}@, size(w_bar, dim=2)))
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  ! ----------------------
  ! f̅ ← uw̅,
  ! v ← v - λ∇⋅f̅.
  ! ----------------------
  call Mul_Outer(mesh, f_bar, w_bar, u)
  call FDM_Divergence_Central(mesh, v, lambda, f_bar)
end subroutine FDM_Convection_Central$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD2_C2(u_l, u, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l, u, u_r
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C2 = u_r - 2.0_dp*u + u_l
end function FD2_C2

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD2_C4(u_ll, u_l, u, u_r, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll, u_l, u, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C4 = &
    & ( (-1.0_dp/12.0_dp)*u_rr + &
    &   (+4.0_dp/03.0_dp)*u_r  + &
    &   (-5.0_dp/02.0_dp)*u    + &
    &   (+4.0_dp/03.0_dp)*u_l  + &
    &   (-1.0_dp/12.0_dp)*u_ll )
end function FD2_C4

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD2_C6(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C6 = &
    & ( (+01.0_dp/90.0_dp)*u_rrr + &
    &   (-03.0_dp/20.0_dp)*u_rr  + &
    &   (+03.0_dp/02.0_dp)*u_r   + &
    &   (-49.0_dp/18.0_dp)*u     + &
    &   (+03.0_dp/02.0_dp)*u_l   + &
    &   (-03.0_dp/20.0_dp)*u_ll  + &
    &   (+01.0_dp/90.0_dp)*u_lll )
end function FD2_C6

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
real(dp) elemental function FD2_C8(u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C8 = &
    & ( (-001.0_dp/560.0_dp)*u_rrrr + &
    &   (+008.0_dp/315.0_dp)*u_rrr  + &
    &   (-001.0_dp/005.0_dp)*u_rr   + &
    &   (+008.0_dp/005.0_dp)*u_r    + &
    &   (-205.0_dp/072.0_dp)*u      + &
    &   (+008.0_dp/005.0_dp)*u_l    + &
    &   (-001.0_dp/005.0_dp)*u_ll   + &
    &   (+008.0_dp/315.0_dp)*u_lll  + &
    &   (-001.0_dp/560.0_dp)*u_llll )
end function FD2_C8

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate Laplacian: v ← v + λΔu.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FDM_Laplacian_Central$rank(mesh, v, lambda, u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), pointer :: u(@:,:), v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells => mesh%NumCells, &
    &   numCellFaces => mesh%NumCellFaces, &
    &     cellToCell => mesh%CellToCell, &
    &     dl_sqr_inv => lambda/(mesh%dl**2))
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    !$omp parallel do schedule(static) &
    !$omp & default(none) private(iCell, iCellFace) shared(u, v)
    do iCell = 1, numCells; block
      integer :: rCell, rrCell, rrrCell, rrrrCell
      integer :: lCell, llCell, lllCell, llllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        associate(rCellFace => iCellFace, &
          &       lCellFace => Flip(iCellFace))
          rCell = cellToCell(rCellFace, iCell)
          lCell = cellToCell(lCellFace, iCell)
          if (gAccuracyOrder >= 3) then
            rrCell = cellToCell(rCellFace, rCell)
            llCell = cellToCell(lCellFace, lCell)
            if (gAccuracyOrder >= 5) then
              rrrCell = cellToCell(rCellFace, rrCell)
              lllCell = cellToCell(lCellFace, llCell)
              if (gAccuracyOrder >= 7) then
                rrrrCell = cellToCell(rCellFace, rrrCell)
                llllCell = cellToCell(lCellFace, lllCell)
              end if
            end if
          end if
        end associate
        ! ----------------------
        ! Compute FDM-approximate Laplacian increment.
        ! ----------------------
        select case(gAccuracyOrder)
          case(1:2)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dl_sqr_inv(iCellFace) * &
              &         FD2_C2(u(@:,lCell), &
              &                u(@:,iCell), &
              &                u(@:,rCell)) )
          ! ----------------------
          case(3:4)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dl_sqr_inv(iCellFace) * &
              &        FD2_C4(u(@:,llCell), &
              &               u(@:, lCell), &
              &               u(@:, iCell), &
              &               u(@:, rCell), &
              &               u(@:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dl_sqr_inv(iCellFace) * &
              &       FD2_C6(u(@:,lllCell), &
              &              u(@:, llCell), &
              &              u(@:,  lCell), &
              &              u(@:,  iCell), &
              &              u(@:,  rCell), &
              &              u(@:, rrCell), &
              &              u(@:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dl_sqr_inv(iCellFace) * &
              &      FD2_C8(u(@:,llllCell), &
              &             u(@:, lllCell), &
              &             u(@:,  llCell), &
              &             u(@:,   lCell), &
              &             u(@:,   iCell), &
              &             u(@:,   rCell), &
              &             u(@:,  rrCell), &
              &             u(@:, rrrCell), &
              &             u(@:,rrrrCell)) )
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Laplacian_Central$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate nonlinear Laplacian: v ← v + λΔf(u).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FDM_LaplacianF_Central$rank(mesh, v, lambda, f, u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), pointer :: u(@:,:), v(@:,:)
  procedure(tMapFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable, target :: w(@:,:)
  allocate(w, mold=u)
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  ! ----------------------
  ! w ← f(u),
  ! v ← v + λΔw.
  ! ----------------------
  call FuncProd(mesh, w, u, f)
  call FDM_Laplacian_Central(mesh, v, lambda, w)
end subroutine FDM_LaplacianF_Central$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate Bilaplacian: v ← v + λΔ²u.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FDM_Bilaplacian_Central$rank(mesh, v, lambda, u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), pointer :: u(@:,:), v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable, target :: w(@:,:)
  allocate(w, mold=u)
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  ! ----------------------
  ! w ← 0,
  ! w ← w + Δw.
  ! v ← v + λΔw.
  ! ----------------------
  call Fill(mesh, w, 0.0_dp)
  call FDM_Laplacian_Central(mesh, w, 1.0_dp, u)
  call FDM_Laplacian_Central(mesh, v, lambda, w)
end subroutine FDM_Bilaplacian_Central$rank
#$end do

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided d(w⋅du/dx)/dx approximation.
!! ----------------------------------------------------------------- !!
real(dp) elemental function WFD2_C2(w_l, u_l, w, u, w_r, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_l, u, u_r
  real(dp), intent(in) :: w_l, w, w_r
  ! >>>>>>>>>>>>>>>>>>>>>>
  WFD2_C2 = 0.5_dp*( (w_r+w)*(u_r-u) - (w+w_l)*(u-u_l) )
end function WFD2_C2

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided d(w⋅du/dx)/dx approximation.
!! ----------------------------------------------------------------- !!
real(dp) elemental function WFD2_C4(w_ll, u_ll, w_l, u_l, w, u, w_r, u_r, w_rr, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_ll, u_l, u, u_r, u_rr
  real(dp), intent(in) :: w_ll, w_l, w, w_r, w_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  WFD2_C4 = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO
end function WFD2_C4

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided d(w⋅du/dx)/dx approximation.
!! ----------------------------------------------------------------- !!
real(dp) elemental function WFD2_C6(w_lll, u_lll, w_ll, u_ll, w_l, u_l, &
  & w, u, w_r, u_r, w_rr, u_rr, w_rrr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr
  real(dp), intent(in) :: w_lll, w_ll, w_l, w, w_r, w_rr, w_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  WFD2_C6 = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO
end function WFD2_C6

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided d(w⋅du/dx)/dx approximation.
!! ----------------------------------------------------------------- !!
real(dp) elemental function WFD2_C8(w_llll, u_llll, w_lll, u_lll, w_ll, u_ll, &
  & w_l, u_l, w, u, w_r, u_r, w_rr, u_rr, w_rrr, u_rrr, w_rrrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  real(dp), intent(in) :: w_llll, w_lll, w_ll, w_l, w, w_r, w_rr, w_rrr, w_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  WFD2_C8 = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO
end function WFD2_C8

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate weighted Laplacian: v ← v + λ∇⋅(w∇u).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FDM_DivWGrad_Central$rank(mesh, v, lambda, w, u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), pointer :: u(@:,:), v(@:,:), w(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells => mesh%NumCells, &
    &   numCellFaces => mesh%NumCellFaces, &
    &     cellToCell => mesh%CellToCell, &
    &     dl_sqr_inv => lambda/(mesh%dl**2))
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    !$omp parallel do schedule(static) &
    !$omp & default(none) private(iCell, iCellFace) shared(u, v, w)
    do iCell = 1, numCells; block
      integer :: rCell, rrCell, rrrCell, rrrrCell
      integer :: lCell, llCell, lllCell, llllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        ! ----------------------
        associate(rCellFace => iCellFace, &
          &       lCellFace => Flip(iCellFace))
          rCell = cellToCell(rCellFace, iCell)
          lCell = cellToCell(lCellFace, iCell)
          if (gAccuracyOrder >= 3) then
            rrCell = cellToCell(rCellFace, rCell)
            llCell = cellToCell(lCellFace, lCell)
            if (gAccuracyOrder >= 5) then
              rrrCell = cellToCell(rCellFace, rrCell)
              lllCell = cellToCell(lCellFace, llCell)
              if (gAccuracyOrder >= 7) then
                rrrrCell = cellToCell(rCellFace, rrrCell)
                llllCell = cellToCell(lCellFace, lllCell)
              end if
            end if
          end if
        end associate
        ! ----------------------
        ! Compute FDM-approximate weighted Laplacian increment.
        ! ----------------------
        select case(gAccuracyOrder)
          case(1:2)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dl_sqr_inv(iCellFace) * &
              &   WFD2_C2(w(@:,lCell), u(@:,lCell), &
              &           w(@:,lCell), u(@:,iCell), &
              &           w(@:,lCell), u(@:,rCell)) )
          ! ----------------------
          case(3:4)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dl_sqr_inv(iCellFace) * &
              &   WFD2_C4(w(@:,llCell), u(@:,llCell), &
              &           w(@:, lCell), u(@:, lCell), &
              &           w(@:, iCell), u(@:, iCell), &
              &           w(@:, rCell), u(@:, rCell), &
              &           w(@:,rrCell), u(@:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dl_sqr_inv(iCellFace) * &
              &   WFD2_C6(w(@:,lllCell), u(@:,lllCell), &
              &           w(@:, llCell), u(@:, llCell), &
              &           w(@:,  lCell), u(@:,  lCell), &
              &           w(@:,  iCell), u(@:,  iCell), &
              &           w(@:,  rCell), u(@:,  rCell), &
              &           w(@:, rrCell), u(@:, rrCell), &
              &           w(@:,rrrCell), u(@:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(@:,iCell) = v(@:,iCell) + &
              & ( dl_sqr_inv(iCellFace) * &
              &   WFD2_C8(w(@:,llllCell), u(@:,llllCell), &
              &           w(@:, lllCell), u(@:, lllCell), &
              &           w(@:,  llCell), u(@:,  llCell), &
              &           w(@:,   lCell), u(@:,   lCell), &
              &           w(@:,   iCell), u(@:,   iCell), &
              &           w(@:,   rCell), u(@:,   rCell), &
              &           w(@:,  rrCell), u(@:,  rrCell), &
              &           w(@:, rrrCell), u(@:, rrrCell), &
              &           w(@:,rrrrCell), u(@:,rrrrCell)) )
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_DivWGrad_Central$rank
#$end do

end module StormRuler_FDM_Operators
