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

implicit none

integer, parameter :: ACCURACY_ORDER = 8

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

interface FDM_Convection_Central
#@do rank = 0, NUM_RANKS-1
  module procedure FDM_Convection_Central$rank
#@end do
end interface FDM_Convection_Central

private :: FD2_C2,FD2_C4,FD2_C6,FD2_C8

interface FDM_Laplacian
#@do rank = 0, NUM_RANKS
  module procedure FDM_Laplacian$rank
#@end do
end interface FDM_Laplacian

interface FDM_LaplacianF
#@do rank = 0, NUM_RANKS
  module procedure FDM_LaplacianF$rank
#@end do
end interface FDM_LaplacianF

interface FDM_Bilaplacian
#@do rank = 0, NUM_RANKS
  module procedure FDM_Bilaplacian$rank
#@end do
end interface FDM_Bilaplacian

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
contains
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

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
!! The central FDM-approximate gradient: v̅ ← v̅ - λ∇u.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Gradient_Central$rank(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda,u(@:,:)
  real(dp), intent(inout) :: v(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell, &
      &        invDn=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(iCell,iCellFace)
        lCell = cellToCell(iCell,iCellFace+1)
        if (ACCURACY_ORDER >= 4) then
          rrCell = cellToCell(rCell,iCellFace)
          llCell = cellToCell(lCell,iCellFace+1)
          if (ACCURACY_ORDER >= 6) then
            rrrCell = cellToCell(rrCell,iCellFace)
            lllCell = cellToCell(llCell,iCellFace+1)
            if (ACCURACY_ORDER >= 8) then
              rrrrCell = cellToCell(rrrCell,iCellFace)
              llllCell = cellToCell(lllCell,iCellFace+1)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate gradient increment.
        select case (ACCURACY_ORDER)
          case (1:2)
            v(:,@:,iCell) -= &
              & Outer(invDn(:,iCellFace), &
              &       FD1_C2(u(@:,rCell),u(@:,lCell)))
          case (3:4)
            v(:,@:,iCell) -= &
              & Outer(invDn(:,iCellFace), &
              &       FD1_C4(u(@:,rCell),u(@:,rrCell), &
              &              u(@:,lCell),u(@:,llCell)))
          case (5:6)
            v(:,@:,iCell) -= &
              & Outer(invDn(:,iCellFace), &
              &       FD1_C6(u(@:,rCell),u(@:,rrCell),u(@:,rrrCell), &
              &              u(@:,lCell),u(@:,llCell),u(@:,lllCell)))
          case (7:8)
            v(:,@:,iCell) -= &
              & Outer(invDn(:,iCellFace), &
              &       FD1_C8(u(@:,rCell),u(@:,rrCell),u(@:,rrrCell),u(@:,rrrrCell), &
              &              u(@:,lCell),u(@:,llCell),u(@:,lllCell),u(@:,llllCell)))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient_Central$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The central FDM-approximate divergence: v ← v - λ∇⋅u̅.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Divergence_Central$rank(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda,u(:,@:,:)
  real(dp), intent(inout) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell, &
      &        invDn=>lambda*SafeInverse(mesh%Dn))
    ! ----------------------
    !$omp parallel do
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(iCell,iCellFace)
        lCell = cellToCell(iCell,iCellFace+1)
        if (ACCURACY_ORDER >= 4) then
          rrCell = cellToCell(rCell,iCellFace)
          llCell = cellToCell(lCell,iCellFace+1)
          if (ACCURACY_ORDER >= 6) then
            rrrCell = cellToCell(rrCell,iCellFace)
            lllCell = cellToCell(llCell,iCellFace+1)
            if (ACCURACY_ORDER >= 8) then
              rrrrCell = cellToCell(rrrCell,iCellFace)
              llllCell = cellToCell(lllCell,iCellFace+1)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate divergence increment.
        select case (ACCURACY_ORDER)
          case (1:2)
            v(@:,iCell) -= &
              & Inner(invDn(:,iCellFace), &
              &       FD1_C2(u(:,@:,rCell),u(:,@:,lCell)))
          case (3:4)
            v(@:,iCell) -= &
              & Inner(invDn(:,iCellFace), &
              &       FD1_C4(u(:,@:,rCell),u(:,@:,rrCell), &
              &              u(:,@:,lCell),u(:,@:,llCell)))
          case (5:6)
            v(@:,iCell) -= &
              & Inner(invDn(:,iCellFace), &
              &       FD1_C6(u(:,@:,rCell),u(:,@:,rrCell),u(:,@:,rrrCell), &
              &              u(:,@:,lCell),u(:,@:,llCell),u(:,@:,lllCell)))
          case (7:8)
            v(@:,iCell) -= &
              & Inner(invDn(:,iCellFace), &
              &       FD1_C8(u(:,@:,rCell),u(:,@:,rrCell),u(:,@:,rrrCell),u(:,@:,rrrrCell), &
              &              u(:,@:,lCell),u(:,@:,llCell),u(:,@:,lllCell),u(:,@:,llllCell)))
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

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! The central FDM-approximate convection: v ← v - λ∇⋅uw̅.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Convection_Central$rank(mesh,v,lambda,u,w)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda,u(@:,:),w(:,:)
  real(dp), intent(inout) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: f(:,@:,:)
  allocate(f(size(w,dim=1), &
    & @{size(u,dim=$$+1)}@,size(w,dim=2)))
  ! ----------------------
  ! Fast exit in case λ=0.
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  ! f̅ ← uw̅,
  ! v ← v - λ∇⋅f̅.
  call Mul_Outer(mesh,f,w,u)
  call FDM_Divergence_Central(mesh,v,lambda,f)
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
!! The FDM-approximate Laplacian: v ← v + λΔu.
#@do rank = 0, NUM_RANKS
subroutine FDM_Laplacian$rank(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda,u(@:,:)
  real(dp), intent(inout) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, iCellFace
  ! ----------------------
  ! Fast exit in case λ=0.
  if (lambda==0.0_dp) then
    return
  end if
  ! ----------------------
  associate(numCells=>mesh%NumCells, &
      & numCellFaces=>mesh%NumCellFaces, &
      &   cellToCell=>mesh%CellToCell, &
      &     invDxSqr=>lambda/mesh%Dx**2)
    ! ----------------------
    !$omp parallel do
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do iCellFace = 1, numCellFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(iCell,iCellFace)
        lCell = cellToCell(iCell,iCellFace+1)
        if (ACCURACY_ORDER >= 4) then
          rrCell = cellToCell(rCell,iCellFace)
          llCell = cellToCell(lCell,iCellFace+1)
          if (ACCURACY_ORDER >= 6) then
            rrrCell = cellToCell(rrCell,iCellFace)
            lllCell = cellToCell(llCell,iCellFace+1)
            if (ACCURACY_ORDER >= 8) then
              rrrrCell = cellToCell(rrrCell,iCellFace)
              llllCell = cellToCell(lllCell,iCellFace+1)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate Laplacian increment.
        select case (ACCURACY_ORDER)
          case (1:2)
            v(@:,iCell) += &
              & invDxSqr(iCellFace)* &
              & FD2_C2(u(@:,iCell), &
              &        u(@:,rCell),u(@:,lCell))
          case (3:4)
            v(@:,iCell) += &
              & invDxSqr(iCellFace)* &
              & FD2_C4(u(@:,iCell), &
              &        u(@:,rCell),u(@:,rrCell), &
              &        u(@:,lCell),u(@:,llCell))
          case (5:6)
            v(@:,iCell) += &
              & invDxSqr(iCellFace)* &
              & FD2_C6(u(@:,iCell), &
              &        u(@:,rCell),u(@:,rrCell),u(@:,rrrCell), &
              &        u(@:,lCell),u(@:,llCell),u(@:,lllCell))
          case (7:8)
            v(@:,iCell) += &
              & invDxSqr(iCellFace)* &
              & FD2_C8(u(@:,iCell), &
              &        u(@:,rCell),u(@:,rrCell),u(@:,rrrCell),u(@:,rrrrCell), &
              &        u(@:,lCell),u(@:,llCell),u(@:,lllCell),u(@:,llllCell))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Laplacian$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate nonlinear Laplacian: v ← v + λΔf(u).
#@do rank = 0, NUM_RANKS
subroutine FDM_LaplacianF$rank(mesh,v,lambda,f,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda,u(@:,:)
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
  call FDM_Laplacian(mesh,v,lambda,w)
end subroutine FDM_LaplacianF$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate Bilaplacian: v ← v + λΔ²u.
#@do rank = 0, NUM_RANKS
subroutine FDM_Bilaplacian$rank(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda,u(@:,:)
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
  call FDM_Laplacian(mesh,w,1.0_dp,u)
  call FDM_Laplacian(mesh,v,lambda,w)
end subroutine FDM_Bilaplacian$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

end module StormRuler_FiniteDifferences
