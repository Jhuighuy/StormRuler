!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! Copyright (C) 2021 Oleg Butakov
!! 
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights  to use, 
!! copy, modify, merge, publish, distribute, sublicense, and/or
!! sell copies of the Software, and to permit persons to whom the  
!! Software is fUrnished to do so, subject to the following 
!! conditions:
!! 
!! The above copyright notice and this permission notice shall be 
!! included in all copies or substantial portions of the Software.
!! 
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
!! OF MERCHANTABILITY, FITNESS FOR A PARTICUlAR PUrPOSE AND 
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module StormRuler_FiniteDifferences

use StormRuler_Helpers
use StormRuler_Arithmetics
use StormRuler_Mesh
#@use 'StormRuler_Parameters.f90'

implicit none

integer, parameter :: ACCURACY_ORDER = 8

private :: FD1_C2, FD1_C4, FD1_C6, FD1_C8

interface FDM_Gradient
#@do rank = 0, NUM_RANKS-1
  module procedure FDM_Gradient$rank
#@end do
end interface FDM_Gradient

interface FDM_Divergence
#@do rank = 0, NUM_RANKS-1
  module procedure FDM_Divergence$rank
#@end do
end interface FDM_Divergence

interface FDM_Convection
#@do rank = 0, 0*(NUM_RANKS-1)
  module procedure FDM_Convection$rank
#@end do
end interface FDM_Convection

private :: FD2_C2, FD2_C4, FD2_C6, FD2_C8

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

contains

!! -----------------------------------------------------------------  
!! Second order accracy central undivided finite difference.
elemental function FD1_C2(Ur, &
                          Ul) result(FD)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: Ur
  real(dp), intent(in) :: Ul
  real(dp) :: FD
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD = (+1.0_dp/2.0_dp)*Ur &
     + (-1.0_dp/2.0_dp)*Ul
end function FD1_C2
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! Fourth order accracy central undivided finite difference.
elemental function FD1_C4(Ur,Urr, &
                          Ul,Ull) result(FD)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: Ur,Urr
  real(dp), intent(in) :: Ul,Ull
  real(dp) :: FD
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD = (-1.0_dp/12.0_dp)*Urr &
     + (+2.0_dp/03.0_dp)*Ur  &
     + (-2.0_dp/03.0_dp)*Ul  &
     + (+1.0_dp/12.0_dp)*Ull
end function FD1_C4
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! Sixth order accracy central undivided finite difference.
elemental function FD1_C6(Ur,Urr,Urrr, &
                          Ul,Ull,Ulll) result(FD)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: Ur,Urr,Urrr
  real(dp), intent(in) :: Ul,Ull,Ulll
  real(dp) :: FD
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD = (+01.0_dp/60.0_dp)*Urrr &
     + (-03.0_dp/20.0_dp)*Urr  &
     + (+03.0_dp/04.0_dp)*Ur   &
     + (-03.0_dp/04.0_dp)*Ul   &
     + (+03.0_dp/20.0_dp)*Ull  &
     + (-01.0_dp/60.0_dp)*Ulll
end function FD1_C6
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! Eighth order accracy central undivided finite difference.
elemental function FD1_C8(Ur,Urr,Urrr,Urrrr, &
                          Ul,Ull,Ulll,Ullll) result(FD)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: Ur,Urr,Urrr,Urrrr
  real(dp), intent(in) :: Ul,Ull,Ulll,Ullll
  real(dp) :: FD
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD = (-001.0_dp/280.0_dp)*Urrrr &
     + (+004.0_dp/105.0_dp)*Urrr  &
     + (-001.0_dp/005.0_dp)*Urr   &
     + (+004.0_dp/005.0_dp)*Ur    &
     + (-004.0_dp/005.0_dp)*Ul    &
     + (+001.0_dp/005.0_dp)*Ull   &
     + (-004.0_dp/105.0_dp)*Ulll  &
     + (+001.0_dp/280.0_dp)*Ullll
end function FD1_C8
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! The central FDM-approximate gradient: V ← V - λ∇U.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Gradient$rank(mesh,V,lambda,U)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: V(:,@:,:),U(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, jFace
  associate(numFaces=>mesh%NumFaces, &
            numCells=>mesh%NumCells, &
          cellToCell=>mesh%CellToCell, &
               invDn=>lambda*SafeInverse(mesh%Dn))
    !$omp parallel do
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do jFace = 1, numFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(iCell,jFace)
        lCell = cellToCell(iCell,jFace+1)
        if (ACCURACY_ORDER >= 4) then
          rrCell = cellToCell(rCell,jFace)
          llCell = cellToCell(lCell,jFace+1)
          if (ACCURACY_ORDER >= 6) then
            rrrCell = cellToCell(rrCell,jFace)
            lllCell = cellToCell(llCell,jFace+1)
            if (ACCURACY_ORDER >= 8) then
              rrrrCell = cellToCell(rrrCell,jFace)
              llllCell = cellToCell(lllCell,jFace+1)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate gradient increment.
        select case (ACCURACY_ORDER)
        case (1:2)
          V(:,@:,iCell) -= &
            Outer(invDn(:,jFace), &
                  FD1_C2(U(@:,rCell),U(@:,lCell)))
        case (3:4)
          V(:,@:,iCell) -= &
            Outer(invDn(:,jFace), &
                  FD1_C4(U(@:,rCell),U(@:,rrCell), &
                         U(@:,lCell),U(@:,llCell)))
        case (5:6)
          V(:,@:,iCell) -= &
            Outer(invDn(:,jFace), &
                  FD1_C6(U(@:,rCell),U(@:,rrCell),U(@:,rrrCell), &
                         U(@:,lCell),U(@:,llCell),U(@:,lllCell)))
        case (7:8)
          V(:,@:,iCell) -= &
            Outer(invDn(:,jFace), &
                  FD1_C8(U(@:,rCell),U(@:,rrCell),U(@:,rrrCell),U(@:,rrrrCell), &
                         U(@:,lCell),U(@:,llCell),U(@:,lllCell),U(@:,llllCell)))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient$rank
#@end do
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! The central FDM-approximate divergence: V ← V - λ∇⋅U.
#@do rank = 0, NUM_RANKS-1
subroutine FDM_Divergence$rank(mesh,V,lambda,U)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: V(@:,:),U(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, jFace
  associate(numFaces=>mesh%NumFaces, &
            numCells=>mesh%NumCells, &
          cellToCell=>mesh%CellToCell, &
               invDn=>lambda*SafeInverse(mesh%Dn))
    !$omp parallel do
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do jFace = 1, numFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(iCell,jFace)
        lCell = cellToCell(iCell,jFace+1)
        if (ACCURACY_ORDER >= 4) then
          rrCell = cellToCell(rCell,jFace)
          llCell = cellToCell(lCell,jFace+1)
          if (ACCURACY_ORDER >= 6) then
            rrrCell = cellToCell(rrCell,jFace)
            lllCell = cellToCell(llCell,jFace+1)
            if (ACCURACY_ORDER >= 8) then
              rrrrCell = cellToCell(rrrCell,jFace)
              llllCell = cellToCell(lllCell,jFace+1)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate divergence increment.
        select case (ACCURACY_ORDER)
        case (1:2)
          V(@:,iCell) -= &
            Inner(invDn(:,jFace), &
                  FD1_C2(U(:,@:,rCell),U(:,@:,lCell)))
        case (3:4)
          V(@:,iCell) -= &
            Inner(invDn(:,jFace), &
                  FD1_C4(U(:,@:,rCell),U(:,@:,rrCell), & 
                         U(:,@:,lCell),U(:,@:,llCell)))
        case (5:6)
          V(@:,iCell) -= &
            Inner(invDn(:,jFace), &
                  FD1_C6(U(:,@:,rCell),U(:,@:,rrCell),U(:,@:,rrrCell), &
                         U(:,@:,lCell),U(:,@:,llCell),U(:,@:,lllCell)))
        case (7:8)
          V(@:,iCell) -= &
            Inner(invDn(:,jFace), &
                  FD1_C8(U(:,@:,rCell),U(:,@:,rrCell),U(:,@:,rrrCell),U(:,@:,rrrrCell), &
                         U(:,@:,lCell),U(:,@:,llCell),U(:,@:,lllCell),U(:,@:,llllCell)))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Divergence$rank
#@end do
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! The central FDM-approximate convection: V ← V - λ∇⋅CU.
#@do rank = 0, 0*(NUM_RANKS-1)
subroutine FDM_Convection$rank(mesh,V,lambda,c,U)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: V(@:,:),C(@:,:),U(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: W(:,@:,:)
  allocate(W, mold=U)
  ! ----------------------
  ! W ← CU,
  ! V ← V - λ∇⋅U.
  call Mul(mesh,W,C,U)
  call FDM_Divergence(mesh,V,lambda,W)
end subroutine FDM_Convection$rank
#@end do
!! -----------------------------------------------------------------  

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
!! Second order accracy central undivided second finite difference.
elemental function FD2_C2(U, &
                          Ur, &
                          Ul) result(FD)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: U
  real(dp), intent(in) :: Ur
  real(dp), intent(in) :: Ul
  real(dp) :: FD
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD = (+1.0_dp)*Ur &
     + (-2.0_dp)*U  &
     + (+1.0_dp)*Ul
end function FD2_C2
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! Fourth order accracy central undivided second finite difference.
elemental function FD2_C4(U, &
                          Ur,Urr, &
                          Ul,Ull) result(FD)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: U
  real(dp), intent(in) :: Ur,Urr
  real(dp), intent(in) :: Ul,Ull
  real(dp) :: FD
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD = (-1.0_dp/12.0_dp)*Urr &
     + (+4.0_dp/03.0_dp)*Ur  &
     + (-5.0_dp/02.0_dp)*U   &
     + (+4.0_dp/03.0_dp)*Ul  &
     + (-1.0_dp/12.0_dp)*Ull
end function FD2_C4
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! Sixth order accracy central undivided second finite difference.
elemental function FD2_C6(U, &
                          Ur,Urr,Urrr, &
                          Ul,Ull,Ulll) result(FD)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: U
  real(dp), intent(in) :: Ur,Urr,Urrr
  real(dp), intent(in) :: Ul,Ull,Ulll
  real(dp) :: FD
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD = (+01.0_dp/90.0_dp)*Urrr &
     + (-03.0_dp/20.0_dp)*Urr  &
     + (+03.0_dp/02.0_dp)*Ur   &
     + (-49.0_dp/18.0_dp)*U    &
     + (+03.0_dp/02.0_dp)*Ul   &
     + (-03.0_dp/20.0_dp)*Ull  &
     + (+01.0_dp/90.0_dp)*Ulll
end function FD2_C6
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! Eighth order accracy central undivided second finite difference.
elemental function FD2_C8(U, &
                          Ur,Urr,Urrr,Urrrr, &
                          Ul,Ull,Ulll,Ullll) result(FD)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: U
  real(dp), intent(in) :: Ur,Urr,Urrr,Urrrr
  real(dp), intent(in) :: Ul,Ull,Ulll,Ullll
  real(dp) :: FD
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD = (-001.0_dp/560.0_dp)*Urrrr &
     + (+008.0_dp/315.0_dp)*Urrr  &
     + (-001.0_dp/005.0_dp)*Urr   &
     + (+008.0_dp/005.0_dp)*Ur    &
     + (-205.0_dp/072.0_dp)*U     &
     + (+008.0_dp/005.0_dp)*Ul    &
     + (-001.0_dp/005.0_dp)*Ull   &
     + (+008.0_dp/315.0_dp)*Ulll  &
     + (-001.0_dp/560.0_dp)*Ullll
end function FD2_C8
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! The FDM-approximate Laplacian: V ← V + λΔu.
#@do rank = 0, NUM_RANKS
subroutine FDM_Laplacian$rank(mesh,V,lambda,U)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: V(@:,:),U(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell, jFace
  associate(numFaces=>mesh%NumFaces, &
            numCells=>mesh%NumCells, &
          cellToCell=>mesh%CellToCell, &
            invDxSqr=>lambda/mesh%Dx**2)
    !$omp parallel do
    do iCell = 1, numCells; block
      integer :: rCell,rrCell,rrrCell,rrrrCell
      integer :: lCell,llCell,lllCell,llllCell
      do jFace = 1, numFaces, 2
        ! ----------------------
        ! Find indices of the adjacent cells.
        rCell = cellToCell(iCell,jFace)
        lCell = cellToCell(iCell,jFace+1)
        if (ACCURACY_ORDER >= 4) then
          rrCell = cellToCell(rCell,jFace)
          llCell = cellToCell(lCell,jFace+1)
          if (ACCURACY_ORDER >= 6) then
            rrrCell = cellToCell(rrCell,jFace)
            lllCell = cellToCell(llCell,jFace+1)
            if (ACCURACY_ORDER >= 8) then
              rrrrCell = cellToCell(rrrCell,jFace)
              llllCell = cellToCell(lllCell,jFace+1)
            end if
          end if
        end if
        ! ----------------------
        ! Compute FDM-approximate Laplacian increment.
        select case (ACCURACY_ORDER)
        case (1:2)
          V(@:,iCell) += &
            invDxSqr(jFace) &
            *FD2_C2(U(@:,iCell), &
                    U(@:,rCell),U(@:,lCell))
        case (3:4)
          V(@:,iCell) += &
            invDxSqr(jFace) &
            *FD2_C4(U(@:,iCell), &
                    U(@:,rCell),U(@:,rrCell), &
                    U(@:,lCell),U(@:,llCell))
        case (5:6)
          V(@:,iCell) += &
            invDxSqr(jFace) &
            *FD2_C4(U(@:,iCell), &
                    U(@:,rCell),U(@:,rrCell), &
                    U(@:,lCell),U(@:,llCell))
        case (7:8)
          V(@:,iCell) += &
            invDxSqr(jFace) &
            *FD2_C8(U(@:,iCell), &
                    U(@:,rCell),U(@:,rrCell),U(@:,rrrCell),U(@:,rrrrCell), &
                    U(@:,lCell),U(@:,llCell),U(@:,lllCell),U(@:,llllCell))
        end select
      end do
    end block; end do
    !$omp end parallel do
  end associate
end subroutine FDM_Laplacian$rank
#@end do
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! The FDM-approximate nonlinear Laplacian: V ← V + λΔf(U).
#@do rank = 0, NUM_RANKS
subroutine FDM_LaplacianF$rank(mesh,V,lambda,f,U)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: V(@:,:),U(@:,:)
  procedure(MathFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: W(@:,:)
  allocate(W, mold=U)
  ! ----------------------
  ! W ← f(U),
  ! V ← V + λΔW.
  call ApplyFunc(mesh,W,U,f)
  call FDM_Laplacian(mesh,V,lambda,W)
end subroutine FDM_LaplacianF$rank
#@end do
!! -----------------------------------------------------------------  

!! -----------------------------------------------------------------  
!! The FDM-approximate Bilaplacian: V ← V + λΔ²U.
#@do rank = 0, NUM_RANKS
subroutine FDM_Bilaplacian$rank(mesh,V,lambda,U)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: V(@:,:),U(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: W(@:,:)
  allocate(W, mold=U)
  ! ----------------------
  ! W ← 0,
  ! W ← W + ΔU.
  ! V ← V + λΔW.
  call Zero(mesh,W)
  call FDM_Laplacian(mesh,W,1.0_dp,U)
  call FDM_Laplacian(mesh,V,lambda,W)
end subroutine FDM_Bilaplacian$rank
#@end do
!! -----------------------------------------------------------------  

end module StormRuler_FiniteDifferences
