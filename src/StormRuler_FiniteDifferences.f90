!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module StormRuler_FiniteDifferences

use StormRuler_Helpers
use StormRuler_Arithmetics
use StormRuler_Mesh

implicit none

integer, parameter :: ACCURACY_ORDER = 8

contains

!! -----------------------------------------------------------------  
!! Second order accracy centeral undivided finite difference.
elemental function FD1_C2(ur &
                         ,ul) result(d)
  real(dp), intent(in) :: ur
  real(dp), intent(in) :: ul
  real(dp) :: d
  d = &
    +1.0_dp/2.0_dp*ur &
    -1.0_dp/2.0_dp*ul
end function FD1_C2
!! Fourth order accracy centeral undivided finite difference.
elemental function FD1_C4(ur,urr &
                         ,ul,ull) result(d)
  real(dp), intent(in) :: ur,urr
  real(dp), intent(in) :: ul,ull
  real(dp) :: d
  d = -1.0_dp/12.0_dp*urr &
      +2.0_dp/03.0_dp*ur  &
      -2.0_dp/03.0_dp*ul  &
      +1.0_dp/12.0_dp*ull
end function FD1_C4
!! Sixth order accracy centeral undivided finite difference.
elemental function FD1_C6(ur,urr,urrr &
                         ,ul,ull,ulll) result(d)
  real(dp), intent(in) :: ur,urr,urrr
  real(dp), intent(in) :: ul,ull,ulll
  real(dp) :: d
  d = +01.0_dp/60.0_dp*urrr &
      -03.0_dp/20.0_dp*urr  &
      +03.0_dp/04.0_dp*ur   &
      -03.0_dp/04.0_dp*ul   &
      +03.0_dp/20.0_dp*ull  &
      -01.0_dp/60.0_dp*ulll
end function FD1_C6
!! Eighth order accracy centeral undivided finite difference.
elemental function FD1_C8(ur,urr,urrr,urrrr &
                         ,ul,ull,ulll,ullll) result(d)
  real(dp), intent(in) :: ur,urr,urrr,urrrr
  real(dp), intent(in) :: ul,ull,ulll,ullll
  real(dp) :: d
  d = -001.0_dp/280.0_dp*urrrr &
      +004.0_dp/105.0_dp*urrr  &
      -001.0_dp/005.0_dp*urr   &
      +004.0_dp/005.0_dp*ur    &
      -004.0_dp/005.0_dp*ul    &
      +001.0_dp/005.0_dp*ull   &
      -004.0_dp/105.0_dp*ulll  &
      +001.0_dp/280.0_dp*ullll
end function FD1_C8

!! -----------------------------------------------------------------  
!! The central FDM-approximate gradient: `v ← v - λ∇u`.
subroutine FDM_Gradient(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: v(:,:),u(:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: kDir
  integer :: iCell 
  integer :: rCell,rrCell,rrrCell,rrrrCell
  integer :: lCell,llCell,lllCell,llllCell 
  associate(numDir=>mesh%NumDir &
         ,numCells=>mesh%NumCells &
       ,cellToCell=>mesh%CellToCell &
            ,invDn=>lambda*SafeInverse(mesh%Dn))
    !$omp parallel do firstprivate(rCell)
    do iCell = 1, numCells
      do kDir = 1, numDir, 2
        ! ----------------------
        ! Second order accracy.
        rCell=cellToCell(iCell,kDir)
        lCell=cellToCell(iCell,kDir+1)
        if (ACCURACY_ORDER <= 2) then
          v(:,iCell) = v(:,iCell) &
            - invDn(:,kDir)*FD1_C2(u(rCell),u(lCell))
        else
          ! ----------------------
          ! Fourth order accracy.
          rrCell=cellToCell(rCell,kDir)
          llCell=cellToCell(lCell,kDir+1)
          if (ACCURACY_ORDER <= 4) then
            v(:,iCell) = v(:,iCell) &
              - invDn(:,kDir)*FD1_C4(u(rCell),u(rrCell)&
                                    ,u(lCell),u(llCell))
          else
            ! ----------------------
            ! Sixth order accracy.
            rrrCell=cellToCell(rrCell,kDir)
            lllCell=cellToCell(llCell,kDir+1)
            if (ACCURACY_ORDER <= 6) then
              v(:,iCell) = v(:,iCell) &
                - invDn(:,kDir)*FD1_C6(u(rCell),u(rrCell),u(rrrCell)&
                                      ,u(lCell),u(llCell),u(lllCell))
            else ! 6
              ! ----------------------
              ! Eighth order accracy.
              rrrrCell=cellToCell(rrrCell,kDir)
              llllCell=cellToCell(lllCell,kDir+1)
              v(:,iCell) = v(:,iCell) &
                - invDn(:,kDir)*FD1_C8(u(rCell),u(rrCell),u(rrrCell),u(rrrrCell)&
                                      ,u(lCell),u(llCell),u(lllCell),u(llllCell))
            end if ! 8 
          end if ! 4
        end if ! 2
      end do
    end do
    !$omp end parallel do
  end associate
end subroutine FDM_Gradient

!! -----------------------------------------------------------------  
!! The central FDM-approximate divergence: v ← v - λ∇⋅u.
subroutine FDM_Divergence(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: v(:),u(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: kDir
  integer :: iCell 
  integer :: rCell,rrCell,rrrCell,rrrrCell
  integer :: lCell,llCell,lllCell,llllCell 
  associate(numDir=>mesh%NumDir &
         ,numCells=>mesh%NumCells &
       ,cellToCell=>mesh%CellToCell &
            ,invDn=>lambda*SafeInverse(mesh%Dn))
    !$omp parallel do firstprivate(rCell)
    do iCell = 1, numCells
      do kDir = 1, numDir, 2
        ! ----------------------
        ! Second order accracy.
        rCell=cellToCell(iCell,kDir)
        lCell=cellToCell(iCell,kDir+1)
        if (ACCURACY_ORDER <= 2) then
          v(iCell) = v(iCell) &
            - dot_product(invDn(:,kDir),FD1_C2(u(:,rCell),u(:,lCell)))
        else
          ! ----------------------
          ! Fourth order accracy.
          rrCell=cellToCell(rCell,kDir)
          llCell=cellToCell(lCell,kDir+1)
          if (ACCURACY_ORDER <= 4) then
            v(iCell) = v(iCell) &
              - dot_product(invDn(:,kDir),FD1_C4(u(:,rCell),u(:,rrCell)& 
                                                ,u(:,lCell),u(:,llCell)))
          else
            ! ----------------------
            ! Sixth order accracy.
            rrrCell=cellToCell(rrCell,kDir)
            lllCell=cellToCell(llCell,kDir+1)
            if (ACCURACY_ORDER <= 6) then
              v(iCell) = v(iCell) &
                - dot_product(invDn(:,kDir),FD1_C6(u(:,rCell),u(:,rrCell),u(:,rrrCell)&
                                                  ,u(:,lCell),u(:,llCell),u(:,lllCell)))
            else ! 6
              ! ----------------------
              ! Eighth order accracy.
              rrrrCell=cellToCell(rrrCell,kDir)
              llllCell=cellToCell(lllCell,kDir+1)
              v(iCell) = v(iCell) &
                - dot_product(invDn(:,kDir),FD1_C8(u(:,rCell),u(:,rrCell),u(:,rrrCell),u(:,rrrrCell)&
                                                  ,u(:,lCell),u(:,llCell),u(:,lllCell),u(:,llllCell)))
            end if ! 8
          end if ! 4
        end if ! 2
      end do
    end do
    !$omp end parallel do
  end associate
end subroutine FDM_Divergence

!! -----------------------------------------------------------------  
!! The central FDM-approximate convection: v ← v - λ∇⋅cu.
subroutine FDM_Convection(mesh,v,lambda,c,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: v(:),c(:),u(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: w(:,:)
  allocate(w, mold=u)
  ! ----------------------
  ! w ← cu,
  ! v ← v - λ∇⋅w.
  call Mul(mesh,w,c,u)
  call FDM_Divergence(mesh,v,lambda,w)
end subroutine FDM_Convection

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!! -----------------------------------------------------------------  
!! Second order accracy centeral undivided second finite difference.
elemental function FD2_C2(u &
                         ,ur &
                         ,ul) result(d)
  real(dp), intent(in) :: u
  real(dp), intent(in) :: ur
  real(dp), intent(in) :: ul
  real(dp) :: d
  d = &
    +1.0_dp*ur &
    -2.0_dp*u  &
    +1.0_dp*ul
end function FD2_C2
!! Fourth order accracy centeral undivided second finite difference.
elemental function FD2_C4(u &
                         ,ur,urr &
                         ,ul,ull) result(d)
  real(dp), intent(in) :: u
  real(dp), intent(in) :: ur,urr
  real(dp), intent(in) :: ul,ull
  real(dp) :: d
  d = -1.0_dp/12.0_dp*urr &
      +4.0_dp/03.0_dp*ur  &
      -5.0_dp/02.0_dp*u   &
      +4.0_dp/03.0_dp*ul  &
      -1.0_dp/12.0_dp*ull
end function FD2_C4
!! Sixth order accracy centeral undivided second finite difference.
elemental function FD2_C6(u &
                         ,ur,urr,urrr &
                         ,ul,ull,ulll) result(d)
  real(dp), intent(in) :: u
  real(dp), intent(in) :: ur,urr,urrr
  real(dp), intent(in) :: ul,ull,ulll
  real(dp) :: d
  d = +01.0_dp/90.0_dp*urrr &
      -03.0_dp/20.0_dp*urr  &
      +03.0_dp/02.0_dp*ur   &
      -49.0_dp/18.0_dp*u    &
      +03.0_dp/02.0_dp*ul   &
      -03.0_dp/20.0_dp*ull  &
      +01.0_dp/90.0_dp*ulll
end function FD2_C6
!! Eighth order accracy centeral undivided second finite difference.
elemental function FD2_C8(u &
                         ,ur,urr,urrr,urrrr &
                         ,ul,ull,ulll,ullll) result(d)
  real(dp), intent(in) :: u
  real(dp), intent(in) :: ur,urr,urrr,urrrr
  real(dp), intent(in) :: ul,ull,ulll,ullll
  real(dp) :: d
  d = -001.0_dp/560.0_dp*urrrr &
      +008.0_dp/315.0_dp*urrr  &
      -001.0_dp/005.0_dp*urr   &
      +008.0_dp/005.0_dp*ur    &
      -205.0_dp/072.0_dp*u     &
      +008.0_dp/005.0_dp*ul    &
      -001.0_dp/005.0_dp*ull   &
      +008.0_dp/315.0_dp*ulll  &
      -001.0_dp/560.0_dp*ullll
end function FD2_C8

!! -----------------------------------------------------------------  
!! The FDM-approximate Laplacian: v ← v + λΔu.
subroutine FDM_Laplacian(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: v(:),u(:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: kDir
  integer :: iCell 
  integer :: rCell,rrCell,rrrCell,rrrrCell
  integer :: lCell,llCell,lllCell,llllCell 
  associate(numDir=>mesh%NumDir &
         ,numCells=>mesh%NumCells &
       ,cellToCell=>mesh%CellToCell &
         ,invDxSqr=>lambda/mesh%Dx**2)
    !$omp parallel do firstprivate(rCell)
    do iCell = 1, numCells
      do kDir = 1, numDir, 2
        ! ----------------------
        ! Second order accracy.
        rCell=cellToCell(iCell,kDir)
        lCell=cellToCell(iCell,kDir+1)
        if (ACCURACY_ORDER <= 2) then
          v(iCell) = v(iCell) &
            + invDxSqr(kDir)*FD2_C2(u(iCell)&
                                   ,u(rCell),u(lCell))
        else
          ! ----------------------
          ! Fourth order accracy.
          rrCell=cellToCell(rCell,kDir)
          llCell=cellToCell(lCell,kDir+1)
          if (ACCURACY_ORDER <= 4) then
            v(iCell) = v(iCell) &
              + invDxSqr(kDir)*FD2_C4(u(iCell)&
                                     ,u(rCell),u(rrCell)&
                                     ,u(lCell),u(llCell))
          else
            ! ----------------------
            ! Sixth order accracy.
            rrrCell=cellToCell(rrCell,kDir)
            lllCell=cellToCell(llCell,kDir+1)
            if (ACCURACY_ORDER <= 6) then
              v(iCell) = v(iCell) &
                + invDxSqr(kDir)*FD2_C6(u(iCell)&
                                       ,u(rCell),u(rrCell),u(rrrCell)&
                                       ,u(lCell),u(llCell),u(lllCell))
            else ! 6
              ! ----------------------
              ! Eighth order accracy.
              rrrrCell=cellToCell(rrrCell,kDir)
              llllCell=cellToCell(lllCell,kDir+1)
              v(iCell) = v(iCell) &
                + invDxSqr(kDir)*FD2_C8(u(iCell)&
                                       ,u(rCell),u(rrCell),u(rrrCell),u(rrrrCell)&
                                       ,u(lCell),u(llCell),u(lllCell),u(llllCell))
            end if ! 8 
          end if ! 4
        end if ! 2
      end do
    end do
    !$omp end parallel do
  end associate
end subroutine FDM_Laplacian

!! -----------------------------------------------------------------  
!! The FDM-approximate nonlinear Laplacian: v ← v + λΔf(u).
subroutine FDM_LaplacianF(mesh,v,lambda,f,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: v(:),u(:)
  procedure(MathFunc) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: g(:)
  allocate(g, mold=u)
  ! ----------------------
  ! g ← f(u),
  ! v ← v + λΔg.
  call ApplyFunc(mesh,g,u,f)
  call FDM_Laplacian(mesh,v,lambda,g)
end subroutine FDM_LaplacianF

!! -----------------------------------------------------------------  
!! The FDM-approximate Bilaplacian: v ← v + λΔ²u.
subroutine FDM_Bilaplacian(mesh,v,lambda,u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: v(:),u(:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: w(:)
  allocate(w, mold=u)
  ! ----------------------
  ! w ← 0,
  ! w ← w + Δu.
  ! v ← v + λΔw.
  call Zero(mesh,w)
  call FDM_Laplacian(mesh,w,1.0_dp,u)
  call FDM_Laplacian(mesh,v,lambda,w)
end subroutine FDM_Bilaplacian

end module StormRuler_FiniteDifferences
  
  