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

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8
use StormRuler_Helpers, only: Assert, Flip, SafeInverse, &
  & Reshape2D, Reshape3D, Reshape4D
use StormRuler_Helpers, only: operator(.inner.), operator(.outer.)
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: Fill, Mul_Outer, FuncProd
use StormRuler_FDM_Base

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

integer(ip), parameter :: FDM_AccuracyOrder = 2

logical, parameter :: FDM_CylCoords = .false.

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate gradient: 𝒗 ← 𝒗 - 𝜆∇𝒖.
!! Shape of 𝒖 is [1, NumVars]×[1, NumAllCells],
!! shape of 𝒗 is [1, Dim]×[1, NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Gradient_Central(mesh, vAny, lambda, uAny, &
    &                           dirAll, dirFace, dirCellFace)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), target :: uAny(:,:)
  real(dp), intent(inout), target :: vAny(:,:,:)
  integer(i8), intent(in), optional :: dirAll, dirFace(:), dirCellFace(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: u(:,:), vVec(:,:,:)

  u => uAny; vVec => vAny

  if (any([ present(dirAll), present(dirFace), present(dirCellFace) ])) then
    call mesh%RunCellKernel(FDM_Gradient_Forward_Kernel)
  else
    call mesh%RunCellKernel(FDM_Gradient_Central_Kernel)
  end if

contains
  subroutine FDM_Gradient_Central_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>
    
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell

    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    do iCellFace = 1, mesh%NumCellFaces, 2

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
        &       lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (FDM_AccuracyOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (FDM_AccuracyOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (FDM_AccuracyOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate gradient increment.
      ! ----------------------
      associate(dr_inv => lambda*SafeInverse(mesh%dr(:,iCellFace)))
        select case(FDM_AccuracyOrder)
          case(1:2)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              &  ( dr_inv.outer.FD1_C2(u(:,lCell), &
              &                        u(:,rCell)) )
          ! ----------------------
          case(3:4)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & ( dr_inv.outer.FD1_C4(u(:,llCell), &
              &                       u(:, lCell), &
              &                       u(:, rCell), &
              &                       u(:,rrCell)) )
          ! ----------------------
          case(5:6)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & ( dr_inv.outer.FD1_C6(u(:,lllCell), &
              &                       u(:, llCell), &
              &                       u(:,  lCell), &
              &                       u(:,  rCell), &
              &                       u(:, rrCell), &
              &                       u(:,rrrCell)) )
          ! ----------------------
          case(7:8)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & ( dr_inv.outer.FD1_C8(u(:,llllCell), &
              &                       u(:, lllCell), &
              &                       u(:,  llCell), &
              &                       u(:,   lCell), &
              &                       u(:,   rCell), &
              &                       u(:,  rrCell), &
              &                       u(:, rrrCell), &
              &                       u(:,rrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_Gradient_Central_Kernel
  subroutine FDM_Gradient_Forward_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell, rrrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell, lllllCell
    integer(i8) :: dir
    
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    do iCellFace = 1, mesh%NumCellFaces, 2

      ! ----------------------
      ! Determine FD direction (default is forward).
      ! ----------------------
      dir = 1_i8
      if (present(dirAll)) dir = dirAll
      if (present(dirFace)) dir = dirFace(iCellFace)
      if (present(dirCellFace)) dir = dirCellFace(iCellFace, iCell)

      ! ----------------------
      ! Find indices of the adjacent cells using the FD direction.
      ! ----------------------
      associate(inc => (1_i8-dir)/2_i8)
        associate(rCellFace => iCellFace+inc, &
          &       lCellFace => Flip(iCellFace+inc))
          rCell = mesh%CellToCell(rCellFace, iCell)
          lCell = mesh%CellToCell(lCellFace, iCell)
          if (FDM_AccuracyOrder >= 2) then
            rrCell = mesh%CellToCell(rCellFace, rCell)
            llCell = mesh%CellToCell(lCellFace, lCell)
            if (FDM_AccuracyOrder >= 4) then
              rrrCell = mesh%CellToCell(rCellFace, rrCell)
              lllCell = mesh%CellToCell(lCellFace, llCell)
              if (FDM_AccuracyOrder >= 6) then
                rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
                llllCell = mesh%CellToCell(lCellFace, lllCell)
                if (FDM_AccuracyOrder >= 8) then
                  rrrrrCell = mesh%CellToCell(rCellFace, rrrrCell)
                  lllllCell = mesh%CellToCell(lCellFace, llllCell)
                end if
              end if
            end if
          end if
        end associate
      end associate

      ! ----------------------
      ! Compute FDM-approximate gradient increment.
      ! ----------------------
      associate(dr_inv => lambda*SafeInverse(mesh%dr(:,iCellFace)))
        select case(FDM_AccuracyOrder)
          case(1)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & dir*( dr_inv.outer.FD1_F1(u(:,iCell), &
              &                           u(:,rCell)) )
          ! ----------------------
          case(2)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & dir*( dr_inv.outer.FD1_F2(u(:, iCell), &
              &                           u(:, rCell), &
              &                           u(:,rrCell)) )
          ! ----------------------
          case(3)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & dir*( dr_inv.outer.FD1_F3(u(:, lCell), &
              &                           u(:, iCell), &
              &                           u(:, rCell), &
              &                           u(:,rrCell)) )
          ! ----------------------
          case(4)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & dir*( dr_inv.outer.FD1_F4(u(:,  lCell), &
              &                           u(:,  iCell), &
              &                           u(:,  rCell), &
              &                           u(:, rrCell), &
              &                           u(:,rrrCell)) )
          ! ----------------------
          case(5)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & dir*( dr_inv.outer.FD1_F5(u(:, llCell), &
              &                           u(:,  lCell), &
              &                           u(:,  iCell), &
              &                           u(:,  rCell), &
              &                           u(:, rrCell), &
              &                           u(:,rrrCell)) )
          ! ----------------------
          case(6)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & dir*( dr_inv.outer.FD1_F6(u(:,  llCell), &
              &                           u(:,   lCell), &
              &                           u(:,   iCell), &
              &                           u(:,   rCell), &
              &                           u(:,  rrCell), &
              &                           u(:, rrrCell), &
              &                           u(:,rrrrCell)) )
          ! ----------------------
          case(7)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & dir*( dr_inv.outer.FD1_F7(u(:, lllCell), &
              &                           u(:,  llCell), &
              &                           u(:,   lCell), &
              &                           u(:,   iCell), &
              &                           u(:,   rCell), &
              &                           u(:,  rrCell), &
              &                           u(:, rrrCell), &
              &                           u(:,rrrrCell)) )
          ! ----------------------
          case(8)
            vVec(:,:,iCell) = vVec(:,:,iCell) - &
              & dir*( dr_inv.outer.FD1_F8(u(:,  lllCell), &
              &                           u(:,   llCell), &
              &                           u(:,    lCell), &
              &                           u(:,    iCell), &
              &                           u(:,    rCell), &
              &                           u(:,   rrCell), &
              &                           u(:,  rrrCell), &
              &                           u(:, rrrrCell), &
              &                           u(:,rrrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_Gradient_Forward_Kernel
end subroutine FDM_Gradient_Central

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate divergence: 𝒗 ← 𝒗 - 𝜆∇⋅𝒖.
!! Shape of 𝒖 is [1,Dim]×[1,NumVars]×[1, NumAllCells],
!! shape of 𝒗 is [1,NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Divergence_Central(mesh, vAny, lambda, uAny, &
    &                             dirAll, dirFace, dirCellFace)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), target :: uAny(:,:,:)
  real(dp), intent(inout), target :: vAny(:,:) 
  integer(i8), intent(in), optional :: dirAll, dirFace(:), dirCellFace(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: uVec(:,:,:), v(:,:)

  uVec => uAny; v => vAny

  if (any([ present(dirAll), present(dirFace), present(dirCellFace) ])) then
    call mesh%RunCellKernel(FDM_Divergence_Backward_Kernel)
  else
    call mesh%RunCellKernel(FDM_Divergence_Central_Kernel)
  end if

contains
  subroutine FDM_Divergence_Central_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>
    
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell

    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    do iCellFace = 1, mesh%NumCellFaces, 2

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
          &     lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (FDM_AccuracyOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (FDM_AccuracyOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (FDM_AccuracyOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate divergence increment.
      ! ----------------------
      associate(dr_inv => lambda*SafeInverse(mesh%dr(:,iCellFace)))

        ! ----------------------
        ! Cylindrical case: 
        ! (1/𝜌)∂(𝜌𝒖₁)/∂𝜌 component, force second order.
        ! ----------------------
        if (FDM_CylCoords.and.(iCellFace == 1)) then
          associate( &
            & rho_l => mesh%CellCenter(1,lCell), &
            & rho_i => mesh%CellCenter(1,iCell), &
            & rho_r => mesh%CellCenter(1,rCell))
            v(:,iCell) = v(:,iCell) - &
              & ( dr_inv.inner.FD1_C2(rho_l*uVec(:,:,lCell), &
              &                       rho_r*uVec(:,:,rCell)) &
              & )/rho_i
          end associate; cycle
        end if

        ! ----------------------
        ! General case.
        ! ----------------------
        select case(FDM_AccuracyOrder)
          case(1:2)
            v(:,iCell) = v(:,iCell) - &
              & ( dr_inv.inner.FD1_C2(uVec(:,:,lCell), &
              &                       uVec(:,:,rCell)) )
          ! ----------------------
          case(3:4)
            v(:,iCell) = v(:,iCell) - &
              & ( dr_inv.inner.FD1_C4(uVec(:,:,llCell), &
              &                       uVec(:,:, lCell), &
              &                       uVec(:,:, rCell), &
              &                       uVec(:,:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(:,iCell) = v(:,iCell) - &
              & ( dr_inv.inner.FD1_C6(uVec(:,:,lllCell), &
              &                       uVec(:,:, llCell), &
              &                       uVec(:,:,  lCell), &
              &                       uVec(:,:,  rCell), &
              &                       uVec(:,:, rrCell), &
              &                       uVec(:,:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(:,iCell) = v(:,iCell) - &
              & ( dr_inv.inner.FD1_C8(uVec(:,:,llllCell), &
              &                       uVec(:,:, lllCell), &
              &                       uVec(:,:,  llCell), &
              &                       uVec(:,:,   lCell), &
              &                       uVec(:,:,   rCell), &
              &                       uVec(:,:,  rrCell), &
              &                       uVec(:,:, rrrCell), &
              &                       uVec(:,:,rrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_Divergence_Central_Kernel
  subroutine FDM_Divergence_Backward_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell, rrrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell, lllllCell
    integer(i8) :: dir
    
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    do iCellFace = 1, mesh%NumCellFaces, 2

      ! ----------------------
      ! Determine FD direction (default is backward).
      ! ----------------------
      dir = -1_i8
      if (present(dirAll)) dir = dirAll
      if (present(dirFace)) dir = dirFace(iCellFace)
      if (present(dirCellFace)) dir = dirCellFace(iCellFace, iCell)

      ! ----------------------
      ! Find indices of the adjacent cells using the FD direction.
      ! ----------------------
      associate(inc => (1_i8-dir)/2_i8)
        associate(rCellFace => iCellFace+inc, &
          &       lCellFace => Flip(iCellFace+inc))
          rCell = mesh%CellToCell(rCellFace, iCell)
          lCell = mesh%CellToCell(lCellFace, iCell)
          if (FDM_AccuracyOrder >= 2) then
            rrCell = mesh%CellToCell(rCellFace, rCell)
            llCell = mesh%CellToCell(lCellFace, lCell)
            if (FDM_AccuracyOrder >= 4) then
              rrrCell = mesh%CellToCell(rCellFace, rrCell)
              lllCell = mesh%CellToCell(lCellFace, llCell)
              if (FDM_AccuracyOrder >= 6) then
                rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
                llllCell = mesh%CellToCell(lCellFace, lllCell)
                if (FDM_AccuracyOrder >= 8) then
                  rrrrrCell = mesh%CellToCell(rCellFace, rrrrCell)
                  lllllCell = mesh%CellToCell(lCellFace, llllCell)
                end if
              end if
            end if
          end if
        end associate
      end associate

      ! ----------------------
      ! Compute FDM-approximate divergence increment.
      ! ----------------------
      associate(dr_inv => lambda*SafeInverse(mesh%dr(:,iCellFace)))

        ! ----------------------
        ! Cylindrical case: 
        ! (1/𝜌)∂(𝜌𝒖₁)/∂𝜌 component, force first order.
        ! ----------------------
        if (FDM_CylCoords.and.(iCellFace == 1)) then
          associate( &
            & rho_i => mesh%CellCenter(1, iCell), &
            & rho_r => mesh%CellCenter(1, rCell) )
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F1(rho_i*uVec(:,:,iCell), &
              &                           rho_r*uVec(:,:,rCell)) &
              &     )/rho_i
          end associate; cycle
        end if

        ! ----------------------
        ! General case.
        ! ----------------------
        select case(FDM_AccuracyOrder)
          case(1)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F1(uVec(:,:,iCell), &
              &                           uVec(:,:,rCell)) )
          ! ----------------------
          case(2)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F2(uVec(:,:, iCell), &
              &                           uVec(:,:, rCell), &
              &                           uVec(:,:,rrCell)) )
          ! ----------------------
          case(3)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F3(uVec(:,:, lCell), &
              &                           uVec(:,:, iCell), &
              &                           uVec(:,:, rCell), &
              &                           uVec(:,:,rrCell)) )
          ! ----------------------
          case(4)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F4(uVec(:,:,  lCell), &
              &                           uVec(:,:,  iCell), &
              &                           uVec(:,:,  rCell), &
              &                           uVec(:,:, rrCell), &
              &                           uVec(:,:,rrrCell)) )
          ! ----------------------
          case(5)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F5(uVec(:,:, llCell), &
              &                           uVec(:,:,  lCell), &
              &                           uVec(:,:,  iCell), &
              &                           uVec(:,:,  rCell), &
              &                           uVec(:,:, rrCell), &
              &                           uVec(:,:,rrrCell)) )
          ! ----------------------
          case(6)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F6(uVec(:,:,  llCell), &
              &                           uVec(:,:,   lCell), &
              &                           uVec(:,:,   iCell), &
              &                           uVec(:,:,   rCell), &
              &                           uVec(:,:,  rrCell), &
              &                           uVec(:,:, rrrCell), &
              &                           uVec(:,:,rrrrCell)) )
          ! ----------------------
          case(7)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F7(uVec(:,:, lllCell), &
              &                           uVec(:,:,  llCell), &
              &                           uVec(:,:,   lCell), &
              &                           uVec(:,:,   iCell), &
              &                           uVec(:,:,   rCell), &
              &                           uVec(:,:,  rrCell), &
              &                           uVec(:,:, rrrCell), &
              &                           uVec(:,:,rrrrCell)) )
          ! ----------------------
          case(8)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dr_inv.inner.FD1_F8(uVec(:,:,  lllCell), &
              &                           uVec(:,:,   llCell), &
              &                           uVec(:,:,    lCell), &
              &                           uVec(:,:,    iCell), &
              &                           uVec(:,:,    rCell), &
              &                           uVec(:,:,   rrCell), &
              &                           uVec(:,:,  rrrCell), &
              &                           uVec(:,:, rrrrCell), &
              &                           uVec(:,:,rrrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_Divergence_Backward_Kernel
end subroutine FDM_Divergence_Central

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate convection: 𝒗 ← 𝒗 - 𝜆∇⋅𝒂𝒖.
!! Shape of 𝒖, 𝒗 is [1, NumVars]×[1, NumAllCells],
!! shape of 𝒂 is [1, Dim]×[1, NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Convection_Central(mesh, v, lambda, u, a)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), target :: u(:,:), a(:,:)
  real(dp), intent(inout) :: v(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), allocatable, target :: q(:,:,:)
  allocate(q(size(a, dim=1), size(u, dim=1), size(u, dim=2)))

  ! ----------------------
  ! Fast exit in case 𝜆 ≡ 0.
  ! ----------------------
  if (lambda == 0.0_dp) then
    return
  end if

  ! ----------------------
  ! 𝒒 ← 𝒂𝒖,
  ! 𝒗 ← 𝒗 - 𝜆∇⋅𝒒.
  ! ----------------------
  call Mul_Outer(mesh, q, a, u)
  call FDM_Divergence_Central(mesh, v, lambda, q)
end subroutine FDM_Convection_Central

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate Laplacian: 𝒗 ← 𝒗 + 𝜆Δ𝒖.
!! • Scalar case:
!!   Shape of 𝒖, 𝒗 is [1, NumVars]×[1, NumAllCells].
!! • Vector case:
!!   Shape of 𝒖, 𝒗 is [1, Dim]×[1, NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Laplacian_Central(mesh, vAny, lambda, uAny)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda
  real(dp), intent(in), target :: uAny(:,:)
  real(dp), intent(inout), target :: vAny(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: u(:,:), v(:,:)

  u => uAny; v => vAny

  call mesh%RunCellKernel(FDM_Laplacian_Central_Kernel)

contains
  subroutine FDM_Laplacian_Central_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell
    
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    do iCellFace = 1, mesh%NumCellFaces, 2

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
        &       lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (FDM_AccuracyOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (FDM_AccuracyOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (FDM_AccuracyOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate Laplacian increment.
      ! ----------------------
      associate(dl_sqr_inv => lambda/(mesh%dl(iCellFace)**2))

        ! ----------------------
        ! Cylindrical case: 
        ! (1/𝜌)∂(𝜌∂𝒖/∂𝜌) component, force second order.
        ! ----------------------
        if (FDM_CylCoords.and.(iCellFace == 1)) then
          associate( &
            & rho_l => mesh%CellCenter(1,lCell), &
            & rho_i => mesh%CellCenter(1,iCell), &
            & rho_r => mesh%CellCenter(1,rCell))
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * WFD2_C2(rho_l, u(:,lCell), &
              &                        rho_i, u(:,iCell), &
              &                        rho_r, u(:,rCell)) &
              & )/rho_i
          end associate; cycle
        end if

        ! ----------------------
        ! General case.
        ! ----------------------
        select case(FDM_AccuracyOrder)
          case(1:2)
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * FD2_C2(u(:,lCell), &
              &                       u(:,iCell), &
              &                       u(:,rCell)) )
          ! ----------------------
          case(3:4)
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * FD2_C4(u(:,llCell), &
              &                       u(:, lCell), &
              &                       u(:, iCell), &
              &                       u(:, rCell), &
              &                       u(:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * FD2_C6(u(:,lllCell), &
              &                       u(:, llCell), &
              &                       u(:,  lCell), &
              &                       u(:,  iCell), &
              &                       u(:,  rCell), &
              &                       u(:, rrCell), &
              &                       u(:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * FD2_C8(u(:,llllCell), &
              &                       u(:, lllCell), &
              &                       u(:,  llCell), &
              &                       u(:,   lCell), &
              &                       u(:,   iCell), &
              &                       u(:,   rCell), &
              &                       u(:,  rrCell), &
              &                       u(:, rrrCell), &
              &                       u(:,rrrrCell)) )
        end select
      end associate
    end do

#$if False
    if (FDM_CylCoords) then
      ! ----------------------
      ! Cylindrical coordinates, vector case.
      ! We have already computed:
      ! 𝒗 ← 𝒗 + 𝜆{Δ𝒖₁, Δ𝒖₂}ᵀ, 
      ! but we need:
      ! 𝒗 ← 𝒗 + 𝜆{Δ𝒖₁ - 𝒖₁/𝑟², Δ𝒖₂}ᵀ.
      ! The correction term is:
      ! 𝒗₁ ← 𝒗₁ - 𝜆𝒖₁/𝑟².
      ! ----------------------
      associate(r => mesh%CellCenter(iCell))
        vVec(1,:,iCell) = vVec(1,:,iCell) - &
          & lambda*uVec(1,:,iCell)/( r(1)**2 )
      end associate
    end if
#$end if
  end subroutine FDM_Laplacian_Central_Kernel
end subroutine FDM_Laplacian_Central

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate variable coefficient Laplacian: 𝒗 ← 𝒗 + 𝜆∇⋅(𝒌∇𝒖).
!! Shape of 𝒖, 𝒗 is [1, NumVars]×[1, NumAllCells].
!! • Scalar coefficient case:
!!   Shape of 𝒌 is [1, NumVars]×[1, NumAllCells].
!! • Tensor coefficient case (not implemented):
!!   Shape of 𝒌 is [1, Dim]×[1, Dim]×[1, NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_DivWGrad_Central(mesh, vAny, lambda, kAny, uAny)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda 
  real(dp), intent(in), target :: uAny(:,:), kAny(:,:)
  real(dp), intent(inout), target :: vAny(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: u(:,:), v(:,:), k(:,:)

  u => uAny; v => vAny; k => kAny

  call mesh%RunCellKernel(FDM_DivWGrad_Central_Kernel)

contains
  subroutine FDM_DivWGrad_Central_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell
  
    ! ----------------------
    ! For each positive cell face do:
    ! ----------------------
    do iCellFace = 1, mesh%NumCellFaces, 2

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
        &       lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (FDM_AccuracyOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (FDM_AccuracyOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (FDM_AccuracyOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate variable coefficient Laplacian increment.
      ! ----------------------
      associate(dl_sqr_inv => lambda/(mesh%dl(iCellFace)**2))

        ! ----------------------
        ! Cylindrical case,
        ! (1/𝜌)∂(𝜌𝒘∂𝒖/∂𝜌) component, force second order.
        ! ----------------------
        if (FDM_CylCoords.and.(iCellFace == 1)) then
          associate( &
            & rho_l => mesh%CellCenter(1,lCell), &
            & rho_i => mesh%CellCenter(1,iCell), &
            & rho_r => mesh%CellCenter(1,rCell))
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * &
              &   WFD2_C2(rho_l*k(:,lCell), u(:,lCell), &
              &           rho_i*k(:,iCell), u(:,iCell), &
              &           rho_r*k(:,rCell), u(:,rCell)) &
              & )/rho_i
          end associate; cycle
        end if

        ! ----------------------
        ! General case.
        ! ----------------------
        select case(FDM_AccuracyOrder)
          case(1:2)
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * &
              &   WFD2_C2(k(:,lCell), u(:,lCell), &
              &           k(:,iCell), u(:,iCell), &
              &           k(:,rCell), u(:,rCell)) )
          ! ----------------------
          case(3:4)
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * &
              &   WFD2_C4(k(:,llCell), u(:,llCell), &
              &           k(:, lCell), u(:, lCell), &
              &           k(:, iCell), u(:, iCell), &
              &           k(:, rCell), u(:, rCell), &
              &           k(:,rrCell), u(:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * &
              &   WFD2_C6(k(:,lllCell), u(:,lllCell), &
              &           k(:, llCell), u(:, llCell), &
              &           k(:,  lCell), u(:,  lCell), &
              &           k(:,  iCell), u(:,  iCell), &
              &           k(:,  rCell), u(:,  rCell), &
              &           k(:, rrCell), u(:, rrCell), &
              &           k(:,rrrCell), u(:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(:,iCell) = v(:,iCell) + &
              & ( dl_sqr_inv * &
              &   WFD2_C8(k(:,llllCell), u(:,llllCell), &
              &           k(:, lllCell), u(:, lllCell), &
              &           k(:,  llCell), u(:,  llCell), &
              &           k(:,   lCell), u(:,   lCell), &
              &           k(:,   iCell), u(:,   iCell), &
              &           k(:,   rCell), u(:,   rCell), &
              &           k(:,  rrCell), u(:,  rrCell), &
              &           k(:, rrrCell), u(:, rrrCell), &
              &           k(:,rrrrCell), u(:,rrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_DivWGrad_Central_Kernel
end subroutine FDM_DivWGrad_Central

end module StormRuler_FDM_Operators
