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

use StormRuler_Parameters, only: dp, ip, i8, gCylCoords
use StormRuler_Helpers, only: Flip

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArrayR

use StormRuler_FDM_Base, only: FD1_C2, FD1_C4, FD1_C6, FD1_C8, &
  & FD1_F1, FD1_F2, FD1_F3, FD1_F4, FD1_F5, FD1_F6, FD1_F7, FD1_F8, &
  & FD2_C2, FD2_C4, FD2_C6, FD2_C8, WFD2_C2, WFD2_C4, WFD2_C6, WFD2_C8

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

integer(ip), parameter :: gTruncErrorOrder = 1

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate gradient: 𝒗⃗ ← 𝒗⃗ - 𝜆∇𝒖.
!! Shape of 𝒖 is [1, NumVars]×[1, NumAllCells],
!! shape of 𝒗⃗ is [1, NumDims]×[1, NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Gradient_Central(mesh, vVecArr, lambda, uArr, &
    &                           dirAll, dirFace, dirCellFace)
  class(tMesh), intent(inout) :: mesh
  class(tArrayR), intent(in) :: uArr
  class(tArrayR), intent(inout) :: vVecArr
  real(dp), intent(in) :: lambda
  integer(i8), intent(in), optional :: dirAll, dirFace(:), dirCellFace(:,:)

  real(dp), pointer :: u(:,:), vVec(:,:,:)

  call uArr%Get(u); call vVecArr%Get(vVec)

  if (present(dirAll).or.present(dirFace).or.present(dirCellFace)) then
    call mesh%RunCellKernel(FDM_Gradient_Forward_Kernel)
  else
    call mesh%RunCellKernel(FDM_Gradient_Central_Kernel)
  end if

contains
  subroutine FDM_Gradient_Central_Kernel(iCell)
    integer(ip), intent(in) :: iCell
    
    integer(ip) :: dim
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell

    ! ----------------------
    ! For each direction do:
    ! ----------------------
    do dim = 1, mesh%NumDims
      iCellFace = 2*(dim - 1) + 1

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
        &       lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (gTruncErrorOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (gTruncErrorOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (gTruncErrorOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate gradient increment.
      ! ----------------------
      associate(dlInv => lambda/mesh%dl(iCellFace))
        select case(gTruncErrorOrder)
          case(1:2)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &          ( dlInv*FD1_C2(u(:,lCell), &
              &                         u(:,rCell)) )
          ! ----------------------
          case(3:4)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &         ( dlInv*FD1_C4(u(:,llCell), &
              &                        u(:, lCell), &
              &                        u(:, rCell), &
              &                        u(:,rrCell)) )
          ! ----------------------
          case(5:6)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &        ( dlInv*FD1_C6(u(:,lllCell), &
              &                       u(:, llCell), &
              &                       u(:,  lCell), &
              &                       u(:,  rCell), &
              &                       u(:, rrCell), &
              &                       u(:,rrrCell)) )
          ! ----------------------
          case(7:8)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &       ( dlInv*FD1_C8(u(:,llllCell), &
              &                      u(:, lllCell), &
              &                      u(:,  llCell), &
              &                      u(:,   lCell), &
              &                      u(:,   rCell), &
              &                      u(:,  rrCell), &
              &                      u(:, rrrCell), &
              &                      u(:,rrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_Gradient_Central_Kernel
  subroutine FDM_Gradient_Forward_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    integer(i8) :: dir
    integer(ip) :: dim
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell, rrrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell, lllllCell
    
    ! ----------------------
    ! For each direction do:
    ! ----------------------
    do dim = 1, mesh%NumDims
      iCellFace = 2*(dim - 1) + 1

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
          if (gTruncErrorOrder >= 2) then
            rrCell = mesh%CellToCell(rCellFace, rCell)
            llCell = mesh%CellToCell(lCellFace, lCell)
            if (gTruncErrorOrder >= 4) then
              rrrCell = mesh%CellToCell(rCellFace, rrCell)
              lllCell = mesh%CellToCell(lCellFace, llCell)
              if (gTruncErrorOrder >= 6) then
                rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
                llllCell = mesh%CellToCell(lCellFace, lllCell)
                if (gTruncErrorOrder >= 8) then
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
      associate(dlInv => lambda/mesh%dl(iCellFace))
        select case(gTruncErrorOrder)
          case(1)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &      dir*( dlInv*FD1_F1(u(:,iCell), &
              &                         u(:,rCell)) )
          ! ----------------------
          case(2)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &     dir*( dlInv*FD1_F2(u(:, iCell), &
              &                        u(:, rCell), &
              &                        u(:,rrCell)) )
          ! ----------------------
          case(3)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &     dir*( dlInv*FD1_F3(u(:, lCell), &
              &                        u(:, iCell), &
              &                        u(:, rCell), &
              &                        u(:,rrCell)) )
          ! ----------------------
          case(4)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &    dir*( dlInv*FD1_F4(u(:,  lCell), &
              &                       u(:,  iCell), &
              &                       u(:,  rCell), &
              &                       u(:, rrCell), &
              &                       u(:,rrrCell)) )
          ! ----------------------
          case(5)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &    dir*( dlInv*FD1_F5(u(:, llCell), &
              &                       u(:,  lCell), &
              &                       u(:,  iCell), &
              &                       u(:,  rCell), &
              &                       u(:, rrCell), &
              &                       u(:,rrrCell)) )
          ! ----------------------
          case(6)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &   dir*( dlInv*FD1_F6(u(:,  llCell), &
              &                      u(:,   lCell), &
              &                      u(:,   iCell), &
              &                      u(:,   rCell), &
              &                      u(:,  rrCell), &
              &                      u(:, rrrCell), &
              &                      u(:,rrrrCell)) )
          ! ----------------------
          case(7)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &   dir*( dlInv*FD1_F7(u(:, lllCell), &
              &                      u(:,  llCell), &
              &                      u(:,   lCell), &
              &                      u(:,   iCell), &
              &                      u(:,   rCell), &
              &                      u(:,  rrCell), &
              &                      u(:, rrrCell), &
              &                      u(:,rrrrCell)) )
          ! ----------------------
          case(8)
            vVec(dim,:,iCell) = vVec(dim,:,iCell) - &
              &  dir*( dlInv*FD1_F8(u(:,  lllCell), &
              &                     u(:,   llCell), &
              &                     u(:,    lCell), &
              &                     u(:,    iCell), &
              &                     u(:,    rCell), &
              &                     u(:,   rrCell), &
              &                     u(:,  rrrCell), &
              &                     u(:, rrrrCell), &
              &                     u(:,rrrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_Gradient_Forward_Kernel
end subroutine FDM_Gradient_Central

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate divergence: 𝒗 ← 𝒗 - 𝜆∇⋅𝒖⃗.
!! Shape of 𝒖⃗ is [1,NumDims]×[1,NumVars]×[1, NumAllCells],
!! shape of 𝒗 is [1,NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Divergence_Central(mesh, vArr, lambda, uVecArr, &
    &                             dirAll, dirFace, dirCellFace)
  class(tMesh), intent(inout) :: mesh
  class(tArrayR), intent(in) :: uVecArr
  class(tArrayR), intent(inout) :: vArr
  real(dp), intent(in) :: lambda
  integer(i8), intent(in), optional :: dirAll, dirFace(:), dirCellFace(:,:)

  real(dp), pointer :: uVec(:,:,:), v(:,:)

  call uVecArr%Get(uVec); call vArr%Get(v)

  if (present(dirAll).or.present(dirFace).or.present(dirCellFace)) then
    call mesh%RunCellKernel(FDM_Divergence_Backward_Kernel)
  else
    call mesh%RunCellKernel(FDM_Divergence_Central_Kernel)
  end if

contains
  subroutine FDM_Divergence_Central_Kernel(iCell)
    integer(ip), intent(in) :: iCell
    
    integer(ip) :: dim
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell

    ! ----------------------
    ! For each direction do:
    ! ----------------------
    do dim = 1, mesh%NumDims
      iCellFace = 2*(dim - 1) + 1

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
          &     lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (gTruncErrorOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (gTruncErrorOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (gTruncErrorOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate divergence increment.
      ! ----------------------
      associate(dlInv => lambda/mesh%dl(iCellFace))

        ! ----------------------
        ! Cylindrical case: 
        ! (1/𝜌)∂(𝜌𝒖₁)/∂𝜌 component, force second order.
        ! ----------------------
        if (gCylCoords.and.(dim == 1)) then
          associate( &
            & rho_l => mesh%CellCenter(1,lCell), &
            & rho_i => mesh%CellCenter(1,iCell), &
            & rho_r => mesh%CellCenter(1,rCell))
            v(:,iCell) = v(:,iCell) - &
              & ( dlInv*FD1_C2(rho_l*uVec(1,:,lCell), &
              &                rho_r*uVec(1,:,rCell)) )/rho_i
          end associate; cycle
        end if

        ! ----------------------
        ! General case.
        ! ----------------------
        select case(gTruncErrorOrder)
          case(1:2)
            v(:,iCell) = v(:,iCell) - &
              & ( dlInv*FD1_C2(uVec(dim,:,lCell), &
              &                uVec(dim,:,rCell)) )
          ! ----------------------
          case(3:4)
            v(:,iCell) = v(:,iCell) - &
              & ( dlInv*FD1_C4(uVec(dim,:,llCell), &
              &                uVec(dim,:, lCell), &
              &                uVec(dim,:, rCell), &
              &                uVec(dim,:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(:,iCell) = v(:,iCell) - &
              & ( dlInv*FD1_C6(uVec(dim,:,lllCell), &
              &                uVec(dim,:, llCell), &
              &                uVec(dim,:,  lCell), &
              &                uVec(dim,:,  rCell), &
              &                uVec(dim,:, rrCell), &
              &                uVec(dim,:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(:,iCell) = v(:,iCell) - &
              & ( dlInv*FD1_C8(uVec(dim,:,llllCell), &
              &                uVec(dim,:, lllCell), &
              &                uVec(dim,:,  llCell), &
              &                uVec(dim,:,   lCell), &
              &                uVec(dim,:,   rCell), &
              &                uVec(dim,:,  rrCell), &
              &                uVec(dim,:, rrrCell), &
              &                uVec(dim,:,rrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_Divergence_Central_Kernel
  subroutine FDM_Divergence_Backward_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    integer(i8) :: dir
    integer(ip) :: dim
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell, rrrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell, lllllCell
    
    ! ----------------------
    ! For each direction do:
    ! ----------------------
    do dim = 1, mesh%NumDims
      iCellFace = 2*(dim - 1) + 1

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
          if (gTruncErrorOrder >= 2) then
            rrCell = mesh%CellToCell(rCellFace, rCell)
            llCell = mesh%CellToCell(lCellFace, lCell)
            if (gTruncErrorOrder >= 4) then
              rrrCell = mesh%CellToCell(rCellFace, rrCell)
              lllCell = mesh%CellToCell(lCellFace, llCell)
              if (gTruncErrorOrder >= 6) then
                rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
                llllCell = mesh%CellToCell(lCellFace, lllCell)
                if (gTruncErrorOrder >= 8) then
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
      associate(dlInv => lambda/mesh%dl(iCellFace))

        ! ----------------------
        ! Cylindrical case: 
        ! (1/𝜌)∂(𝜌𝒖₁)/∂𝜌 component, force first order.
        ! ----------------------
        if (gCylCoords.and.(dim == 1)) then
          associate( &
            & rho_i => mesh%CellCenter(1, iCell), &
            & rho_r => mesh%CellCenter(1, rCell) )
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F1(rho_i*uVec(dim,:,iCell), &
              &                    rho_r*uVec(dim,:,rCell)) )/rho_i
          end associate; cycle
        end if

        ! ----------------------
        ! General case.
        ! ----------------------
        select case(gTruncErrorOrder)
          case(1)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F1(uVec(dim,:,iCell), &
              &                    uVec(dim,:,rCell)) )
          ! ----------------------
          case(2)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F2(uVec(dim,:, iCell), &
              &                    uVec(dim,:, rCell), &
              &                    uVec(dim,:,rrCell)) )
          ! ----------------------
          case(3)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F3(uVec(dim,:, lCell), &
              &                    uVec(dim,:, iCell), &
              &                    uVec(dim,:, rCell), &
              &                    uVec(dim,:,rrCell)) )
          ! ----------------------
          case(4)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F4(uVec(dim,:,  lCell), &
              &                    uVec(dim,:,  iCell), &
              &                    uVec(dim,:,  rCell), &
              &                    uVec(dim,:, rrCell), &
              &                    uVec(dim,:,rrrCell)) )
          ! ----------------------
          case(5)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F5(uVec(dim,:, llCell), &
              &                    uVec(dim,:,  lCell), &
              &                    uVec(dim,:,  iCell), &
              &                    uVec(dim,:,  rCell), &
              &                    uVec(dim,:, rrCell), &
              &                    uVec(dim,:,rrrCell)) )
          ! ----------------------
          case(6)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F6(uVec(dim,:,  llCell), &
              &                    uVec(dim,:,   lCell), &
              &                    uVec(dim,:,   iCell), &
              &                    uVec(dim,:,   rCell), &
              &                    uVec(dim,:,  rrCell), &
              &                    uVec(dim,:, rrrCell), &
              &                    uVec(dim,:,rrrrCell)) )
          ! ----------------------
          case(7)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F7(uVec(dim,:, lllCell), &
              &                    uVec(dim,:,  llCell), &
              &                    uVec(dim,:,   lCell), &
              &                    uVec(dim,:,   iCell), &
              &                    uVec(dim,:,   rCell), &
              &                    uVec(dim,:,  rrCell), &
              &                    uVec(dim,:, rrrCell), &
              &                    uVec(dim,:,rrrrCell)) )
          ! ----------------------
          case(8)
            v(:,iCell) = v(:,iCell) - &
              & dir*( dlInv*FD1_F8(uVec(dim,:,  lllCell), &
              &                    uVec(dim,:,   llCell), &
              &                    uVec(dim,:,    lCell), &
              &                    uVec(dim,:,    iCell), &
              &                    uVec(dim,:,    rCell), &
              &                    uVec(dim,:,   rrCell), &
              &                    uVec(dim,:,  rrrCell), &
              &                    uVec(dim,:, rrrrCell), &
              &                    uVec(dim,:,rrrrrCell)) )
        end select
      end associate
    end do
  end subroutine FDM_Divergence_Backward_Kernel
end subroutine FDM_Divergence_Central

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate Laplacian: 𝒗 ← 𝒗 + 𝜆Δ𝒖 or 𝒗⃗ ← 𝒗⃗ + 𝜆Δ𝒖⃗.
!! • Scalar case:
!!   Shape of 𝒖, 𝒗 is [1, NumVars]×[1, NumAllCells].
!! • Vector case:
!!   Shape of 𝒖⃗, 𝒗⃗ is [1, NumDims]×[1, NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Laplacian_Central(mesh, vArr, lambda, uArr)
  class(tMesh), intent(inout) :: mesh
  class(tArrayR), intent(in) :: uArr
  class(tArrayR), intent(inout) :: vArr
  real(dp), intent(in) :: lambda

  real(dp), pointer :: u(:,:), v(:,:)
  real(dp), pointer :: uVec(:,:,:), vVec(:,:,:)

  call uArr%Get(u); call vArr%Get(v)
  if (gCylCoords.and.(uArr%Rank() == 3)) then
    call uArr%Get(uVec); call vArr%Get(vVec)
  end if

  call mesh%RunCellKernel(FDM_Laplacian_Central_Kernel)

contains
  subroutine FDM_Laplacian_Central_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    integer(ip) :: dim
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell
    
    ! ----------------------
    ! For each direction do:
    ! ----------------------
    do dim = 1, mesh%NumDims
      iCellFace = 2*(dim - 1) + 1

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
        &       lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (gTruncErrorOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (gTruncErrorOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (gTruncErrorOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate Laplacian increment.
      ! ----------------------
      associate(dlSqrInv => lambda/(mesh%dl(iCellFace)**2))

        ! ----------------------
        ! Cylindrical case: 
        ! (1/𝜌)∂(𝜌∂𝒖/∂𝜌) component, force second order.
        ! ----------------------
        if (gCylCoords.and.(iCellFace == 1)) then
          associate( &
            & rho_l => mesh%CellCenter(1,lCell), &
            & rho_i => mesh%CellCenter(1,iCell), &
            & rho_r => mesh%CellCenter(1,rCell))
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*WFD2_C2(rho_l, u(:,lCell), &
              &                    rho_i, u(:,iCell), &
              &                    rho_r, u(:,rCell)) )/rho_i
          end associate; cycle
        end if

        ! ----------------------
        ! General case.
        ! ----------------------
        select case(gTruncErrorOrder)
          case(1:2)
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*FD2_C2(u(:,lCell), &
              &                   u(:,iCell), &
              &                   u(:,rCell)) )
          ! ----------------------
          case(3:4)
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*FD2_C4(u(:,llCell), &
              &                   u(:, lCell), &
              &                   u(:, iCell), &
              &                   u(:, rCell), &
              &                   u(:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*FD2_C6(u(:,lllCell), &
              &                   u(:, llCell), &
              &                   u(:,  lCell), &
              &                   u(:,  iCell), &
              &                   u(:,  rCell), &
              &                   u(:, rrCell), &
              &                   u(:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*FD2_C8(u(:,llllCell), &
              &                   u(:, lllCell), &
              &                   u(:,  llCell), &
              &                   u(:,   lCell), &
              &                   u(:,   iCell), &
              &                   u(:,   rCell), &
              &                   u(:,  rrCell), &
              &                   u(:, rrrCell), &
              &                   u(:,rrrrCell)) )
        end select
      end associate
    end do

    if (gCylCoords.and.(uArr%Rank() == 3)) then
      ! ----------------------
      ! Cylindrical coordinates, vector case.
      ! We have already computed:
      ! 𝒗⃗ ← 𝒗⃗ + 𝜆{Δ𝒖₁, Δ𝒖₂}ᵀ, 
      ! but we need:
      ! 𝒗⃗ ← 𝒗⃗ + 𝜆{Δ𝒖₁ - 𝒖₁/𝜌², Δ𝒖₂}ᵀ.
      ! ----------------------
      associate(rho => mesh%CellCenter(1,iCell))
        vVec(1,:,iCell) = vVec(1,:,iCell) - lambda*uVec(1,:,iCell)/( rho**2 )
      end associate
    end if
  end subroutine FDM_Laplacian_Central_Kernel
end subroutine FDM_Laplacian_Central

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The FDM-approximate variable coefficient 
!! Laplacian: 𝒗 ← 𝒗 + 𝜆∇⋅(𝒘∇𝒖) or 𝒗 ← 𝒗 + 𝜆∇⋅(̂𝒘∇𝒖).
!! Shape of 𝒖, 𝒗 is [1, NumVars]×[1, NumAllCells].
!! • Scalar coefficient case:
!!   Shape of 𝒘 is [1, NumVars]×[1, NumAllCells].
!! • Tensor coefficient case (not implemented):
!!   Shape of ̂𝒘 is [1, NumDims]×[1, NumDims]×[1, NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_DivWGrad_Central(mesh, vArr, lambda, wArr, uArr)
  class(tMesh), intent(inout) :: mesh
  class(tArrayR), intent(in) :: uArr, wArr
  class(tArrayR), intent(inout) :: vArr
  real(dp), intent(in) :: lambda

  real(dp), pointer :: u(:,:), v(:,:), w(:,:)

  call uArr%Get(u); call vArr%Get(v); call wArr%Get(w)

  call mesh%RunCellKernel(FDM_DivWGrad_Central_Kernel)

contains
  subroutine FDM_DivWGrad_Central_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer(ip) :: dim
    integer(ip) :: iCellFace
    integer(ip) :: rCell, rrCell, rrrCell, rrrrCell
    integer(ip) :: lCell, llCell, lllCell, llllCell
  
    ! ----------------------
    ! For each direction do:
    ! ----------------------
    do dim = 1, mesh%NumDims
      iCellFace = 2*(dim - 1) + 1

      ! ----------------------
      ! Find indices of the adjacent cells.
      ! ----------------------
      associate(rCellFace => iCellFace, &
        &       lCellFace => Flip(iCellFace))
        rCell = mesh%CellToCell(rCellFace, iCell)
        lCell = mesh%CellToCell(lCellFace, iCell)
        if (gTruncErrorOrder >= 3) then
          rrCell = mesh%CellToCell(rCellFace, rCell)
          llCell = mesh%CellToCell(lCellFace, lCell)
          if (gTruncErrorOrder >= 5) then
            rrrCell = mesh%CellToCell(rCellFace, rrCell)
            lllCell = mesh%CellToCell(lCellFace, llCell)
            if (gTruncErrorOrder >= 7) then
              rrrrCell = mesh%CellToCell(rCellFace, rrrCell)
              llllCell = mesh%CellToCell(lCellFace, lllCell)
            end if
          end if
        end if
      end associate

      ! ----------------------
      ! Compute FDM-approximate variable coefficient Laplacian increment.
      ! ----------------------
      associate(dlSqrInv => lambda/(mesh%dl(iCellFace)**2))

        ! ----------------------
        ! Cylindrical case,
        ! (1/𝜌)∂(𝜌𝒘∂𝒖/∂𝜌) component, force second order.
        ! ----------------------
        if (gCylCoords.and.(iCellFace == 1)) then
          associate( &
            & rho_l => mesh%CellCenter(1,lCell), &
            & rho_i => mesh%CellCenter(1,iCell), &
            & rho_r => mesh%CellCenter(1,rCell))
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*WFD2_C2(rho_l*w(:,lCell), u(:,lCell), &
              &                    rho_i*w(:,iCell), u(:,iCell), &
              &                    rho_r*w(:,rCell), u(:,rCell)) )/rho_i
          end associate; cycle
        end if

        ! ----------------------
        ! General case.
        ! ----------------------
        select case(gTruncErrorOrder)
          case(1:2)
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*WFD2_C2(w(:,lCell), u(:,lCell), &
              &                    w(:,iCell), u(:,iCell), &
              &                    w(:,rCell), u(:,rCell)) )
          ! ----------------------
          case(3:4)
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*WFD2_C4(w(:,llCell), u(:,llCell), &
              &                    w(:, lCell), u(:, lCell), &
              &                    w(:, iCell), u(:, iCell), &
              &                    w(:, rCell), u(:, rCell), &
              &                    w(:,rrCell), u(:,rrCell)) )
          ! ----------------------
          case(5:6)
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*WFD2_C6(w(:,lllCell), u(:,lllCell), &
              &                    w(:, llCell), u(:, llCell), &
              &                    w(:,  lCell), u(:,  lCell), &
              &                    w(:,  iCell), u(:,  iCell), &
              &                    w(:,  rCell), u(:,  rCell), &
              &                    w(:, rrCell), u(:, rrCell), &
              &                    w(:,rrrCell), u(:,rrrCell)) )
          ! ----------------------
          case(7:8)
            v(:,iCell) = v(:,iCell) + &
              & ( dlSqrInv*WFD2_C8(w(:,llllCell), u(:,llllCell), &
              &                    w(:, lllCell), u(:, lllCell), &
              &                    w(:,  llCell), u(:,  llCell), &
              &                    w(:,   lCell), u(:,   lCell), &
              &                    w(:,   iCell), u(:,   iCell), &
              &                    w(:,   rCell), u(:,   rCell), &
              &                    w(:,  rrCell), u(:,  rrCell), &
              &                    w(:, rrrCell), u(:, rrrCell), &
              &                    w(:,rrrrCell), u(:,rrrrCell)) )
        end select
      end associate
    end do

  end subroutine FDM_DivWGrad_Central_Kernel
end subroutine FDM_DivWGrad_Central

end module StormRuler_FDM_Operators
