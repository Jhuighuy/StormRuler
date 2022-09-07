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

use StormRuler_Consts, only: ip, dp

use StormRuler_Helpers, only: Flip

use StormRuler_Mesh, only: tMesh

#$use 'StormRuler_Macros.fi'

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
subroutine FDM_ApplyBCs$rank(mesh, mark, u, alpha, beta, gamma)!, f)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: mark
  real(dp), intent(in) :: alpha, beta, gamma
  real(dp), intent(inout) :: u(@:,:)
  !procedure(tSMapFuncR$rank), optional :: f
  
  integer(ip) :: cell, bndCell, bndCellConn, gCell

  associate(mLambda => (0.5_dp*alpha - beta/mesh%dl), &
    &    pLambdaInv => 1.0_dp/(0.5_dp*alpha + beta/mesh%dl))

    !$omp parallel do private(cell, bndCell, bndCellConn, gCell)
    do bndCell = mesh%BndCellAddrs(mark), mesh%BndCellAddrs(mark + 1) - 1
      bndCellConn = mesh%BndCellConns(bndCell)
      cell = mesh%CellToCell(Flip(bndCellConn), bndCell)
  
      u(@:,bndCell) = pLambdaInv(bndCellConn) * &
        & (gamma - mLambda(bndCellConn)*u(@:,cell))
  
      ! ----------------------
      ! Propagate the boundary condition towards the ghost cells.
      ! ----------------------
      gCell = mesh%CellToCell(bndCellConn, bndCell)
      do while(gCell /= 0)
        u(@:,gCell) = u(@:,bndCell)
        gCell = mesh%CellToCell(bndCellConn, gCell)
      end do
  
    end do
    !$omp end parallel do

  end associate
end subroutine FDM_ApplyBCs$rank
#$end do

subroutine FDM_ApplyBCs_SlipWall(mesh, mark, v)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: mark
  real(dp), intent(inout) :: v(:,:)

  integer(ip) :: dim
  integer(ip) :: cell, bndCell, bndCellConn, gCell

  !$omp parallel do private(cell, bndCell, bndCellConn, gCell)
  do bndCell = mesh%BndCellAddrs(mark), mesh%BndCellAddrs(mark + 1) - 1
    bndCellConn = mesh%BndCellConns(bndCell)
    cell = mesh%CellToCell(Flip(bndCellConn), bndCell)

    dim = (bndCellConn - 1)/2 + 1

    v(:,bndCell) = v(:,cell)  
    v(dim,bndCell) = -v(dim,cell)  

    ! ----------------------
    ! Propagate the boundary condition towards the ghost cells.
    ! ----------------------
    gCell = mesh%CellToCell(bndCellConn, bndCell)
    do while(gCell /= 0)
      v(:,gCell) = v(:,bndCell)
      gCell = mesh%CellToCell(bndCellConn, gCell)
    end do

  end do
  !$omp end parallel do

end subroutine FDM_ApplyBCs_SlipWall

subroutine FDM_ApplyBCs_CosWall(mesh, mark, phi_hat, phi, a)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: mark
  real(dp), intent(inout) :: phi_hat(:), phi(:)
  real(dp) :: a

  integer(ip) :: dim
  integer(ip) :: cell, bndCell, bndCellConn, gCell

  !$omp parallel do private(cell, bndCell, bndCellConn, gCell)
  do bndCell = mesh%BndCellAddrs(mark), mesh%BndCellAddrs(mark + 1) - 1
    bndCellConn = mesh%BndCellConns(bndCell)
    cell = mesh%CellToCell(Flip(bndCellConn), bndCell)

    dim = (bndCellConn - 1)/2 + 1

    !phi_hat(bndCell) = phi_hat(cell) - a*( phi(cell)*( 1.0_dp - phi(cell) ) ) 
    
    phi_hat(bndCell) = phi_hat(cell) - &
      & a*0.5_dp*( phi_hat(cell)*( 1.0_dp - phi(cell) ) + phi(cell)*( 1.0_dp - phi_hat(cell) ) )

    ! ----------------------
    ! Propagate the boundary condition towards the ghost cells.
    ! ----------------------
    gCell = mesh%CellToCell(bndCellConn, bndCell)
    do while(gCell /= 0)
      phi_hat(gCell) = phi_hat(bndCell)
      gCell = mesh%CellToCell(bndCellConn, gCell)
    end do

  end do
  !$omp end parallel do

end subroutine FDM_ApplyBCs_CosWall

subroutine FDM_ApplyBCs_InOutLet(mesh, mark, v, Q)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: mark
  real(dp), intent(in) :: Q !Flux on the boundary
  real(dp), intent(inout) :: v(:,:)

  integer(ip) :: dim
  integer(ip) :: cell, bndCell, bndCellConn, gCell

  real(dp) :: R, RR, RRR, Pi

  R = 0.0_dp

  !$omp parallel do  private(cell, bndCell, bndCellConn, gCell) reduction(max: R)
  do bndCell = mesh%BndCellAddrs(mark), mesh%BndCellAddrs(mark + 1) - 1
    bndCellConn = mesh%BndCellConns(bndCell)
    cell = mesh%CellToCell(Flip(bndCellConn), bndCell)

    dim = (bndCellConn - 1)/2 + 1

    R = max(R, mesh%CellCenter(1, cell))

  end do
  !$omp end parallel do

  R = R + 0.5_dp*mesh%dl(1)

  !$omp parallel do private(cell, bndCell, bndCellConn, gCell)
  do bndCell = mesh%BndCellAddrs(mark), mesh%BndCellAddrs(mark + 1) - 1
    bndCellConn = mesh%BndCellConns(bndCell)
    cell = mesh%CellToCell(Flip(bndCellConn), bndCell)

    RR = mesh%CellCenter(1, cell)

    !RRR = merge(0.1_dp, 10.0_dp, mark == 2)
    !RRR = merge(0.1_dp, 0.101_dp, mark == 2)
    Pi = 3.14159265358979323846
    RRR = 2.0_dp * Q / Pi / (R**4)

    !RRR = 10 = dp/dz / (4 eta)?

    v(:,bndCell) = 0.0_dp
    v(2,bndCell) = -RRR*( R**2 - RR**2 )

    ! ----------------------
    ! Propagate the boundary condition towards the ghost cells.
    ! ----------------------
    gCell = mesh%CellToCell(bndCellConn, bndCell)
    do while(gCell /= 0)
      v(:,gCell) = v(:,bndCell)
      gCell = mesh%CellToCell(bndCellConn, gCell)
    end do

  end do
  !$omp end parallel do

end subroutine FDM_ApplyBCs_InOutLet

end module StormRuler_FDM_BCs
