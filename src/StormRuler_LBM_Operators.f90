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
module StormRuler_LBM_Operators

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

use StormRuler_Helpers, only: Flip

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray

use StormRuler_LBM_Base, only:

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface LBM_Stream
  module procedure LBM_Stream
end interface LBM_Stream

interface LBM_Macroscopics
  module procedure LBM_Macroscopics
end interface LBM_Macroscopics

interface LBM_Collision_BGK
  module procedure LBM_Collision_BGK
end interface LBM_Collision_BGK

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! LBM streaming phase: 𝒈ᵢ(𝒓 + 𝜏⋅𝒗ᵢ,𝑡 + 𝜏) ← 𝒇ᵢ(𝒓,𝑡). 
!! Shape of 𝒇, 𝒈 is [1, NumVars]×[1, NumCellConns]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LBM_Stream(mesh, gArr, fArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: fArr
  class(tArray), intent(inout) :: gArr

  real(dp), pointer :: f(:,:), g(:,:)
  logical, save :: first = .true.
  
  call fArr%Get(f); call gArr%Get(g)

  if (first) then
    first = .false.
    f(:,mesh%NumCells+1:) = 1.0
  end if

  call mesh%RunCellKernel(LBM_Stream_Kernel)

  block
    integer(ip) :: mark, cell, bndCell, bndCellConn

    do mark = 1, mesh%NumBndMarks
      do bndCell = mesh%BndCellAddrs(mark), mesh%BndCellAddrs(mark + 1) - 1

        if (mark == 1) then
          associate(ff => f(:,bndCell))
            f(:,bndCell) = [ff(2), ff(1), ff(4), ff(3), ff(6), ff(5), ff(8), ff(7), ff(9)] 
          end associate
        else
          bndCellConn = mesh%BndCellConns(bndCell)
          cell = mesh%CellToCell(Flip(bndCellConn), bndCell)
          f(:,bndCell) = f(:,cell)
        end if
  
        do bndCellConn = 1, mesh%NumCellConns

          ! ----------------------
          ! Index of the adjacent cell.
          ! ----------------------
          cell = mesh%CellToCell(bndCellConn, bndCell)
          if (cell == 0) cycle
    
          ! ----------------------
          ! Stream the distribution function component:
          ! ----------------------
          g(bndCellConn,cell) = f(bndCellConn,bndCell)
    
        end do

      end do
    end do

  end block
  
contains
  subroutine LBM_Stream_Kernel(cell)
    integer(ip), intent(in) :: cell

    integer(ip) :: cellConn
    integer(ip) :: cellCell

    ! ----------------------
    ! For each connection do:
    ! ----------------------
    do cellConn = 1, mesh%NumCellConns

      ! ----------------------
      ! Index of the adjacent cell.
      ! ----------------------
      cellCell = mesh%CellToCell(cellConn, cell)

      ! ----------------------
      ! Stream the distribution function component:
      ! ----------------------
      g(cellConn,cellCell) = f(cellConn,cell)

    end do

  end subroutine LBM_Stream_Kernel
end subroutine LBM_Stream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute LBM macroscopic variables: 
!! • Density 𝛒: 𝛒(𝒓,𝑡) ← 𝑚⋅∑𝒇ᵢ(𝒓,𝑡), 
!! • Velocity 𝒗: 𝛒(𝒓,𝑡)𝒗(𝒓,𝑡) ← 𝑚⋅∑𝒗ᵢ⋅𝒇ᵢ(𝒓,𝑡).
!!
!! Shape of 𝒇 is [1, NumVars]×[1, NumCellConns]×[1, NumAllCells].
!! Shape of 𝑚 is [1, NumVars].
!! Shape of 𝛒 is [1, NumAllCells].
!! Shape of 𝒗 is [1, NumDims]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LBM_Macroscopics(mesh, rhoArr, vArr, m, fArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: fArr
  class(tArray), intent(inout) :: rhoArr, vArr
  real(dp), intent(in) :: m(:)

  real(dp), pointer :: f(:,:), rho(:), v(:,:)

  call fArr%Get(f); call rhoArr%Get(rho); call vArr%Get(v)

  call mesh%RunCellKernel(LBM_Macroscopics_Kernel)

contains
  subroutine LBM_Macroscopics_Kernel(cell)
    integer(ip), intent(in) :: cell

    integer(ip) :: cellConn

    rho(cell) = 0.0_dp; v(:,cell) = 0.0_dp

    ! ----------------------
    ! For each connection do:
    ! ----------------------
    do cellConn = 1, mesh%NumCellConns

      ! ----------------------
      ! Compute macroscopic variables increment:
      ! ----------------------
      associate(mf => f(cellConn,cell))
        rho(cell) = rho(cell) + mf
        v(:,cell) = v(:,cell) + mf*mesh%dr(:,cellConn)
      end associate

    end do

    ! ----------------------
    ! Convert momentum to velocity:
    ! ----------------------
    v(:,cell) = v(:,cell)/rho(cell)

  end subroutine LBM_Macroscopics_Kernel
end subroutine LBM_Macroscopics

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute LBM collision integral, Bhatnagar Gross and Krook model: 
!! 𝒇ᵢ ← 𝒇ᵢ + 𝜴ᵢ[𝒇(𝒓,𝑡)], where 𝜴ᵢ[𝒇] = [𝓕ᵢ(𝛒(𝒓,𝑡),𝒗(𝒓,𝑡)) - 𝒇ᵢ(𝒓,𝑡)]/𝜏.
!!
!! Shape of 𝒇 is [1, NumVars]×[1, NumCellConns]×[1, NumAllCells].
!! Shape of 𝛒 is [1, NumAllCells].
!! Shape of 𝒗 is [1, NumDims]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LBM_Collision_BGK(mesh, fArr, tau, rhoArr, vArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: fArr
  class(tArray), intent(inout) :: rhoArr, vArr
  real(dp) :: tau

  !! TODO: select true `w` based on the mesh extended connectivity.
  real(dp), pointer :: f(:,:), rho(:), v(:,:)
  real(dp), parameter :: w(*) = [ &
    & 1/9.0_dp, 1/9.0_dp, 1/9.0_dp, 1/9.0_dp, &
    & 1/36.0_dp, 1/36.0_dp, 1/36.0_dp, 1/36.0_dp, 4/9.0_dp ]

  call fArr%Get(f); call rhoArr%Get(rho); call vArr%Get(v)

  call mesh%RunCellKernel(LBM_Collision_BGK_Kernel)
  
contains
  subroutine LBM_Collision_BGK_Kernel(cell)
    integer(ip), intent(in) :: cell

    integer(ip) :: cellConn
    real(dp) :: vSqr, fEq

    vSqr = dot_product(v(:,cell), v(:,cell))

    ! ----------------------
    ! For each connection do:
    ! ----------------------
    do cellConn = 1, mesh%NumCellConns

      ! ----------------------
      ! Compute equilibrium distribution: 
      ! 𝓕ᵢ ← 𝛒𝑤ᵢ⋅(1 + 3⋅𝒗ᵢ⋅𝒗 + 4.5⋅(𝒗ᵢ⋅𝒗)² - 1.5⋅(𝒗⋅𝒗)²).
      ! ----------------------
      associate(d => dot_product(v(:,cell), mesh%dr(:,cellConn)))
        fEq = rho(cell)*w(cellConn)*(1.0_dp + 3.0_dp*d + 4.5_dp*d**2 - 1.5_dp*vSqr)
      end associate

      ! ----------------------
      ! Update collistion integral:
      ! 𝒇ᵢ ← 𝒇ᵢ + [𝓕ᵢ - 𝒇ᵢ]/𝜏.
      ! ----------------------
      f(cellConn,cell) = f(cellConn,cell) + (fEq - f(cellConn,cell))/tau 

    end do

  end subroutine LBM_Collision_BGK_Kernel
end subroutine LBM_Collision_BGK

end module StormRuler_LBM_Operators
