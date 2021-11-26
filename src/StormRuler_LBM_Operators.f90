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
!! Shape of 𝒇, 𝒈 is [1, NumVars]×[1, NumExtDirs]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LBM_Stream(mesh, gArr, fArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: fArr
  class(tArray), intent(inout) :: gArr

  real(dp), pointer :: f(:,:,:), g(:,:,:)

  call fArr%Get(f); call gArr%Get(g)

  call mesh%RunCellKernel(LBM_Stream_Kernel)

contains
  subroutine LBM_Stream_Kernel(cell)
    integer(ip), intent(in) :: cell

    integer(ip) :: extDir
    integer(ip) :: cellCell

    ! ----------------------
    ! For each extended direction do:
    ! ----------------------
    do extDir = 1, mesh%NumExtDirs - 1

      ! ----------------------
      ! Index of the adjacent cell.
      ! ----------------------
      cellCell = mesh%CellToCell(Flip(extDir), cell)

      ! ----------------------
      ! Stream the distribution function:
      ! ----------------------
      g(:,extDir,cell) = f(:,extDir,cellCell)

    end do

    ! ----------------------
    ! Stream the distribution function stationary component:
    ! ----------------------
    g(:,extDir,cell) = f(:,extDir,cell)

  end subroutine LBM_Stream_Kernel
end subroutine LBM_Stream

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute LBM macroscopic variables: 
!! • Density 𝛒: 𝛒(𝒓,𝑡) ← 𝑚⋅∑𝒇ᵢ(𝒓,𝑡), 
!! • Velocity 𝒗: 𝛒(𝒓,𝑡)𝒗(𝒓,𝑡) ← 𝑚⋅∑𝒗ᵢ⋅𝒇ᵢ(𝒓,𝑡).
!! Shape of 𝒇 is [1, NumVars]×[1, NumExtDirs]×[1, NumAllCells].
!! Shape of 𝑚 is [1, NumVars].
!! Shape of 𝛒 is [1, NumAllCells].
!! Shape of 𝒗 is [1, NumDims]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LBM_Macroscopics(mesh, rhoArr, m, vArr, fArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: fArr
  class(tArray), intent(inout) :: rhoArr, vArr
  real(dp), intent(in) :: m(:)

  real(dp), pointer :: f(:,:,:), rho(:), v(:,:)

  call fArr%Get(f); call rhoArr%Get(rho); call vArr%Get(v)

  call mesh%RunCellKernel(LBM_Macroscopics_Kernel)

contains
  subroutine LBM_Macroscopics_Kernel(cell)
    integer(ip), intent(in) :: cell

    integer(ip) :: extDir

    rho(cell) = 0.0_dp; v(:,cell) = 0.0_dp

    ! ----------------------
    ! For each extended direction do:
    ! ----------------------
    do extDir = 1, mesh%NumExtDirs - 1

      ! ----------------------
      ! Compute macroscopic variables increment:
      ! ----------------------
      associate(mf => dot_product(m, f(:,extDir,cell)))
        rho(cell) = rho(cell) + mf
        v(:,cell) = v(:,cell) + mf*mesh%dr(:,extDir)
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
!! 𝜴ᵢ ← 𝒇ᵢ + 𝜴ᵢ[𝒇(𝒓,𝑡)], where 𝜴[𝒇] = [𝓕ᵢ(𝛒(𝒓,𝑡),𝒗(𝒓,𝑡)) - 𝒇ᵢ(𝒓,𝑡)]/𝜏.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine LBM_Collision_BGK(mesh, omgArr, tau, fArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: fArr
  class(tArray), intent(inout) :: omgArr
  real(dp) :: tau

end subroutine LBM_Collision_BGK

end module StormRuler_LBM_Operators
