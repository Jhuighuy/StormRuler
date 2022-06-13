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
module StormRuler_BLAS

use StormRuler_Consts, only: ip, dp
use StormRuler_Parameters, only: gCylCoords

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Mul
  module procedure Mul
end interface Mul

interface SpFuncProd
  module procedure SpFuncProd
end interface SpFuncProd

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Mathematical function: 𝑦 ← 𝑓(𝒓,𝒙), 𝒓 ∊ 𝛺.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  pure function tSpMapFunc(r, x) result(y)
    import dp
    real(dp), intent(in) :: r(:), x(:)
    real(dp) :: y(size(x))
  end function tSpMapFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a product: 𝒛 ← 𝒚𝒙.
!! • Scalar multiplier case:
!!   Shape of 𝒙, 𝒛 is [1, NumVars]×[1, NumAllCells],
!!   Shape of 𝒚 is [1, NumAllCells].
!! • Diagonal multiplier case:
!!   Shape of 𝒙, 𝒚, 𝒛 is [1, NumVars]×[1, NumAllCells].
!! • Matrix multiplier case:
!!   Shape of 𝒙, 𝒛 is [1, NumVars]×[1, NumAllCells],
!!   Shape of 𝒚 is [1, NumVars]×[1, NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Mul(mesh, zArr, yArr, xArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr, yArr
  class(tArray), intent(inout) :: zArr

  real(dp), pointer :: x(:,:), z(:,:)
  real(dp), pointer :: yScal(:), yDiag(:,:), yMat(:,:,:)

  call xArr%Get(x); call zArr%Get(z)
  if (yArr%Rank() == xArr%Rank() - 1) then
    call yArr%Get(yScal)

    call mesh%RunCellKernel(MulScal_Kernel)

  else if (yArr%Rank() == xArr%Rank()) then
    call yArr%Get(yDiag)

    call mesh%RunCellKernel(Mul_Kernel)
  
  else if (yArr%Rank() == xArr%Rank() + 1) then
    call yArr%Get(yMat)

    call mesh%RunCellKernel(MulMat_Kernel)
  
  end if

contains
  subroutine MulScal_Kernel(cell)
    integer(ip), intent(in) :: cell

    z(:,cell) = yScal(cell)*x(:,cell)
    
  end subroutine MulScal_Kernel
  subroutine Mul_Kernel(cell)
    integer(ip), intent(in) :: cell

    z(:,cell) = yDiag(:,cell)*x(:,cell)
    
  end subroutine Mul_Kernel
  subroutine MulMat_Kernel(cell)
    integer(ip), intent(in) :: cell

    z(:,cell) = matmul(yMat(:,:,cell), x(:,cell))
    
  end subroutine MulMat_Kernel
end subroutine Mul

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: 𝒚 ← 𝑓(𝒓,𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SpFuncProd(mesh, yArr, xArr, f)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  procedure(tSpMapFunc) :: f

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(SpFuncProd_Kernel)
  
contains
  subroutine SpFuncProd_Kernel(cell)
    integer(ip), intent(in) :: cell

    y(:,cell) = f(mesh%CellCenter(cell), x(:,cell))

  end subroutine SpFuncProd_Kernel
end subroutine SpFuncProd

end module StormRuler_BLAS
