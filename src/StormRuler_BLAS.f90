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

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Parameters, only: gCylCoords

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Dot
  module procedure Dot
end interface Dot

interface Norm_1
  module procedure Norm_1
end interface Norm_1

interface Norm_2
  module procedure Norm_2
end interface Norm_2

interface Norm_C
  module procedure Norm_C
end interface Norm_C

interface Fill
  module procedure Fill
end interface Fill

interface Fill_Random
  module procedure Fill_Random
end interface Fill_Random

interface Set
  module procedure Set
end interface Set

interface Scale
  module procedure Scale
end interface Scale

interface Add
  module procedure Add
end interface Add

interface Sub
  module procedure Sub
end interface Sub

interface Mul
  module procedure Mul
end interface Mul

interface SpFuncProd
  module procedure SpFuncProd
end interface SpFuncProd

interface Integrate
  module procedure Integrate
end interface Integrate

interface FuncProd
  module procedure FuncProd
end interface FuncProd

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Mathematical function: 𝑦 ← 𝑓(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  pure function tMapFunc(x) result(y)
    import dp
    real(dp), intent(in) :: x(:)
    real(dp) :: y(size(x))
  end function tMapFunc
end interface

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

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Matrix-vector product function: 𝒚 ← 𝓐(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tMatVecFunc(mesh, yArr, xArr)
    import :: tMesh, tArray
    class(tMesh), intent(in), target :: mesh
    class(tArray), intent(inout), target :: xArr, yArr
  end subroutine tMatVecFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Matrix-vector product function with a parameter: 𝒚 ← 𝓙(𝒙,𝒙̃).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tBiMatVecFunc(mesh, yArr, xArr, xTildeArr)
    import :: tMesh, tArray
    class(tMesh), intent(in), target :: mesh
    class(tArray), intent(inout), target :: xArr, xTildeArr, yArr
  end subroutine tBiMatVecFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute dot product: 𝑑 ← <𝒙⋅𝒚> = 𝒙ᵀ𝒚.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Dot(mesh, xArr, yArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr, yArr

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)
  
  Dot = mesh%RunCellKernel_Sum(Dot_Kernel)

contains
  real(dp) function Dot_Kernel(cell)
    integer(ip), intent(in) :: cell

    ! ----------------------
    ! 𝑑 ← 𝑑 + 𝒙ᵢ𝒚ᵢ. 
    ! ----------------------
    Dot_Kernel = sum(x(:,cell)*y(:,cell))
    if (gCylCoords) Dot_Kernel = Dot_Kernel*mesh%CellCenter(1,cell)
    
  end function Dot_Kernel
end function Dot

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ₁-norm: 𝑑 ← ‖𝒙‖₁.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Norm_1(mesh, xArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr

  real(dp), pointer :: x(:,:)

  call xArr%Get(x)

  Norm_1 = mesh%RunCellKernel_Sum(Norm_1_Kernel)

contains
  real(dp) function Norm_1_Kernel(cell)
    integer(ip), intent(in) :: cell

    Norm_1_Kernel = sum(abs(x(:,cell)))
    if (gCylCoords) Norm_1_Kernel = Norm_1_Kernel*mesh%CellCenter(1,cell)

  end function Norm_1_Kernel
end function Norm_1

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ₂-norm: 𝑑 ← ‖𝒙‖₂.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Norm_2(mesh, xArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr

  Norm_2 = sqrt( Dot(mesh, xArr, xArr) )

end function Norm_2

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ∞-norm: 𝑑 ← ‖𝒙‖∞.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Norm_C(mesh, xArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr

  real(dp), pointer :: x(:,:)

  call xArr%Get(x)

  Norm_C = mesh%RunCellKernel_Max(Norm_C_Kernel)

contains
  real(dp) function Norm_C_Kernel(cell)
    integer(ip), intent(in) :: cell

    Norm_C_Kernel = maxval(abs(x(:,cell)))
    
  end function Norm_C_Kernel
end function Norm_C

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components: 𝒚 ← 𝛼, 𝛼 ∊ ℝ.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Fill(mesh, yArr, alpha)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(inout) :: yArr
  real(dp), intent(in) :: alpha

  real(dp), pointer :: y(:,:)

  call yArr%Get(y)

  call mesh%RunCellKernel(Fill_Kernel)

contains
  subroutine Fill_Kernel(cell)
    integer(ip), intent(in) :: cell

    y(:,cell) = alpha
    
  end subroutine Fill_Kernel
end subroutine Fill

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components randomly: 𝒚 ← {𝛼ᵢ}ᵀ, where: 
!! • 𝛼ᵢ ~ 𝘜(𝑎,𝑏), 𝒚 ∊ ℝⁿ,
!! • 𝛼ᵢ = 𝜙ᵢ + 𝑐⋅𝘦𝘹𝘱(𝑖𝜓ᵢ), 𝜙ᵢ ~ 𝘜(𝑎,𝑏), 𝜓ᵢ ~ 𝘜(0,2𝜋), 𝒚 ∊ ℂⁿ.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Fill_Random(mesh, yArr, a, b)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(inout) :: yArr
  real(dp), intent(in), optional :: a, b

  real(dp), pointer :: y(:,:)

  call yArr%Get(y)

  call mesh%RunCellKernel_Block(Fill_Random_Kernel)

contains
  subroutine Fill_Random_Kernel(firstCell, lastCell)
    integer(ip), intent(in) :: firstCell, lastCell

    integer :: cell

    do cell = firstCell, lastCell
      call random_number(y(:,cell))
      if (present(a).and.present(b)) then
        y(:,cell) = min(a, b) + abs(b - a)*y(:,cell)
      end if
    end do
    
  end subroutine Fill_Random_Kernel
end subroutine Fill_Random

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒚 ← 𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Set(mesh, yArr, xArr)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  
  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(Set_Kernel)

contains
  subroutine Set_Kernel(cell)
    integer(ip), intent(in) :: cell

    y(:,cell) = x(:,cell)

  end subroutine Set_Kernel
end subroutine Set

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Scale: 𝒚 ← 𝛼𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Scale(mesh, yArr, xArr, alpha)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  real(dp), intent(in) :: alpha

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(Scale_Kernel)

contains
  subroutine Scale_Kernel(cell)
    integer(ip), intent(in) :: cell

    y(:,cell) = alpha*x(:,cell)

  end subroutine Scale_Kernel
end subroutine Scale

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: 𝒛 ← [[𝛽]]𝒚 + [𝛼]𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Add(mesh, zArr, yArr, xArr, alpha, beta)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr, yArr
  class(tArray), intent(inout) :: zArr
  real(dp), intent(in), optional :: alpha, beta

  real(dp), pointer :: x(:,:), y(:,:), z(:,:)
  real(dp) :: a, b

  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present( beta)) b =  beta

  call xArr%Get(x); call yArr%Get(y); call zArr%Get(z)

  call mesh%RunCellKernel(Add_Kernel)

contains
  subroutine Add_Kernel(cell)
    integer(ip), intent(in) :: cell

    z(:,cell) = b*y(:,cell) + a*x(:,cell)

  end subroutine Add_Kernel
end subroutine Add

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: 𝒛 ← 𝛽𝒚 - 𝛼𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Sub(mesh, zArr, yArr, xArr, alpha, beta)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr, yArr
  class(tArray), intent(inout) :: zArr
  real(dp), intent(in), optional :: alpha, beta

  real(dp), pointer :: x(:,:), y(:,:), z(:,:)
  real(dp) :: a, b

  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present( beta)) b =  beta

  call xArr%Get(x); call yArr%Get(y); call zArr%Get(z)

  call mesh%RunCellKernel(Sub_Kernel)

contains
  subroutine Sub_Kernel(cell)
    integer(ip), intent(in) :: cell

    z(:,cell) = b*y(:,cell) - a*x(:,cell)
    
  end subroutine Sub_Kernel
end subroutine Sub

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
!! Compute integral average: 𝑖 ← ∫𝑓(𝒙(𝒓))𝑑𝛺/∫1𝑑𝛺.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Integrate(mesh, xArr, f) result(integral)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr
  procedure(tMapFunc) :: f
  
  real(dp), pointer :: x(:,:)

  call xArr%Get(x)

  integral = mesh%RunCellKernel_Sum(Integrate_Kernel)/mesh%NumCells

contains
  real(dp) function Integrate_Kernel(cell)
    integer(ip), intent(in) :: cell

    associate(y => f(x(:,cell)))
      Integrate_Kernel = y(1)
    end associate
    if (gCylCoords) Integrate_Kernel = Integrate_Kernel*mesh%CellCenter(1,cell)

  end function Integrate_Kernel
end function Integrate

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: 𝒚 ← 𝑓(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FuncProd(mesh, yArr, xArr, f)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  procedure(tMapFunc) :: f

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(FuncProd_Kernel)
  
contains
  subroutine FuncProd_Kernel(cell)
    integer(ip), intent(in) :: cell

    y(:,cell) = f(x(:,cell))

  end subroutine FuncProd_Kernel
end subroutine FuncProd

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
