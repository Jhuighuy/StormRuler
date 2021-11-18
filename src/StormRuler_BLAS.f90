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

use StormRuler_Parameters, only: dp, ip, gCylCoords
use StormRuler_Helpers, only: Re, Im, R2C

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

interface MatVecProd_Diagonal
  module procedure MatVecProd_Diagonal
end interface MatVecProd_Diagonal

interface MatVecProd_Triangular
  module procedure MatVecProd_Triangular
end interface MatVecProd_Triangular

interface Solve_Triangular
  module procedure Solve_Triangular
end interface Solve_Triangular

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Mathematical function: ℳ𝒙 ← 𝑓(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  pure function tMapFunc(x) result(Mx)
    import dp
    real(dp), intent(in) :: x(:)
    real(dp) :: Mx(size(x))
  end function tMapFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Mathematical function: ℳ𝒙 ← 𝑓(𝒓,𝒙), 𝒓 ∊ 𝛺.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  pure function tSpMapFunc(r, x) result(Mx)
    import dp
    real(dp), intent(in) :: r(:)
    real(dp), intent(in) :: x(:)
    real(dp) :: Mx(size(x))
  end function tSpMapFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Matrix-vector product function: 𝓐𝒙 ← 𝓐(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tMatVecFunc(mesh, Ax, x)
    import :: tMesh, tArray
    class(tMesh), intent(inout), target :: mesh
    class(tArray), intent(inout), target :: x, Ax
  end subroutine tMatVecFunc
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Matrix-vector product function with a parameter: 𝓙𝒙 ← 𝓙(𝒙,𝒙₀).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tBiMatVecFunc(mesh, Jx, x, x0)
    import :: tMesh, tArray
    class(tMesh), intent(inout), target :: mesh
    class(tArray), intent(inout), target :: x, x0, Jx
  end subroutine tBiMatVecFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute dot product: 𝑑 ← <𝒙⋅𝒚> = 𝒙ᵀ𝒚.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Dot(mesh, xArr, yArr)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr, yArr

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)
  
  Dot = mesh%RunCellKernel_Sum(Dot_Kernel)

contains
  real(dp) function Dot_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    ! ----------------------
    ! 𝑑 ← 𝑑 + 𝒙ᵢ𝒚ᵢ. 
    ! ----------------------
    Dot_Kernel = sum(x(:,iCell)*y(:,iCell))
    if (gCylCoords) Dot_Kernel = Dot_Kernel*mesh%CellCenter(1,iCell)
    
  end function Dot_Kernel
end function Dot

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ₁-norm: 𝑑 ← ‖𝒙‖₁.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Norm_1(mesh, xArr)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr

  real(dp), pointer :: x(:,:)

  call xArr%Get(x)

  Norm_1 = mesh%RunCellKernel_Sum(Norm_1_Kernel)

contains
  real(dp) function Norm_1_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    Norm_1_Kernel = sum(abs(x(:,iCell)))
    if (gCylCoords) Norm_1_Kernel = Norm_1_Kernel*mesh%CellCenter(1,iCell)

  end function Norm_1_Kernel
end function Norm_1

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ₂-norm: 𝑑 ← ‖𝒙‖₂.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Norm_2(mesh, xArr)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr

  Norm_2 = sqrt( Dot(mesh, xArr, xArr) )

end function Norm_2

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ∞-norm: 𝑑 ← ‖𝒙‖∞.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Norm_C(mesh, xArr)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr

  real(dp), pointer :: x(:,:)

  call xArr%Get(x)

  Norm_C = mesh%RunCellKernel_Max(Norm_C_Kernel)

contains
  real(dp) function Norm_C_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    Norm_C_Kernel = maxval(abs(x(:,iCell)))
    
  end function Norm_C_Kernel
end function Norm_C

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components: 𝒚 ← 𝛼 + [𝛽], 𝛼 ∊ ℝ, [𝛽 ∊ ℝ or ℂ].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Fill(mesh, yArr, alpha, beta)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(inout) :: yArr
  real(dp), intent(in) :: alpha
  real(dp), intent(in), optional :: beta

  real(dp), pointer :: y(:,:)
  real(dp) :: gamma

  gamma = alpha
  if (present(beta)) gamma = gamma + beta

  call yArr%Get(y)

  call mesh%RunCellKernel(Fill_Kernel)

contains
  subroutine Fill_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = gamma
    
  end subroutine Fill_Kernel
end subroutine Fill

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components randomly: 𝒚 ← {𝛼ᵢ}ᵀ, where: 
!! • 𝛼ᵢ ~ 𝘜(𝑎,𝑏), 𝒚 ∊ ℝⁿ,
!! • 𝛼ᵢ = 𝜙ᵢ + 𝑐⋅𝘦𝘹𝘱(𝑖𝜓ᵢ), 𝜙ᵢ ~ 𝘜(𝑎,𝑏), 𝜓ᵢ ~ 𝘜(0,2𝜋), 𝒚 ∊ ℂⁿ.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Fill_Random(mesh, yArr, a, b)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(inout) :: yArr
  real(dp), intent(in), optional :: a, b

  real(dp), pointer :: y(:,:)

  call yArr%Get(y)

  ! TODO: not very parallel..
  call mesh%SetRange(parallel=.false.)
  call mesh%RunCellKernel_Block(Fill_Random_Kernel)
  call mesh%SetRange()

contains
  subroutine Fill_Random_Kernel(mesh, firstCell, lastCell)
    class(tMesh), intent(inout), target :: mesh
    integer(ip), intent(in) :: firstCell, lastCell

    integer :: iCell

    do iCell = firstCell, lastCell
      call random_number(y(:,iCell))
      if (present(a).and.present(b)) then
        y(:,iCell) = min(a, b) + abs(b - a)*y(:,iCell)
      end if
    end do
    
  end subroutine Fill_Random_Kernel
end subroutine Fill_Random

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒚 ← 𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Set(mesh, yArr, xArr)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  
  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(Set_Kernel)

contains
  subroutine Set_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = x(:,iCell)

  end subroutine Set_Kernel
end subroutine Set

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Scale: 𝒚 ← 𝛼𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Scale(mesh, yArr, xArr, alpha)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  real(dp), intent(in) :: alpha

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(Scale_Kernel)

contains
  subroutine Scale_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = alpha*x(:,iCell)

  end subroutine Scale_Kernel
end subroutine Scale

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: 𝒛 ← [[𝛽]]𝒚 + [𝛼]𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Add(mesh, zArr, yArr, xArr, alpha, beta)
  class(tMesh), intent(inout) :: mesh
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
  subroutine Add_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,iCell) = b*y(:,iCell) + a*x(:,iCell)

  end subroutine Add_Kernel
end subroutine Add

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: 𝒛 ← 𝛽𝒚 - 𝛼𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Sub(mesh, zArr, yArr, xArr, alpha, beta)
  class(tMesh), intent(inout) :: mesh
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
  subroutine Sub_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,iCell) = b*y(:,iCell) - a*x(:,iCell)
    
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
  class(tMesh), intent(inout) :: mesh
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
  subroutine MulScal_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,iCell) = yScal(iCell)*x(:,iCell)
    
  end subroutine MulScal_Kernel
  subroutine Mul_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,iCell) = yDiag(:,iCell)*x(:,iCell)
    
  end subroutine Mul_Kernel
  subroutine MulMat_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,iCell) = matmul(yMat(:,:,iCell), x(:,iCell))
    
  end subroutine MulMat_Kernel
end subroutine Mul

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute integral average: 𝑖 ← ∫𝑓(𝒙(𝒓))𝑑𝛺/∫1𝑑𝛺.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
real(dp) function Integrate(mesh, xArr, f) result(integral)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr
  procedure(tMapFunc) :: f
  
  real(dp), pointer :: x(:,:)

  call xArr%Get(x)

  integral = mesh%RunCellKernel_Sum(Integrate_Kernel)/mesh%NumCells

contains
  real(dp) function Integrate_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    associate(y => f(x(:,iCell)))
      Integrate_Kernel = y(1)
    end associate
    if (gCylCoords) Integrate_Kernel = Integrate_Kernel*mesh%CellCenter(1,iCell)

  end function Integrate_Kernel
end function Integrate

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: 𝒚 ← 𝑓(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FuncProd(mesh, yArr, xArr, f)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  procedure(tMapFunc) :: f

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(FuncProd_Kernel)
  
contains
  subroutine FuncProd_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = f(x(:,iCell))

  end subroutine FuncProd_Kernel
end subroutine FuncProd

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: 𝒚 ← 𝑓(𝒓,𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SpFuncProd(mesh, yArr, xArr, f)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: yArr
  procedure(tSpMapFunc) :: f

  real(dp), pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(SpFuncProd_Kernel)
  
contains
  subroutine SpFuncProd_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = f(mesh%CellCenter(iCell), x(:,iCell))

  end subroutine SpFuncProd_Kernel
end subroutine SpFuncProd

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Multiply a vector by diagonal of the matrix: 𝓓𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine MatVecProd_Diagonal(mesh, DxArr, xArr, MatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: DxArr
  procedure(tMatVecFunc) :: MatVec

  real(dp), pointer :: x(:,:), Dx(:,:)

  call xArr%Get(x); call DxArr%Get(Dx)

  call mesh%RunCellKernel_Block(MatVecProd_Diagonal_BlockKernel)
  call mesh%SetRange()

contains
  subroutine MatVecProd_Diagonal_BlockKernel(mesh, firstCell, lastCell)
    class(tMesh), intent(inout), target :: mesh
    integer(ip), intent(in) :: firstCell, lastCell

    integer(ip) :: iCell
    type(tArray) :: eArr
    real(dp), pointer :: e(:,:)

    call AllocArray(eArr, mold=xArr)
    call eArr%Get(e)

    e(:,:) = 0.0_dp

    do iCell = firstCell, lastCell
      e(:,iCell) = x(:,iCell)
      call mesh%SetRange(iCell)
      call MatVec(mesh, DxArr, eArr)
      e(:,iCell) = 0.0_dp
    end do

  end subroutine MatVecProd_Diagonal_BlockKernel  
end subroutine MatVecProd_Diagonal

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Multiply a vector by a lower/upper triangular 
!! part of the matrix: 𝓣𝒙 ← 𝘵𝘳𝘪𝘶(𝓐)𝒙 or 𝓣𝒙 ← 𝘵𝘳𝘪𝘭(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine MatVecProd_Triangular(mesh, upLo, TxArr, xArr, MatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: xArr
  class(tArray), intent(inout) :: TxArr
  procedure(tMatVecFunc) :: MatVec
  character, intent(in) :: upLo

  real(dp), pointer :: x(:,:), Tx(:,:)

  call xArr%Get(x); call TxArr%Get(Tx)

  if (upLo == 'U') then
    call mesh%RunCellKernel_Block(MatVecProd_UpperTriangular_BlockKernel)
  else if (upLo == 'L') then
    call mesh%RunCellKernel_Block(MatVecProd_LowerTriangular_BlockKernel)
  end if
  call mesh%SetRange()

contains
  subroutine MatVecProd_UpperTriangular_BlockKernel(mesh, firstCell, lastCell)
    class(tMesh), intent(inout), target :: mesh
    integer(ip), intent(in) :: firstCell, lastCell

    integer(ip) :: iCell
    type(tArray) :: eArr
    real(dp), pointer :: e(:,:)

    call AllocArray(eArr, mold=xArr)
    call eArr%Get(e)

    e(:,:firstCell-1) = 0.0_dp
    e(:,firstCell:) = x(:,firstCell:)

    do iCell = firstCell, lastCell
      call mesh%SetRange(iCell)
      call MatVec(mesh, TxArr, eArr)
      e(:,iCell) = 0.0_dp
    end do

  end subroutine MatVecProd_UpperTriangular_BlockKernel 
  subroutine MatVecProd_LowerTriangular_BlockKernel(mesh, firstCell, lastCell)
    class(tMesh), intent(inout), target :: mesh
    integer(ip), intent(in) :: firstCell, lastCell

    integer(ip) :: iCell
    type(tArray) :: eArr
    real(dp), pointer :: e(:,:)

    call AllocArray(eArr, mold=xArr)
    call eArr%Get(e)

    e(:,:lastCell) = x(:,:lastCell)
    e(:,lastCell+1:) = 0.0_dp
    
    do iCell = lastCell, firstCell, -1
      call mesh%SetRange(iCell)
      call MatVec(mesh, TxArr, eArr)
      e(:,iCell) = 0.0_dp
    end do

  end subroutine MatVecProd_LowerTriangular_BlockKernel 
end subroutine MatVecProd_Triangular

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve a linear system with a lower/upper 
!! triangular part of the matrix: 𝘵𝘳𝘪𝘶(𝓐)𝒙 = 𝒃 or 𝘵𝘳𝘪𝘭(𝓐)𝒙 = 𝒃.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Solve_Triangular(mesh, upLo, xArr, bArr, diagArr, MatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray), intent(in) :: bArr, diagArr
  class(tArray), intent(inout) :: xArr
  procedure(tMatVecFunc) :: MatVec
  character, intent(in) :: upLo

  real(dp), pointer :: x(:,:), b(:,:), diag(:,:)

  call xArr%Get(x); call bArr%Get(b); call diagArr%Get(diag)

  ! TODO: not very parallel..
  call mesh%SetRange(parallel=.false.)
  if (upLo == 'U') then
    call mesh%RunCellKernel_Block(Solve_UpperTriangular_BlockKernel)
  else if (upLo == 'L') then
    call mesh%RunCellKernel_Block(Solve_LowerTriangular_BlockKernel)
  end if
  call mesh%SetRange()

contains
  subroutine Solve_UpperTriangular_BlockKernel(mesh, firstCell, lastCell)
    class(tMesh), intent(inout), target :: mesh
    integer(ip), intent(in) :: firstCell, lastCell

    integer(ip) :: iCell

    type(tArray) :: AxArr
    real(dp), pointer :: Ax(:,:)

    call AllocArray(xArr, mold=AxArr)
    call AxArr%Get(Ax)

    Ax(:,:) = 0.0_dp
    x(:,:) = 0.0_dp
    x(:,lastCell) = b(:,lastCell)/diag(:,lastCell)

    do iCell = lastCell - 1, firstCell, -1
      call mesh%SetRange(iCell)
      call MatVec(mesh, AxArr, xArr)
      ! TODO: this is not a correct diagonal solution in block case!
      x(:,iCell) = (b(:,iCell) - Ax(:,iCell))/diag(:,iCell)
    end do

  end subroutine Solve_UpperTriangular_BlockKernel 
  subroutine Solve_LowerTriangular_BlockKernel(mesh, firstCell, lastCell)
    class(tMesh), intent(inout), target :: mesh
    integer(ip), intent(in) :: firstCell, lastCell

    integer(ip) :: iCell
    type(tArray) :: AxArr
    real(dp), pointer :: Ax(:,:)

    call AllocArray(xArr, mold=AxArr)
    call AxArr%Get(Ax)

    Ax(:,:) = 0.0_dp
    x(:,:) = 0.0_dp
    x(:,firstCell) = b(:,firstCell)/diag(:,firstCell)

    do iCell = firstCell + 1, lastCell
      call mesh%SetRange(iCell)
      call MatVec(mesh, AxArr, xArr)
      ! TODO: this is not a correct diagonal solution in block case!
      x(:,iCell) = (b(:,iCell) - Ax(:,iCell))/diag(:,iCell)
    end do

  end subroutine Solve_LowerTriangular_BlockKernel 
end subroutine Solve_Triangular

end module StormRuler_BLAS
