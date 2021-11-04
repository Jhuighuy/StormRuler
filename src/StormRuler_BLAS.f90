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
use StormRuler_Helpers, only: Re, Im, R2C, operator(.inner.), operator(.outer.)

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArrayR

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Dot
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Dot$T
#$end for
end interface Dot

interface Norm_1
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Norm_1$T
#$end for
end interface Norm_1

interface Norm_2
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Norm_2$T
#$end for
end interface Norm_2

interface Norm_C
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Norm_C$T
#$end for
end interface Norm_C

interface Fill
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Fill$T
#$end for
end interface Fill

interface Fill_Random
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Fill_Random$T
#$end for
end interface Fill_Random

interface Set
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Set$T
#$end for
end interface Set

interface Set_Real
  module procedure Set_Real
end interface Set_Real

interface Set_Imag
  module procedure Set_Imag
end interface Set_Imag

interface Set_Complex
  module procedure Set_Complex
end interface Set_Complex

interface Scale
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Scale$T
#$end for
end interface Scale

interface Add
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Add$T
#$end for
end interface Add

interface Sub
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Sub$T
#$end for
end interface Sub

interface Mul
#$do rank = 0, NUM_RANKS-3
  module procedure Mul$rank
#$end do
end interface Mul

interface Mul_Inner
#$do rank = 0, NUM_RANKS-3
  module procedure Mul_Inner$rank
#$end do
end interface Mul_Inner

interface Mul_Outer
#$do rank = 0, NUM_RANKS-3
  module procedure Mul_Outer$rank
#$end do
end interface Mul_Outer

interface SFuncProd
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure SFuncProd$T
#$end for
end interface SFuncProd

interface Integrate
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Integrate$T
#$end for
end interface Integrate

interface FuncProd
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure FuncProd$T
#$end for
end interface FuncProd

interface MatVecProd_Diagonal
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure MatVecProd_Diagonal$T
#$end for
end interface MatVecProd_Diagonal

interface MatVecProd_Triangular
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure MatVecProd_Triangular$T
#$end for
end interface MatVecProd_Triangular

interface Solve_Triangular
#$for T, _ in [SCALAR_TYPES[0]]
  module procedure Solve_Triangular$T
#$end for
end interface Solve_Triangular

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Mathematical function: ℳ𝒙 ← 𝑓(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for T, typename in [SCALAR_TYPES[0]]
  pure function tMapFunc$T(x) result(Mx)
    import dp
    $typename, intent(in) :: x(:)
    $typename :: Mx(size(x))
  end function tMapFunc$T
#$end for
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Mathematical function: ℳ𝒙 ← 𝑓(𝒓,𝒙), 𝒓 ∊ 𝛺.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for T, typename in SCALAR_TYPES
  pure function tSMapFunc$T(r, x) result(SMx)
    import dp
    real(dp), intent(in) :: r(:)
    $typename, intent(in) :: x(:)
    $typename :: SMx(size(x))
  end function tSMapFunc$T
#$end for
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Matrix-vector product function: 𝓐𝒙 ← 𝓐(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for T, typename in [SCALAR_TYPES[0]]
  subroutine tMatVecFunc$T(mesh, Ax, x)
    import :: tMesh, tArray$T
    class(tMesh), intent(inout), target :: mesh
    class(tArray$T), intent(inout), target :: x, Ax
  end subroutine tMatVecFunc$T
#$end for
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Matrix-vector product function with a parameter: 𝓙𝒙 ← 𝓙(𝒙,𝒙₀).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for T, typename in [SCALAR_TYPES[0]]
  subroutine tBiMatVecFunc$T(mesh, Jx, x, x0)
    import :: tMesh, tArray$T
    class(tMesh), intent(inout), target :: mesh
    class(tArray$T), intent(inout), target :: x, x0, Jx
  end subroutine tBiMatVecFunc$T
#$end for
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute dot product: 
!! • 𝑑 ← <𝒙⋅𝒚> = 𝒙ᴴ𝒚 (default), or 
!! • 𝑑 ← [𝒙⋅𝒚] = 𝒙ᵀ𝒚 (do_conjg = false).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
$typename function Dot$T(mesh, xArr, yArr, do_conjg) result(Dot)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr, yArr
  logical, intent(in), optional :: do_conjg

  $typename, pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)
  
#$if T == 'C'
  if (present(do_conjg)) then
    if (.not.do_conjg) then
      Dot = mesh%RunCellKernel_Sum(Dot_Kernel)
      return
    end if
  end if
  Dot = mesh%RunCellKernel_Sum(Dot_Conjg_Kernel)
#$else
  Dot = mesh%RunCellKernel_Sum(Dot_Kernel)
#$end if

contains
  $typename function Dot_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    ! ----------------------
    ! 𝗼𝘂𝘁 ← 𝗼𝘂𝘁 + 𝒙ᵢ𝒚ᵢ. 
    ! ----------------------
    Dot_Kernel = sum(x(:,iCell)*y(:,iCell))
    if (gCylCoords) Dot_Kernel = Dot_Kernel*mesh%CellCenter(1,iCell)
    
  end function Dot_Kernel
#$if T == 'C'
  $typename function Dot_Conjg_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    ! ----------------------
    ! 𝗼𝘂𝘁 ← 𝗼𝘂𝘁 + 𝒙̅ᵢ𝒚ᵢ.
    ! ----------------------
    Dot_Conjg_Kernel = sum(conjg(x(:,iCell))*y(:,iCell))
    if (gCylCoords) Dot_Conjg_Kernel = Dot_Conjg_Kernel*mesh%CellCenter(1,iCell)
    
  end function Dot_Conjg_Kernel
#$end if
end function Dot$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ₁-norm: 𝑑 ← ‖𝒙‖₁.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
real(dp) function Norm_1$T(mesh, xArr) result(Norm_1)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr

  $typename, pointer :: x(:,:)

  call xArr%Get(x)

  Norm_1 = mesh%RunCellKernel_Sum(Norm_1_Kernel)

contains
  real(dp) function Norm_1_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    Norm_1_Kernel = sum(abs(x(:,iCell)))
    if (gCylCoords) Norm_1_Kernel = Norm_1_Kernel*mesh%CellCenter(1,iCell)

  end function Norm_1_Kernel
end function Norm_1$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ₂-norm: 𝑑 ← ‖𝒙‖₂.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
real(dp) function Norm_2$T(mesh, xArr) result(Norm_2)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr

#$if T == 'C'
  Norm_2 = sqrt( Re(Dot(mesh, xArr, xArr)) )
#$else
  Norm_2 = sqrt( Dot(mesh, xArr, xArr) )
#$end if

end function Norm_2$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ∞-norm: 𝑑 ← ‖𝒙‖∞.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
real(dp) function Norm_C$T(mesh, xArr) result(Norm_C)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr

  $typename, pointer :: x(:,:)

  call xArr%Get(x)

  Norm_C = mesh%RunCellKernel_Max(Norm_C_Kernel)

contains
  real(dp) function Norm_C_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    Norm_C_Kernel = maxval(abs(x(:,iCell)))
    
  end function Norm_C_Kernel
end function Norm_C$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components: 𝒚 ← 𝛼 + [𝛽], 𝛼 ∊ ℝ, [𝛽 ∊ ℝ or ℂ].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Fill$T(mesh, yArr, alpha, beta)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(inout) :: yArr
  real(dp), intent(in) :: alpha
  $typename, intent(in), optional :: beta

  $typename, pointer :: y(:,:)
  $typename :: gamma

  gamma = alpha
  if (present(beta)) gamma = gamma + beta

  call yArr%Get(y)

  call mesh%RunCellKernel(Fill_Kernel)

contains
  subroutine Fill_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = gamma
    
  end subroutine Fill_Kernel
end subroutine Fill$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components randomly: 𝒚 ← {𝛼ᵢ}ᵀ, where: 
!! • 𝛼ᵢ ~ 𝘜(𝑎,𝑏), 𝒚 ∊ ℝⁿ,
!! • 𝛼ᵢ = 𝜙ᵢ + 𝑐⋅𝘦𝘹𝘱(𝑖𝜓ᵢ), 𝜙ᵢ ~ 𝘜(𝑎,𝑏), 𝜓ᵢ ~ 𝘜(0,2𝜋), 𝒚 ∊ ℂⁿ.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Fill_Random$T(mesh, yArr, a, b)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(inout) :: yArr
  real(dp), intent(in), optional :: a, b

  $typename, pointer :: y(:,:)

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

#$if T == 'R'
    do iCell = firstCell, lastCell
      call random_number(y(:,iCell))
      if (present(a).and.present(b)) then
        y(:,iCell) = min(a, b) + abs(b - a)*y(:,iCell)
      end if
    end do
#$else
    error stop 'complex Fill_Random is not implemented yet!', "$type" 
#$end if
    
  end subroutine Fill_Random_Kernel
end subroutine Fill_Random$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒚 ← 𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Set$T(mesh, yArr, xArr)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr
  class(tArray$T), intent(inout) :: yArr
  
  $typename, pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(Set_Kernel)

contains
  subroutine Set_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = x(:,iCell)

  end subroutine Set_Kernel
end subroutine Set$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒚 ← 𝕹𝖊(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Set_Real(mesh, y, x)
  class(tMesh), intent(inout) :: mesh
  complex(dp), intent(in) :: x(:,:)
  real(dp), intent(inout) :: y(:,:)
  
  call mesh%RunCellKernel(Set_Real_Kernel)

contains
  subroutine Set_Real_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = Re(x(:,iCell))

  end subroutine Set_Real_Kernel
end subroutine Set_Real

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒚 ← 𝕹𝖊(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Set_Imag(mesh, y, x)
  class(tMesh), intent(inout) :: mesh
  complex(dp), intent(in) :: x(:,:)
  real(dp), intent(inout) :: y(:,:)
  
  call mesh%RunCellKernel(Set_Imag_Kernel)

contains
  subroutine Set_Imag_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = Im(x(:,iCell))

  end subroutine Set_Imag_Kernel
end subroutine Set_Imag

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒛 ← 𝒚 + 𝑖𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Set_Complex(mesh, z, y, x)
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: x(:,:), y(:,:)
  complex(dp), intent(inout) :: z(:,:)
  
  call mesh%RunCellKernel(Set_Complex_Kernel)

contains
  subroutine Set_Complex_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,iCell) = R2C(y(:,iCell), x(:,iCell))

  end subroutine Set_Complex_Kernel
end subroutine Set_Complex

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Scale: 𝒚 ← 𝛼𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Scale$T(mesh, yArr, xArr, alpha)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr
  class(tArray$T), intent(inout) :: yArr
  $typename, intent(in) :: alpha

  $typename, pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(Scale_Kernel)

contains
  subroutine Scale_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = alpha*x(:,iCell)

  end subroutine Scale_Kernel
end subroutine Scale$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: 𝒛 ← [[𝛽]]𝒚 + [𝛼]𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Add$T(mesh, zArr, yArr, xArr, alpha, beta)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr, yArr
  class(tArray$T), intent(inout) :: zArr
  $typename, intent(in), optional :: alpha, beta

  $typename, pointer :: x(:,:), y(:,:), z(:,:)
  $typename :: a, b

  call xArr%Get(x); call yArr%Get(y); call zArr%Get(z)
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present( beta)) b =  beta

  call mesh%RunCellKernel(Add_Kernel)

contains
  subroutine Add_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,iCell) = b*y(:,iCell) + a*x(:,iCell)

  end subroutine Add_Kernel
end subroutine Add$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: 𝒛 ← 𝛽𝒚 - 𝛼𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Sub$T(mesh, zArr, yArr, xArr, alpha, beta)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr, yArr
  class(tArray$T), intent(inout) :: zArr
  $typename, intent(in), optional :: alpha, beta

  $typename, pointer :: x(:,:), y(:,:), z(:,:)
  $typename :: a, b

  call xArr%Get(x); call yArr%Get(y); call zArr%Get(z)
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present( beta)) b =  beta

  call mesh%RunCellKernel(Sub_Kernel)

contains
  subroutine Sub_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,iCell) = b*y(:,iCell) - a*x(:,iCell)
    
  end subroutine Sub_Kernel
end subroutine Sub$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute product: u̅ ← vw̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-3
subroutine Mul$rank(mesh, u, v, w, power)
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: v(:), w(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  integer(ip), intent(in), optional :: power
  
  integer(ip) :: p
  p = 1; if (present(power)) p = power
  
  call mesh%RunCellKernel(Mul_Kernel)

contains
  subroutine Mul_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    u(@:,iCell) = (v(iCell)**p)*w(@:,iCell)

  end subroutine Mul_Kernel
end subroutine Mul$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute an inner product: 𝒛 ← 𝒚⋅𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-3
subroutine Mul_Inner$rank(mesh, z, y, x)
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: y(:,:), x(:,@:,:)
  real(dp), intent(inout) :: z(@:,:)
  
  call mesh%RunCellKernel(Mul_Inner_Kernel)

contains
  subroutine Mul_Inner_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(@:,iCell) = y(:,iCell).inner.x(:,@:,iCell)

  end subroutine Mul_Inner_Kernel
end subroutine Mul_Inner$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute an outer product: 𝒛 ← 𝒚⊗𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-3
subroutine Mul_Outer$rank(mesh, z, y, x)
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: y(:,:), x(@:,:)
  real(dp), intent(inout) :: z(:,@:,:)
  
  call mesh%RunCellKernel(Mul_Outer_Kernel)

contains
  subroutine Mul_Outer_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    z(:,@:,iCell) = y(:,iCell).outer.x(@:,iCell)

  end subroutine Mul_Outer_Kernel
end subroutine Mul_Outer$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute integral average: 𝑖 ← ∫𝑓(𝒙(𝒓))𝑑𝛺/∫1𝑑𝛺.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
real(dp) function Integrate$T(mesh, xArr, f) result(integral)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr
  procedure(tMapFunc$T) :: f
  
  $typename, pointer :: x(:,:)

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
end function Integrate$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: 𝒚 ← 𝑓(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine FuncProd$T(mesh, yArr, xArr, f)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr
  class(tArray$T), intent(inout) :: yArr
  procedure(tMapFunc$T) :: f

  $typename, pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(FuncProd_Kernel)
  
contains
  subroutine FuncProd_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = f(x(:,iCell))

  end subroutine FuncProd_Kernel
end subroutine FuncProd$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: 𝒚 ← 𝑓(𝒓,𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine SFuncProd$T(mesh, yArr, xArr, f)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr
  class(tArray$T), intent(inout) :: yArr
  procedure(tSMapFunc$T) :: f

  $typename, pointer :: x(:,:), y(:,:)

  call xArr%Get(x); call yArr%Get(y)

  call mesh%RunCellKernel(SFuncProd_Kernel)
  
contains
  subroutine SFuncProd_Kernel(iCell)
    integer(ip), intent(in) :: iCell

    y(:,iCell) = f(mesh%CellCenter(iCell), x(:,iCell))

  end subroutine SFuncProd_Kernel
end subroutine SFuncProd$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Multiply a vector by diagonal of the matrix: 𝓓𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine MatVecProd_Diagonal$T(mesh, DxArr, xArr, MatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr
  class(tArray$T), intent(inout) :: DxArr
  procedure(tMatVecFunc$T) :: MatVec

  $typename, pointer :: x(:,:), Dx(:,:)

  call xArr%Get(x); call DxArr%Get(Dx)

  call mesh%RunCellKernel_Block(MatVecProd_Diagonal_BlockKernel)
  call mesh%SetRange()

contains
  subroutine MatVecProd_Diagonal_BlockKernel(mesh, firstCell, lastCell)
    class(tMesh), intent(inout), target :: mesh
    integer(ip), intent(in) :: firstCell, lastCell

    integer(ip) :: iCell
    type(tArray$T) :: eArr
    $typename, pointer :: e(:,:)

    call eArr%AllocMold(xArr)
    call eArr%Get(e)

    e(:,:) = 0.0_dp

    do iCell = firstCell, lastCell
      e(:,iCell) = x(:,iCell)
      call mesh%SetRange(iCell)
      call MatVec(mesh, DxArr, eArr)
      e(:,iCell) = 0.0_dp
    end do

  end subroutine MatVecProd_Diagonal_BlockKernel  
end subroutine MatVecProd_Diagonal$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Multiply a vector by a lower/upper triangular 
!! part of the matrix: 𝓣𝒙 ← 𝘵𝘳𝘪𝘶(𝓐)𝒙 or 𝓣𝒙 ← 𝘵𝘳𝘪𝘭(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine MatVecProd_Triangular$T(mesh, upLo, TxArr, xArr, MatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: xArr
  class(tArray$T), intent(inout) :: TxArr
  procedure(tMatVecFunc$T) :: MatVec
  character, intent(in) :: upLo

  $typename, pointer :: x(:,:), Tx(:,:)

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
    type(tArray$T) :: eArr
    $typename, pointer :: e(:,:)

    call eArr%AllocMold(xArr)
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
    type(tArray$T) :: eArr
    $typename, pointer :: e(:,:)

    call eArr%AllocMold(xArr)
    call eArr%Get(e)

    e(:,:lastCell) = x(:,:lastCell)
    e(:,lastCell+1:) = 0.0_dp
    
    do iCell = lastCell, firstCell, -1
      call mesh%SetRange(iCell)
      call MatVec(mesh, TxArr, eArr)
      e(:,iCell) = 0.0_dp
    end do

  end subroutine MatVecProd_LowerTriangular_BlockKernel 
end subroutine MatVecProd_Triangular$T
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve a linear system with a lower/upper 
!! triangular part of the matrix: 𝘵𝘳𝘪𝘶(𝓐)𝒙 = 𝒃 or 𝘵𝘳𝘪𝘭(𝓐)𝒙 = 𝒃.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for T, typename in [SCALAR_TYPES[0]]
subroutine Solve_Triangular$T(mesh, upLo, xArr, bArr, diagArr, MatVec)
  class(tMesh), intent(inout) :: mesh
  class(tArray$T), intent(in) :: bArr, diagArr
  class(tArray$T), intent(inout) :: xArr
  procedure(tMatVecFunc$T) :: MatVec
  character, intent(in) :: upLo

  $typename, pointer :: x(:,:), b(:,:), diag(:,:)

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
    type(tArray$T) :: AxArr
    $typename, pointer :: Ax(:,:)

    call AxArr%AllocMold(xArr)
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
    type(tArray$T) :: AxArr
    $typename, pointer :: Ax(:,:)

    call AxArr%AllocMold(xArr)
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
end subroutine Solve_Triangular$T
#$end for

end module StormRuler_BLAS
