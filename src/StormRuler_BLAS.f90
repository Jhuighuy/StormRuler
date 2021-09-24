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

use StormRuler_Parameters, only: dp, ip, not_implemented_code
use StormRuler_Helpers, only: Re, Im, R2C, &
  & operator(.inner.), operator(.outer.)
use StormRuler_Mesh, only: tMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Dot
#$for type, _ in SCALAR_TYPES
  module procedure Dot$type
#$end for
end interface Dot

interface Norm_1
#$for type, _ in SCALAR_TYPES
  module procedure Norm_1$type
#$end for
end interface Norm_1

interface Norm_2
#$for type, _ in SCALAR_TYPES
  module procedure Norm_2$type
#$end for
end interface Norm_2

interface Norm_C
#$for type, _ in SCALAR_TYPES
  module procedure Norm_C$type
#$end for
end interface Norm_C

interface Fill
#$for type, _ in SCALAR_TYPES
  module procedure Fill$type
#$end for
end interface Fill

interface Fill_Random
#$for type, _ in SCALAR_TYPES
  module procedure Fill_Random$type
#$end for
end interface Fill_Random

interface Set
#$for type, _ in SCALAR_TYPES
  module procedure Set$type
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
#$for type, _ in SCALAR_TYPES
  module procedure Scale$type
#$end for
end interface Scale

interface Add
#$for type, _ in SCALAR_TYPES
  module procedure Add$type
#$end for
end interface Add

interface Sub
#$for type, _ in SCALAR_TYPES
  module procedure Sub$type
#$end for
end interface Sub

interface Mul
#$do rank = 0, NUM_RANKS
  module procedure Mul$rank
#$end do
end interface Mul

interface Mul_Inner
#$do rank = 0, NUM_RANKS-1
  module procedure Mul_Inner$rank
#$end do
end interface Mul_Inner

interface Mul_Outer
#$do rank = 0, NUM_RANKS-1
  module procedure Mul_Outer$rank
#$end do
end interface Mul_Outer

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Mathematical function: ℳ𝒙 ← 𝑓(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for type, typename in SCALAR_TYPES
  pure function tMapFunc$type(x) result(Mx)
    import dp
    ! <<<<<<<<<<<<<<<<<<<<<<
    $typename, intent(in) :: x(:)
    $typename :: Mx(size(x))
    ! >>>>>>>>>>>>>>>>>>>>>>
  end function tMapFunc$type
#$end for
end interface

interface FuncProd
#$for type, _ in SCALAR_TYPES
  module procedure FuncProd$type
#$end for
end interface FuncProd

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Mathematical function: ℳ𝒙 ← 𝑓(𝒓,𝒙), 𝒓 ∊ 𝛺.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for type, typename in SCALAR_TYPES
  pure function tSMapFunc$type(r, x) result(SMx)
    import dp
    ! <<<<<<<<<<<<<<<<<<<<<<
    real(dp), intent(in) :: r(:)
    $typename, intent(in) :: x(:)
    $typename :: SMx(size(x))
    ! >>>>>>>>>>>>>>>>>>>>>>
  end function tSMapFunc$type
#$end for
end interface

interface SFuncProd
#$for type, _ in SCALAR_TYPES
  module procedure SFuncProd$type
#$end for
end interface SFuncProd

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Matrix-vector product function.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
#$for type, typename in SCALAR_TYPES
  subroutine tMatVecFunc$type(mesh, Au, u, env)
    import :: dp, tMesh
    ! <<<<<<<<<<<<<<<<<<<<<<
    class(tMesh), intent(in) :: mesh
    $typename, intent(in), target :: u(:,:)
    $typename, intent(inout), target :: Au(:,:)
    class(*), intent(inout) :: env
    ! >>>>>>>>>>>>>>>>>>>>>>
  end subroutine tMatVecFunc$type
#$end for
end interface

interface MatVecProd_Diagonal
#$for type, _ in SCALAR_TYPES
  module procedure MatVecProd_Diagonal$type
#$end for
end interface MatVecProd_Diagonal

interface MatVecProd_Triangular
#$for type, _ in SCALAR_TYPES
  module procedure MatVecProd_Triangular$type
#$end for
end interface MatVecProd_Triangular

interface Solve_Triangular
#$for type, _ in SCALAR_TYPES
  module procedure Solve_Triangular$type
#$end for
end interface Solve_Triangular

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute dot product: 
!! • 𝑑 ← <𝒙⋅𝒚> = 𝒙ᴴ𝒚 (default), or 
!! • 𝑑 ← [𝒙⋅𝒚] = 𝒙ᵀ𝒚 (do_conjg = false).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
$typename function Dot$type(mesh, x, y, do_conjg) result(Dot)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:), y(:,:)
  logical, intent(in), optional :: do_conjg
  ! >>>>>>>>>>>>>>>>>>>>>>
  
#$if type == 'c'
  if (present(do_conjg)) then
    if (.not.do_conjg) then
      Dot = mesh%RunCellKernel_Sum(Dot_Kernel)
      return
    end if
  end if
  Dot = mesh%RunCellKernel_Sum(Dot_Kernel_Conjg)
#$else
  Dot = mesh%RunCellKernel_Sum(Dot_Kernel)
#$end if

contains
  $typename function Dot_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    ! ----------------------
    ! 𝗼𝘂𝘁 ← 𝗼𝘂𝘁 + 𝒙ᵢ𝒚ᵢ. 
    ! ----------------------
    Dot_Kernel = sum(x(:,iCell) * y(:,iCell))
    
  end function Dot_Kernel
#$if type == 'c'
  $typename function Dot_Kernel_Conjg(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    ! ----------------------
    ! 𝗼𝘂𝘁 ← 𝗼𝘂𝘁 + 𝒙̅ᵢ𝒚ᵢ.
    ! ----------------------
    Dot_Kernel_Conjg = sum(conjg(x(:,iCell)) * y(:,iCell))
    
  end function Dot_Kernel_Conjg
#$end if
end function Dot$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ₁-norm: 𝑑 ← ‖𝒙‖₁.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
real(dp) function Norm_1$type(mesh, x) result(Norm_1)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  Norm_1 = mesh%RunCellKernel_Sum(Norm_1_Kernel)

contains
  real(dp) function Norm_1_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    Norm_1_Kernel = sum(abs(x(:,iCell)))
    
  end function Norm_1_Kernel
end function Norm_1$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ₂-norm: 𝑑 ← ‖𝒙‖₂.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
real(dp) function Norm_2$type(mesh, x) result(Norm_2)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

#$if type == 'c'
  Norm_2 = sqrt( Re(Dot(mesh, x, x)) )
#$else
  Norm_2 = sqrt( Dot(mesh, x, x) )
#$end if

end function Norm_2$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute ℒ∞-norm: 𝑑 ← ‖𝒙‖∞.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
real(dp) function Norm_C$type(mesh, x) result(Norm_C)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  Norm_C = mesh%RunCellKernel_Max(Norm_C_Kernel)

contains
  real(dp) function Norm_C_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    Norm_C_Kernel = maxval(abs(x(:,iCell)))
    
  end function Norm_C_Kernel
end function Norm_C$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components: 𝒚 ← 𝛼 + [𝛽], 𝛼 ∊ ℝ, [𝛽 ∊ ℝ or ℂ].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine Fill$type(mesh, y, alpha, beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: alpha
  $typename, intent(inout) :: y(:,:)
  $typename, intent(in), optional :: beta
  ! >>>>>>>>>>>>>>>>>>>>>>

  $typename :: gamma

  gamma = alpha
  if (present(beta)) gamma = gamma + beta

  call mesh%RunCellKernel(Fill_Kernel)

contains
  subroutine Fill_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(:,iCell) = gamma
    
  end subroutine Fill_Kernel
end subroutine Fill$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components randomly: 𝒚 ← {𝛼ᵢ}ᵀ, where: 
!! • 𝛼ᵢ ~ 𝘜(𝑎,𝑏), 𝒚 ∊ ℝⁿ,
!! • 𝕹𝖊(𝛼ᵢ) ~ ???, 𝕴𝖒(𝛼ᵢ) ~ ???, 𝒚 ∊ ℂⁿ.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine Fill_Random$type(mesh, y, a, b)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  $typename, intent(inout) :: y(:,:)
  real(dp), intent(in), optional :: a, b
  ! >>>>>>>>>>>>>>>>>>>>>>

  ! TODO: not very parallel..
  call mesh%SetRange(parallel=.false.)
  call mesh%RunCellKernel_Block(Fill_Random_Kernel)
  call mesh%SetRange()

contains
  subroutine Fill_Random_Kernel(firstCell, lastCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: firstCell, lastCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer :: iCell

#$if type == 'r'
    do iCell = firstCell, lastCell
      call random_number(y(@:,iCell))
      if (present(a).and.present(b)) then
        y(:,iCell) = min(a, b) + abs(b - a)*y(:,iCell)
      end if
    end do
#$else
    print *, 'complex Fill_Random is not implemented yet!'
    error stop not_implemented_code
#$end if
    
  end subroutine Fill_Random_Kernel
end subroutine Fill_Random$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒚 ← 𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine Set$type(mesh, y, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:)
  $typename, intent(inout) :: y(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  call mesh%RunCellKernel(Set_Kernel)

contains
  subroutine Set_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(:,iCell) = x(:,iCell)

  end subroutine Set_Kernel
end subroutine Set$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒚 ← 𝕹𝖊(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Set_Real(mesh, y, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  complex(dp), intent(in) :: x(:,:)
  real(dp), intent(inout) :: y(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  call mesh%RunCellKernel(Set_Real_Kernel)

contains
  subroutine Set_Real_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(:,iCell) = Re(x(:,iCell))

  end subroutine Set_Real_Kernel
end subroutine Set_Real

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒚 ← 𝕹𝖊(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Set_Imag(mesh, y, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  complex(dp), intent(in) :: x(:,:)
  real(dp), intent(inout) :: y(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  call mesh%RunCellKernel(Set_Imag_Kernel)

contains
  subroutine Set_Imag_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(:,iCell) = Im(x(:,iCell))

  end subroutine Set_Imag_Kernel
end subroutine Set_Imag

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: 𝒛 ← 𝒚 + 𝑖𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine Set_Complex(mesh, z, y, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: x(:,:), y(:,:)
  complex(dp), intent(inout) :: z(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  call mesh%RunCellKernel(Set_Complex_Kernel)

contains
  subroutine Set_Complex_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    z(:,iCell) = R2C(y(:,iCell), x(:,iCell))

  end subroutine Set_Complex_Kernel
end subroutine Set_Complex

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Scale: 𝒚 ← 𝛼𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine Scale$type(mesh, y, x, alpha)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:), alpha
  $typename, intent(inout) :: y(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  call mesh%RunCellKernel(Scale_Kernel)

contains
  subroutine Scale_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(:,iCell) = alpha*x(:,iCell)

  end subroutine Scale_Kernel
end subroutine Scale$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: 𝒛 ← [[𝛽]]𝒚 + [𝛼]𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine Add$type(mesh, z, y, x, alpha, beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:), y(:,:)
  $typename, intent(inout) :: z(:,:)
  $typename, intent(in), optional :: alpha, beta
  ! >>>>>>>>>>>>>>>>>>>>>>

  $typename :: a, b
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present(beta)) b = beta

  call mesh%RunCellKernel(Add_Kernel)

contains
  subroutine Add_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    z(:,iCell) = b*y(:,iCell) + a*x(:,iCell)

  end subroutine Add_Kernel
end subroutine Add$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: 𝒛 ← 𝛽𝒚 - 𝛼𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine Sub$type(mesh, z, y, x, alpha, beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:), y(:,:)
  $typename, intent(inout) :: z(:,:)
  $typename, intent(in), optional :: alpha, beta
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  $typename :: a, b
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present(beta)) b = beta

  call mesh%RunCellKernel(Sub_Kernel)

contains
  subroutine Sub_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    z(:,iCell) = b*y(:,iCell) - a*x(:,iCell)
    
  end subroutine Sub_Kernel
end subroutine Sub$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute product: u̅ ← vw̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Mul$rank(mesh, u, v, w, power)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: v(:), w(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  integer(ip), intent(in), optional :: power
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: p
  p = 1; if (present(power)) p = power
  
  call mesh%RunCellKernel(Mul_Kernel)

contains
  subroutine Mul_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    u(@:,iCell) = (v(iCell)**p)*w(@:,iCell)

  end subroutine Mul_Kernel
end subroutine Mul$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute an inner product: 𝒛 ← 𝒚⋅𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine Mul_Inner$rank(mesh, z, y, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: y(:,:), x(:,@:,:)
  real(dp), intent(inout) :: z(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  call mesh%RunCellKernel(Mul_Inner_Kernel)

contains
  subroutine Mul_Inner_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    z(@:,iCell) = y(:,iCell).inner.x(:,@:,iCell)

  end subroutine Mul_Inner_Kernel
end subroutine Mul_Inner$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute an outer product: 𝒛 ← 𝒚⊗𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine Mul_Outer$rank(mesh, z, y, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: y(:,:), x(@:,:)
  real(dp), intent(inout) :: z(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  call mesh%RunCellKernel(Mul_Outer_Kernel)

contains
  subroutine Mul_Outer_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    z(:,@:,iCell) = y(:,iCell).outer.x(@:,iCell)

  end subroutine Mul_Outer_Kernel
end subroutine Mul_Outer$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: 𝒚 ← 𝑓(𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine FuncProd$type(mesh, y, x, f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:)
  $typename, intent(inout) :: y(:,:)
  procedure(tMapFunc$type) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>

  call mesh%RunCellKernel(FuncProd_Kernel)
  
contains
  subroutine FuncProd_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(:,iCell) = f(x(:,iCell))

  end subroutine FuncProd_Kernel
end subroutine FuncProd$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: 𝒚 ← 𝑓(𝒓,𝒙).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine SFuncProd$type(mesh, y, x, f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  $typename, intent(in) :: x(:,:)
  $typename, intent(inout) :: y(:,:)
  procedure(tSMapFunc$type) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>

  call mesh%RunCellKernel(SFuncProd_Kernel)
  
contains
  subroutine SFuncProd_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(:,iCell) = f(mesh%CellCenter(iCell), x(:,iCell))

  end subroutine SFuncProd_Kernel
end subroutine SFuncProd$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Multiply a vector by diagonal of the matrix: 𝓓𝒙 ← 𝘥𝘪𝘢𝘨(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine MatVecProd_Diagonal$type(mesh, Dx, x, MatVec, env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  $typename, intent(in) :: x(:,:)
  $typename, intent(inout) :: Dx(:,:)
  procedure(tMatVecFunc$type) :: MatVec
  class(*), intent(inout) :: env
  ! >>>>>>>>>>>>>>>>>>>>>>

  call mesh%RunCellKernel_Block(MatVecProd_Diagonal_BlockKernel)
  call mesh%SetRange()

contains
  subroutine MatVecProd_Diagonal_BlockKernel(firstCell, lastCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: firstCell, lastCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer :: iCell

    $typename, allocatable :: e(:,:)
    allocate(e, mold=x)

    e(:,:) = 0.0_dp

    do iCell = firstCell, lastCell
      e(:,iCell) = x(:,iCell)
      call mesh%SetRange(iCell)
      call MatVec(mesh, Dx, e, env)
      e(:,iCell) = 0.0_dp
    end do

  end subroutine MatVecProd_Diagonal_BlockKernel  
end subroutine MatVecProd_Diagonal$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Multiply a vector by a lower/upper triangular 
!! part of the matrix: 𝓣𝒙 ← 𝘵𝘳𝘪𝘶(𝓐)𝒙 or 𝓣𝒙 ← 𝘵𝘳𝘪𝘭(𝓐)𝒙.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine MatVecProd_Triangular$type(mesh, Tx, x, UpLo, MatVec, env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  $typename, intent(in) :: x(:,:)
  $typename, intent(inout) :: Tx(:,:)
  character :: UpLo
  procedure(tMatVecFunc$type) :: MatVec
  class(*), intent(inout) :: env
  ! >>>>>>>>>>>>>>>>>>>>>>

  if (UpLo == 'U') then
    call mesh%RunCellKernel_Block(MatVecProd_UpperTriangular_BlockKernel)
  else if (UpLo == 'L') then
    call mesh%RunCellKernel_Block(MatVecProd_LowerTriangular_BlockKernel)
  end if
  call mesh%SetRange()

contains
  subroutine MatVecProd_UpperTriangular_BlockKernel(firstCell, lastCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: firstCell, lastCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer :: iCell

    $typename, allocatable :: e(:,:)
    allocate(e, mold=x)

    e(:,:firstCell-1) = 0.0_dp
    e(:,firstCell:) = x(:,firstCell:)

    do iCell = firstCell, lastCell
      call mesh%SetRange(iCell)
      call MatVec(mesh, Tx, e, env)
      e(:,iCell) = 0.0_dp
    end do

  end subroutine MatVecProd_UpperTriangular_BlockKernel 
  subroutine MatVecProd_LowerTriangular_BlockKernel(firstCell, lastCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: firstCell, lastCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer :: iCell

    $typename, allocatable :: e(:,:)
    allocate(e, mold=x)

    e(:,:lastCell) = x(:,:lastCell)
    e(:,lastCell+1:) = 0.0_dp
    
    do iCell = lastCell, firstCell, -1
      call mesh%SetRange(iCell)
      call MatVec(mesh, Tx, e, env)
      e(:,iCell) = 0.0_dp
    end do

  end subroutine MatVecProd_LowerTriangular_BlockKernel 
end subroutine MatVecProd_Triangular$type
#$end for

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve a linear system with a lower/upper 
!! triangular part of the matrix: 𝘵𝘳𝘪𝘶(𝓐)𝒙 = 𝒃 or 𝘵𝘳𝘪𝘭(𝓐)𝒙 = 𝒃.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$for type, typename in SCALAR_TYPES
subroutine Solve_Triangular$type(mesh, x, b, diag, UpLo, MatVec, env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  $typename, intent(in) :: b(:,:), diag(:,:)
  $typename, intent(inout) :: x(:,:)
  character :: UpLo
  procedure(tMatVecFunc$type) :: MatVec
  class(*), intent(inout) :: env
  ! >>>>>>>>>>>>>>>>>>>>>>

  ! TODO: not very parallel..
  call mesh%SetRange(parallel=.false.)
  if (UpLo == 'U') then
    call mesh%RunCellKernel_Block(Solve_UpperTriangular_BlockKernel)
  else if (UpLo == 'L') then
    call mesh%RunCellKernel_Block(Solve_LowerTriangular_BlockKernel)
  end if
  call mesh%SetRange()

contains
  subroutine Solve_UpperTriangular_BlockKernel(firstCell, lastCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: firstCell, lastCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer :: iCell

    $typename, allocatable :: Ax(:,:)
    allocate(Ax, mold=x)

    x(:,:) = 0.0_dp
    Ax(:,:) = 0.0_dp

    x(:,lastCell) = b(:,lastCell)/diag(:,lastCell)
    do iCell = lastCell - 1, firstCell, -1
      call mesh%SetRange(iCell)
      call MatVec(mesh, Ax, x, env)
      ! TODO: this is not a correct diagonal solution in block case!
      x(:,iCell) = (b(:,iCell) - Ax(:,iCell))/diag(:,iCell)
    end do

  end subroutine Solve_UpperTriangular_BlockKernel 
  subroutine Solve_LowerTriangular_BlockKernel(firstCell, lastCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: firstCell, lastCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    integer :: iCell

    $typename, allocatable :: Ax(:,:)
    allocate(Ax, mold=x)

    x(:,:) = 0.0_dp
    Ax(:,:) = 0.0_dp

    x(:,firstCell) = b(:,firstCell)/diag(:,firstCell)
    do iCell = firstCell + 1, lastCell
      call mesh%SetRange(iCell)
      call MatVec(mesh, Ax, x, env)
      ! TODO: this is not a correct diagonal solution in block case!
      x(:,iCell) = (b(:,iCell) - Ax(:,iCell))/diag(:,iCell)
    end do

  end subroutine Solve_LowerTriangular_BlockKernel 
end subroutine Solve_Triangular$type
#$end for

end module StormRuler_BLAS
