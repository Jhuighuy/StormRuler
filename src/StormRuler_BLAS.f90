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
use StormRuler_Helpers, only: &
  & @{tMapFunc$$, tSMapFunc$$@|@0, NUM_RANKS}@, &
  & operator(.inner.), operator(.outer.)
use StormRuler_Mesh, only: tMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Dot
#$do rank = 0, NUM_RANKS
  module procedure Dot$rank
#$end do
end interface Dot

interface Norm1
#$do rank = 0, NUM_RANKS
  module procedure Norm1$rank
#$end do
end interface Norm1

interface Norm2
#$do rank = 0, NUM_RANKS
  module procedure Norm2$rank
#$end do
end interface Norm2

interface NormC
#$do rank = 0, NUM_RANKS
  module procedure NormC$rank
#$end do
end interface NormC

interface Fill
#$do rank = 0, NUM_RANKS
  module procedure Fill$rank
#$end do
end interface Fill

interface Set
#$do rank = 0, NUM_RANKS
  module procedure Set$rank
#$end do
end interface Set

interface Add
#$do rank = 0, NUM_RANKS
  module procedure Add$rank
#$end do
end interface Add

interface Sub
#$do rank = 0, NUM_RANKS
  module procedure Sub$rank
#$end do
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

interface FuncProd
#$do rank = 0, NUM_RANKS
  module procedure FuncProd$rank
#$end do 
end interface FuncProd

interface SFuncProd
#$do rank = 0, NUM_RANKS
  module procedure SFuncProd$rank
#$end do 
end interface SFuncProd

abstract interface
#$do rank = 0, NUM_RANKS
  subroutine tMatVecFunc$rank(mesh, Au, u, env)
    import :: dp, tMesh
    ! <<<<<<<<<<<<<<<<<<<<<<
    class(tMesh), intent(in) :: mesh
    real(dp), intent(in), target :: u(@:,:)
    real(dp), intent(inout), target :: Au(@:,:)
    class(*), intent(inout) :: env
    ! >>>>>>>>>>>>>>>>>>>>>>
  end subroutine tMatVecFunc$rank
#$end do
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute dot product: d ← <x⋅y>.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
real(dp) function Dot$rank(mesh, x, y)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: x(@:,:), y(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  Dot$rank = mesh%RunCellKernel_Sum(Dot_Kernel)

contains
  real(dp) function Dot_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

#$if rank == 0
    Dot_Kernel = x(iCell) * y(iCell)
#$else
    Dot_Kernel = sum(x(@:,iCell) * y(@:,iCell))
#$end if
    
  end function Dot_Kernel
end function Dot$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute L₁-norm: d ← ‖x‖₁.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
real(dp) function Norm1$rank(mesh, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: x(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  Norm1$rank = mesh%RunCellKernel_Sum(Norm1_Kernel)

contains
  real(dp) function Norm1_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

#$if rank == 0
    Norm1_Kernel = abs(x(iCell))
#$else
    Norm1_Kernel = sum(abs(x(@:,iCell)))
#$end if
    
  end function Norm1_Kernel
end function Norm1$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute L₂-norm: d ← ‖x‖₂.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
real(dp) function Norm2$rank(mesh, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: x(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  Norm2$rank = sqrt(Dot(mesh, x, x))

end function Norm2$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute L∞-norm: d ← ‖x‖∞.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
real(dp) function NormC$rank(mesh, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: x(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  NormC$rank = mesh%RunCellKernel_Max(NormC_Kernel)

contains
  real(dp) function NormC_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

#$if rank == 0
    NormC_Kernel = abs(x(iCell))
#$else
    NormC_Kernel = maxval(abs(x(@:,iCell)))
#$end if
    
  end function NormC_Kernel
end function NormC$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components: y ← α.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Fill$rank(mesh, y, alpha)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: alpha
  real(dp), intent(inout) :: y(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  call mesh%RunCellKernel(Fill_Kernel)

contains
  subroutine Fill_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(@:,iCell) = alpha
    
  end subroutine Fill_Kernel
end subroutine Fill$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Set: y ← x.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Set$rank(mesh, y, x)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: x(@:,:)
  real(dp), intent(inout) :: y(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  call mesh%RunCellKernel(Set_Kernel)

contains
  subroutine Set_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    y(@:,iCell) = x(@:,iCell)

  end subroutine Set_Kernel
end subroutine Set$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: z ← βy + αx.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Add$rank(mesh, z, y, x, alpha, beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: x(@:,:), y(@:,:)
  real(dp), intent(inout) :: z(@:,:)
  real(dp), intent(in), optional :: alpha, beta
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp) :: a, b
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present(beta)) b = beta

  call mesh%RunCellKernel(Add_Kernel)

contains
  subroutine Add_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    z(@:,iCell) = b*y(@:,iCell) + a*x(@:,iCell)

  end subroutine Add_Kernel
end subroutine Add$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: z ← βy - αx.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Sub$rank(mesh, z, y, x, alpha, beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: x(@:,:), y(@:,:)
  real(dp), intent(inout) :: z(@:,:)
  real(dp), intent(in), optional :: alpha, beta
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  real(dp) :: a, b
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present(beta)) b = beta

  call mesh%RunCellKernel(Sub_Kernel)

contains
  subroutine Sub_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    z(@:,iCell) = b*y(@:,iCell) - a*x(@:,iCell)
    
  end subroutine Sub_Kernel
end subroutine Sub$rank
#$end do

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
!! Compute an inner product: z ← y⋅x.
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
!! Compute an outer product: z ← y⊗x.
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
!! Compute a function product: v ← f(u).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FuncProd$rank(mesh, v, u, f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:)
  real(dp), intent(inout) :: v(@:,:)
  procedure(tMapFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>

  call mesh%RunCellKernel(FuncProd_Kernel)
  
contains
  subroutine FuncProd_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    v(@:,iCell) = f(u(@:,iCell))

  end subroutine FuncProd_Kernel
end subroutine FuncProd$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: v ← f(r,u).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine SFuncProd$rank(mesh, v, u, f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:)
  real(dp), intent(inout) :: v(@:,:)
  procedure(tSMapFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>

  call mesh%RunCellKernel(SFuncProd_Kernel)
  
contains
  subroutine SFuncProd_Kernel(iCell)
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: iCell
    ! >>>>>>>>>>>>>>>>>>>>>>

    v(@:,iCell) = f(mesh%CellCenter(iCell), u(@:,iCell))

  end subroutine SFuncProd_Kernel
end subroutine SFuncProd$rank
#$end do

end module StormRuler_BLAS
