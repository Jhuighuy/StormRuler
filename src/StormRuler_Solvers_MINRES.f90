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
module StormRuler_Solvers_MINRES

use StormRuler_Consts, only: ip, dp
use StormRuler_Parameters, only: gMaxIterGMRES

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArray, AllocArray

use StormRuler_BLAS, only: Dot, Norm_2, Fill, Set, Scale, Add, Sub
use StormRuler_BLAS, only: tMatVecFunc

use StormRuler_ConvParams, only: tConvParams
use StormRuler_Preconditioner, only: tPreconditioner

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Solve_GMRES
  module procedure Solve_GMRES
end interface Solve_GMRES

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Generate Givens rotation.
!! ----------------------------------------------------------------- !!
subroutine SymOrtho(a, b, cs, sn, rr)
  real(dp), intent(in) :: a, b
  real(dp), intent(out) :: cs, sn, rr

  ! ----------------------
  ! 𝑟𝑟 ← (𝑎² + 𝑏²)¹ᐟ²,
  ! 𝑐𝑠 ← 𝑎/𝑟𝑟, 𝑠𝑛 ← 𝑏/𝑟𝑟. 
  ! ----------------------
  rr = hypot(a, b)
  if (rr > 0.0_dp) then
    cs = a/rr; sn = b/rr
  else
    cs = 1.0_dp; sn = 0.0_dp
  end if

end subroutine SymOrtho

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
!! Solve a linear operator equation: 𝓐[𝓟]𝒙 = 𝒃, 𝒙 = [𝓟]𝒚, using 
!! the monstrous Generalized Minimal Residual method (GMRES).
!! 
!! The classical GMRES(𝑚) implementation with restarts
!! after 𝑚 iterations is used.
!! 
!! GMRES may be applied to the singular problems, and the square
!! least squares problems: ‖(𝓐[𝓟]𝒚 - 𝒃)‖₂ → 𝘮𝘪𝘯, 𝒙 = [𝓟]𝒚, 
!! although convergeance to minimum norm solution is not guaranteed 
!! (is this true?).
!!
!! References:
!! [1] Saad and M.H. Schultz, 
!!     "GMRES: A generalized minimal residual algorithm for solving 
!!      nonsymmetric linear systems", 
!!     SIAM J. Sci. Stat. Comput., 7:856–869, 1986.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
subroutine Solve_GMRES(mesh, xArr, bArr, MatVec, params, pre)
  class(tMesh), intent(in) :: mesh
  class(tArray), intent(in) :: bArr
  class(tArray), intent(inout) :: xArr
  class(tConvParams), intent(inout) :: params
  class(tPreconditioner), intent(inout), optional :: pre
  procedure(tMatVecFunc) :: MatVec

  logical :: converged
  integer(ip) :: i, k
  real(dp) :: chi, phi, phiTilde
  real(dp), pointer :: beta(:), cs(:), sn(:), H(:,:)
  type(tArray) :: QArr, sArr, tArr, rArr

  call AllocArray(rArr, mold=xArr)
  associate(m => gMaxIterGMRES)
    call AllocArray(QArr, shape=[xArr%mShape, m + 1])
    allocate(beta(m + 1), cs(m), sn(m), H(m + 1,m))
  end associate

  if (present(pre)) then
    error stop 'preconditioned GMRES solver is not implemented.'
  end if

  ! ----------------------
  ! Pre-initialize:
  ! 𝜑̃ ← 𝟢,
  ! 𝘤𝘰𝘯𝘷𝘦𝘳𝘨𝘦𝘥 ← 𝙛𝙖𝙡𝙨𝙚.
  ! ----------------------
  phiTilde = 0.0_dp
  converged = .false.

  do while(.not.converged)
    ! ----------------------
    ! Initialize:
    ! 𝒓 ← 𝓐𝒙,
    ! 𝒓 ← 𝒃 - 𝒓,
    ! 𝜑 ← ‖𝒓‖,
    ! 𝗶𝗳 𝜑̃ = 𝟢: // first non-restarted pass.
    !   𝜑̃ ← 𝜑,
    !   Check convergence for 𝜑̃.
    ! 𝗲𝗻𝗱 𝗶𝗳
    ! ----------------------
    call MatVec(mesh, rArr, xArr)
    call Sub(mesh, rArr, bArr, rArr)
    phi = Norm_2(mesh, rArr)
    if (phiTilde <= 0.0) then
      phiTilde = phi
      if (params%Check(phiTilde)) return
    end if

    ! ----------------------
    ! 𝒄𝒔 ← {𝟢}ᵀ, 𝒔𝒏 ← {𝟢}ᵀ,
    ! 𝜷 ← {𝜑,𝟢,…,𝟢}ᵀ,
    ! 𝓠₁ ← 𝒓/𝜑. 
    ! ----------------------
    cs(:) = 0.0_dp; sn(:) = 0.0_dp
    beta(1) = phi; beta(2:) = 0.0_dp
    tArr = QArr%Slice(1); call Scale(mesh, tArr, rArr, 1.0_dp/phi)

    do k = 1, gMaxIterGMRES
      ! ----------------------
      ! Arnoldi iteration:
      ! 𝓠ₖ₊₁ ← 𝓐𝓠ₖ,
      ! 𝗳𝗼𝗿 𝑖 = 𝟣, 𝑘 𝗱𝗼:
      !   𝓗ᵢₖ ← <𝓠ₖ₊₁⋅𝓠ᵢ>,
      !   𝓠ₖ₊₁ ← 𝓠ₖ₊₁ - 𝓗ᵢₖ𝓠ᵢ,
      ! 𝗲𝗻𝗱 𝗳𝗼𝗿
      ! 𝓗ₖ₊₁,ₖ ← ‖𝓠ₖ₊₁‖, 𝓠ₖ₊₁ ← 𝓠ₖ₊₁/𝓗ₖ₊₁,ₖ.  
      ! ----------------------
      sArr = QArr%Slice(k); tArr = QArr%Slice(k+1)
      call MatVec(mesh, tArr, sArr)
      do i = 1, k
        sArr = QArr%Slice(i); H(i,k) = Dot(mesh, tArr, sArr)
        call Sub(mesh, tArr, tArr, sArr, H(i,k))
      end do
      H(k+1,k) = Norm_2(mesh, tArr); call Scale(mesh, tArr, tArr, 1.0_dp/H(k+1,k))

      ! ----------------------
      ! Eliminate the last element in 𝓗
      ! and and update the rotation matrix:
      ! 𝗳𝗼𝗿 𝑖 = 𝟣, 𝑘 - 𝟣 𝗱𝗼:
      !   𝜒 ← 𝒄𝒔ᵢ⋅𝓗ᵢₖ + 𝒔𝒏ᵢ⋅𝓗ᵢ₊₁,ₖ,
      !   𝓗ᵢ₊₁,ₖ ← -𝒔𝒏ᵢ⋅𝓗ᵢₖ + 𝒄𝒔ᵢ⋅𝓗ᵢ₊₁,ₖ 
      !   𝓗ᵢₖ ← 𝜒,
      ! 𝗲𝗻𝗱 𝗳𝗼𝗿
      ! 𝒄𝒔ₖ, 𝒔𝒏ₖ ← 𝘚𝘺𝘮𝘖𝘳𝘵𝘩𝘰(𝓗ₖₖ, 𝓗ₖ₊₁,ₖ),
      ! 𝓗ₖₖ ← 𝒄𝒔ₖ⋅𝓗ₖₖ + 𝒔𝒏ₖ⋅𝓗ₖ₊₁,ₖ,
      ! 𝓗ₖ₊₁,ₖ ← 𝟢.
      ! ----------------------
      do i = 1, k - 1
        chi = cs(i)*H(i,k) + sn(i)*H(i+1,k)
        H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
        H(i,k) = chi
      end do
      call SymOrtho(H(k,k), H(k+1,k), cs(k), sn(k), chi)
      H(k,k) = cs(k)*H(k,k) + sn(k)*H(k+1,k)
      H(k+1,k) = 0.0_dp

      ! ----------------------
      ! Update the residual vector:
      ! 𝜷ₖ₊₁ ← -𝒔𝒏ₖ⋅𝜷ₖ, 𝜷ₖ ← 𝒄𝒔ₖ⋅𝜷ₖ,
      ! 𝜑 ← |𝜷ₖ₊₁|,
      ! Check convergence for 𝜑 and 𝜑/𝜑̃.
      ! ----------------------
      beta(k+1) = -sn(k)*beta(k); beta(k) = cs(k)*beta(k)
      phi = abs(beta(k+1))
      if (params%Check(phi, phi/phiTilde)) exit

    end do

    ! ----------------------
    ! Check if restart is required.
    ! ----------------------
    converged = k <= gMaxIterGMRES
    if (.not.converged) k = gMaxIterGMRES

    ! ----------------------
    ! Compute 𝒙-solution:
    ! 𝜷₁:ₖ ← (𝓗₁:ₖ,₁:ₖ)⁻¹𝜷₁:ₖ, 
    ! 𝗳𝗼𝗿 𝑖 = 1, 𝑘 𝗱𝗼:
    !   𝒙 ← 𝒙 + 𝜷ᵢ𝓠ᵢ.
    ! 𝗲𝗻𝗱 𝗳𝗼𝗿
    ! // Since 𝓗₁:ₖ is upper triangular, 
    ! // operations can be combined:
    ! 𝗳𝗼𝗿 𝑖 = 𝑘, 𝟣, -𝟣 𝗱𝗼:
    !   𝜷ᵢ ← (𝜷ᵢ - <𝓗ᵢ,ᵢ₊₁:ₖ⋅𝜷ᵢ₊₁:ₖ>)/𝓗ᵢᵢ,
    !   𝒙 ← 𝒙 + 𝜷ᵢ𝓠ᵢ.
    ! 𝗲𝗻𝗱 𝗳𝗼𝗿
    ! ----------------------
    do i = k, 1, -1
      beta(i) = beta(i) - dot_product(H(i,(i + 1):k), beta((i + 1):k))
      beta(i) = beta(i)/H(i,i)
      tArr = QArr%Slice(i); call Add(mesh, xArr, xArr, tArr, beta(i))
    end do
  end do

end subroutine Solve_GMRES

end module StormRuler_Solvers_MINRES
