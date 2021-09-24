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
module StormRuler_Tridiag_LAPACK

#$use 'StormRuler_Params.fi'
#$if HAS_LAPACK

use StormRuler_Parameters, only: dp, ip, &
  & error_code, not_implemented_code
use StormRuler_Tridiag, only: tTridiagMatrix

use, intrinsic :: iso_fortran_env, only: error_unit

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

#$if HAS_MKL
include 'mkl_lapack.fi'
#$else
external dstebz, dstein, dsterf, dsteqr
#$end if

interface ComputeEigenpairs_Symm_LAPACK
  module procedure ComputeEigenpairs_Symm_LAPACK
  module procedure ComputeEigenpairs_All_Symm_LAPACK
end interface ComputeEigenpairs_Symm_LAPACK

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute the 𝑚 eigenpairs (𝜃ⱼ,𝒚ⱼ) of the symmetric tridiagonal
!! matrix: 𝜽 = {𝜃ⱼ}, 𝓣𝒚ⱼ = 𝜃ⱼ𝒚ⱼ, using bisection for the eigenvalue
!! computation and inverse iterations for the eigenvectors (from LAPACK).
!!
!! Which eigenvalues are computed:
!! • `SM` (smallest magnitude) case: 
!!   ( 𝑗 ∊ [1,𝑚], 𝜃₁ ≤ … ≤ 𝜃ₘ ≤ … ),
!! • `LM` (largest magnitude) case: 
!!   ( 𝑗 ∊ [𝑛-𝑚,𝑛], … ≤ 𝜃ₙ₋ₘ ≤ … ≤ 𝜃ₙ ),
!! • `BE` (both ends, 𝑚 = 2𝑝):
!!   ( 𝑗 ∊ [1,𝑝]∪[𝑛-𝑝,𝑛], 𝜃₁ ≤ … ≤ 𝜃ₚ ≤ … ≤ 𝜃ₙ₋ₚ ≤ … ≤ 𝜃ₙ ).
!!
!! See: `dstebz`, https://intel.ly/3z34f2f
!!      `dstein`, https://intel.ly/3tvtqJE
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ComputeEigenpairs_Symm_LAPACK(T, m, which, theta, y)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tTridiagMatrix), intent(in) :: T
  integer(ip), intent(in) :: m
  character(len=*), intent(in) :: which
  real(dp), intent(inout) :: theta(:)
  real(dp), intent(inout), optional, target :: y(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  character :: range, order
  integer(ip) :: pass, info, m_out, j_min, j_max, nsplit
  integer(ip), pointer :: iwork(:), iblock(:), isplit(:), ifailv(:)
  real(dp), pointer :: work(:), y_out(:,:), theta_out(:)

  ! ----------------------
  ! Prepair the input arguments.
  ! ----------------------
  allocate(iwork(3*T%Dim), &
    & iblock(T%Dim), isplit(T%Dim), theta_out(T%Dim))
  if (present(y)) then
    order = 'B'
    allocate(work(5*T%Dim), ifailv(m))
  else
    order = 'E'
    allocate(work(4*T%Dim))
  end if

  do pass = 1, merge(2, 1, which == 'BE')
    ! ----------------------
    ! Select pass eigenpair ranges.
    ! ----------------------
    select case(which)
      case('SM')
        j_min = 1; j_max = m
      case('LM')
        j_min = T%Dim-m+1; j_max = T%Dim
      case('BE')
        if (pass == 1) then
          j_min = 1; j_max = m/2
        else
          j_min = T%Dim-m/2+1; j_max = T%Dim
        end if
    end select

    ! ----------------------
    ! Compute {𝜃ⱼ} with LAPACK.
    ! ----------------------
    call dstebz('I', order, T%Dim, 0.0_dp, 0.0_dp, j_min, j_max, 0.0_dp, &
      & T%Diag, T%Subdiag(2:), m_out, nsplit, theta_out, iblock, isplit, &
      & work, iwork, info)
    if (info /= 0) then
      write(error_unit, *) 'LAPACK dstebz failed, INFO=', info
      error stop error_code
    end if
    if (which /= 'BE') then
      theta(1:m) = theta_out(1:m_out)
    else
      theta((pass-1)*m/2+1:pass*m/2) = theta_out(1:m_out)
    end if

    ! ----------------------
    ! Compute {𝒚ⱼ} with LAPACK.
    ! ----------------------
    if (present(y)) then
      y_out => y(:, (pass-1)*m/2+1:)
      call dstein(T%Dim, T%Diag, T%Subdiag(2:), m_out, theta_out, &
        & iblock, isplit, y_out, T%Dim, work, iwork, ifailv, info)
      if (info /= 0) then
        write(error_unit, *) &
          & 'LAPACK dstein failed, INFO=', info, ' IFAILV=', ifailv
        error stop error_code
      end if
    end if

    ! ----------------------
    ! Merge-sort (𝜃ⱼ,𝒚ⱼ) from blocks.
    ! ----------------------
    if (nsplit /= 1) then
      ! TODO: implement me.
      write(error_unit, *) 'LAPACK dstebz NSPLIT/=1, =', nsplit
      write(error_unit, *) 'MERGE SORT NOT IMPLEMENTED'
      error stop not_implemented_code
    end if
  end do

end subroutine ComputeEigenpairs_Symm_LAPACK

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute all eigenpairs (𝜃ⱼ,𝒚ⱼ) of the symmetric tridiagonal matrix:
!! 𝜽 = {𝜃ⱼ}, 𝒀 = {𝒚ⱼ}, 𝓣𝒚ⱼ = 𝜃ⱼ𝒚ⱼ, using the QR algorithm (from LAPACK).
!!
!! See: `dsterf`, https://intel.ly/38Wdyqb
!!      `dsteqr`, https://intel.ly/37N95pe
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ComputeEigenpairs_All_Symm_LAPACK(T, theta, y)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tTridiagMatrix), intent(in) :: T
  real(dp), intent(inout) :: theta(:)
  real(dp), intent(inout), optional :: y(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  character :: compute_y
  integer(ip) :: dim_y, info
  real(dp), allocatable :: subdiag(:), work(:)

  ! ----------------------
  ! Prepair the input arguments.
  ! ----------------------
  theta = T%Diag(:)
  subdiag = T%Subdiag(2:)
  if (present(y)) then
    compute_y = 'I'
    allocate(work(2*T%Dim-2))
  end if

  ! ----------------------
  ! Compute {(𝜃ⱼ,𝒚ⱼ)} with LAPACK.
  ! ----------------------
  if (present(y)) then
    ! Use QR iterations.
    call dsteqr(compute_y, T%Dim, theta, subdiag, y, T%Dim, work, info)
    if (info /= 0) then
      write(error_unit, *) 'LAPACK dsteqr failed, INFO=', info
      error stop error_code
    end if
  else
    ! Use root-free QR iterations.
    call dsterf(T%Dim, theta, subdiag, info)
    if (info /= 0) then
      write(error_unit, *) 'LAPACK dsterf failed, INFO=', info
      error stop error_code
    end if
  end if

end subroutine ComputeEigenpairs_All_Symm_LAPACK

#$end if

end module StormRuler_Tridiag_LAPACK
