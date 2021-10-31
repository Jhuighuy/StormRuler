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
module StormRuler_FDM_Base_Flux

use StormRuler_Parameters, only: dp

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Reconstruct a value of the flux 𝑓 = 𝒂ₙ𝒖.
!! ----------------------------------------------------------------- !!
elemental function FR_C2(a_i, u_i, a_r, u_r)
  real(dp), intent(in) :: u_i, u_r
  real(dp), intent(in) :: a_i, a_r
  real(dp) :: FR_C2

  FR_C2 = 0.5_dp*(a_i + a_r)*merge(u_i, u_r, (a_i + a_r) > 0.0_dp)
  !FR_C2 = 0.25_dp*(a_i + a_r)*(u_i + u_r)
end function FR_C2

!! ----------------------------------------------------------------- !!
!! Reconstruct a value of the flux 𝑓 = 𝒂ₙ𝒖.
!! ----------------------------------------------------------------- !!
elemental function FR_C4(a_l, u_l, a_i, u_i, a_r, u_r, a_rr, u_rr)
  real(dp), intent(in) :: u_l, u_i, u_r, u_rr
  real(dp), intent(in) :: a_l, a_i, a_r, a_rr
  real(dp) :: FR_C4

  FR_C4 = FR_C2(a_i, u_i, a_r, u_r) ! TODO
end function FR_C4

!! ----------------------------------------------------------------- !!
!! Reconstruct a value of the flux 𝑓 = 𝒂ₙ𝒖.
!! ----------------------------------------------------------------- !!
elemental function FR_C6(a_ll, u_ll, a_l, u_l, a_i, u_i, &
    &                a_r, u_r, a_rr, u_rr, a_rrr, u_rrr)
  real(dp), intent(in) :: u_ll, u_l, u_i, u_r, u_rr, u_rrr
  real(dp), intent(in) :: a_ll, a_l, a_i, a_r, a_rr, a_rrr
  real(dp) :: FR_C6

  FR_C6 = FR_C2(a_i, u_i, a_r, u_r) ! TODO
end function FR_C6

!! ----------------------------------------------------------------- !!
!! Reconstruct a value of the flux 𝑓 = 𝒂ₙ𝒖.
!! ----------------------------------------------------------------- !!
elemental function FR_C8(a_lll, u_lll, a_ll, u_ll, a_l, u_l, a_i, u_i, &
    &              a_r, u_r, a_rr, u_rr, a_rrr, u_rrr, a_rrrr, u_rrrr)
  real(dp), intent(in) :: u_lll, u_ll, u_l, u_i, u_r, u_rr, u_rrr, u_rrrr
  real(dp), intent(in) :: a_lll, a_ll, a_l, a_i, a_r, a_rr, a_rrr, a_rrrr
  real(dp) :: FR_C8

  FR_C8 = FR_C2(a_i, u_i, a_r, u_r) ! TODO
end function FR_C8

end module StormRuler_FDM_Base_Flux
