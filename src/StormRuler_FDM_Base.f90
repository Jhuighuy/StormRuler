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
module StormRuler_FDM_Base

use StormRuler_Parameters, only: dp

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_C2(u_l, u_r)
  real(dp), intent(in) :: u_r, u_l
  real(dp) :: FD1_C2

  FD1_C2 = 0.5_dp*(u_r - u_l)

end function FD1_C2

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_C4(u_ll, u_l, u_r, u_rr)
  real(dp), intent(in) :: u_ll, u_l, u_r, u_rr
  real(dp) :: FD1_C4

  FD1_C4 = &
    & ( (-01.0_dp/12.0_dp)*u_rr + &
    &   (+02.0_dp/03.0_dp)*u_r  + &
    &   (-02.0_dp/03.0_dp)*u_l  + &
    &   (+01.0_dp/12.0_dp)*u_ll )

end function FD1_C4

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_C6(u_lll, u_ll, u_l, u_r, u_rr, u_rrr)
  real(dp), intent(in) :: u_lll, u_ll, u_l, u_r, u_rr, u_rrr
  real(dp) :: FD1_C6

  FD1_C6 = &
    & ( (+01.0_dp/60.0_dp)*u_rrr + &
    &   (-03.0_dp/20.0_dp)*u_rr  + &
    &   (+03.0_dp/04.0_dp)*u_r   + &
    &   (-03.0_dp/04.0_dp)*u_l   + &
    &   (+03.0_dp/20.0_dp)*u_ll  + &
    &   (-01.0_dp/60.0_dp)*u_lll )

end function FD1_C6

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_C8(u_llll, u_lll, u_ll, u_l, u_r, u_rr, u_rrr, u_rrrr)
  real(dp), intent(in) :: u_llll, u_lll, u_ll, u_l, u_r, u_rr, u_rrr, u_rrrr
  real(dp) :: FD1_C8

  FD1_C8 = &
    & ( (-001.0_dp/280.0_dp)*u_rrrr + &
    &   (+004.0_dp/105.0_dp)*u_rrr  + &
    &   (-001.0_dp/005.0_dp)*u_rr   + &
    &   (+004.0_dp/005.0_dp)*u_r    + &
    &   (-004.0_dp/005.0_dp)*u_l    + &
    &   (+001.0_dp/005.0_dp)*u_ll   + &
    &   (-004.0_dp/105.0_dp)*u_lll  + &
    &   (+001.0_dp/280.0_dp)*u_llll )

end function FD1_C8

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! First order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_F1(u, u_r)
  real(dp), intent(in) :: u, u_r
  real(dp) :: FD1_F1

  FD1_F1 = u_r - u
end function FD1_F1

!! ----------------------------------------------------------------- !!
!! Second order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_F2(u, u_r, u_rr)
  real(dp), intent(in) :: u, u_r, u_rr
  real(dp) :: FD1_F2

  FD1_F2 = &
    & ( (-1.5_dp)*u   + &
    &   (+2.0_dp)*u_r + &
    &   (-0.5_dp)*u_rr )

end function FD1_F2

!! ----------------------------------------------------------------- !!
!! Third order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_F3(u_l, u, u_r, u_rr)
  real(dp), intent(in) :: u_l, u, u_r, u_rr
  real(dp) :: FD1_F3

  FD1_F3 = &
    & ( (-1.0_dp/3.0_dp)*u_l + &
    &   (-1.0_dp/2.0_dp)*u   + &
    &                    u_r + &
    &   (-1.0_dp/6.0_dp)*u_rr )

end function FD1_F3

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_F4(u_l, u, u_r, u_rr, u_rrr)
  real(dp), intent(in) :: u_l, u, u_r, u_rr, u_rrr
  real(dp) :: FD1_F4

  FD1_F4 = &
    & ( (-01.0_dp/04.0_dp)*u_l  + &
    &   (-05.0_dp/06.0_dp)*u    + &
    &   (+03.0_dp/02.0_dp)*u_r  + &
    &   (-01.0_dp/02.0_dp)*u_rr + &
    &   (+01.0_dp/12.0_dp)*u_rrr )

end function FD1_F4

!! ----------------------------------------------------------------- !!
!! Fifth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_F5(u_ll, u_l, u, u_r, u_rr, u_rrr)
  real(dp), intent(in) :: u_ll, u_l, u, u_r, u_rr, u_rrr
  real(dp) :: FD1_F5

  FD1_F5 = &
    & ( (+01.0_dp/20.0_dp)*u_ll + &
    &   (-01.0_dp/02.0_dp)*u_l  + &
    &   (-01.0_dp/03.0_dp)*u    + &
    &                      u_r  + &
    &   (-01.0_dp/04.0_dp)*u_rr + &
    &   (+01.0_dp/30.0_dp)*u_rrr )

end function FD1_F5

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_F6(u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  real(dp), intent(in) :: u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  real(dp) :: FD1_F6

  FD1_F6 = &
    & ( (+01.0_dp/30.0_dp)*u_ll  + &
    &   (-02.0_dp/05.0_dp)*u_l   + &
    &   (-07.0_dp/12.0_dp)*u     + &
    &   (+04.0_dp/03.0_dp)*u_r   + &
    &   (-01.0_dp/02.0_dp)*u_rr  + &
    &   (+02.0_dp/15.0_dp)*u_rrr + &
    &   (-01.0_dp/60.0_dp)*u_rrrr )

end function FD1_F6

!! ----------------------------------------------------------------- !!
!! Seventh order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_F7(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  real(dp), intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  real(dp) :: FD1_F7

  FD1_F7 = &
    & ( (-001.0_dp/105.0_dp)*u_lll + &
    &   (+001.0_dp/010.0_dp)*u_ll  + &
    &   (-003.0_dp/005.0_dp)*u_l   + &
    &   (-001.0_dp/004.0_dp)*u     + &
    &                        u_r   + &
    &   (-003.0_dp/010.0_dp)*u_rr  + &
    &   (+001.0_dp/015.0_dp)*u_rrr + &
    &   (-001.0_dp/140.0_dp)*u_rrrr )

end function FD1_F7

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD1_F8(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr, u_rrrrr)
  real(dp), intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr, u_rrrrr
  real(dp) :: FD1_F8

  FD1_F8 = &
    & ( (-001.0_dp/168.0_dp)*u_lll  + &
    &   (+001.0_dp/014.0_dp)*u_ll   + &
    &   (-001.0_dp/002.0_dp)*u_l    + &
    &   (-009.0_dp/020.0_dp)*u      + &
    &   (+005.0_dp/004.0_dp)*u_r    + &
    &   (-001.0_dp/002.0_dp)*u_rr   + &
    &   (+001.0_dp/006.0_dp)*u_rrr  + &
    &   (-001.0_dp/028.0_dp)*u_rrrr + &
    &   (+001.0_dp/280.0_dp)*u_rrrrr )

end function FD1_F8

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD2_C2(u_l, u, u_r)
  real(dp), intent(in) :: u_l, u, u_r
  real(dp) :: FD2_C2

  FD2_C2 = u_r - 2.0_dp*u + u_l

end function FD2_C2

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD2_C4(u_ll, u_l, u, u_r, u_rr)
  real(dp), intent(in) :: u_ll, u_l, u, u_r, u_rr
  real(dp) :: FD2_C4

  FD2_C4 = &
    & ( (-1.0_dp/12.0_dp)*u_rr + &
    &   (+4.0_dp/03.0_dp)*u_r  + &
    &   (-5.0_dp/02.0_dp)*u    + &
    &   (+4.0_dp/03.0_dp)*u_l  + &
    &   (-1.0_dp/12.0_dp)*u_ll )

end function FD2_C4

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD2_C6(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr)
  real(dp), intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr
  real(dp) :: FD2_C6

  FD2_C6 = &
    & ( (+01.0_dp/90.0_dp)*u_rrr + &
    &   (-03.0_dp/20.0_dp)*u_rr  + &
    &   (+03.0_dp/02.0_dp)*u_r   + &
    &   (-49.0_dp/18.0_dp)*u     + &
    &   (+03.0_dp/02.0_dp)*u_l   + &
    &   (-03.0_dp/20.0_dp)*u_ll  + &
    &   (+01.0_dp/90.0_dp)*u_lll )

end function FD2_C6

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
elemental function FD2_C8(u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  real(dp), intent(in) :: u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  real(dp) :: FD2_C8

  FD2_C8 = &
    & ( (-001.0_dp/560.0_dp)*u_rrrr + &
    &   (+008.0_dp/315.0_dp)*u_rrr  + &
    &   (-001.0_dp/005.0_dp)*u_rr   + &
    &   (+008.0_dp/005.0_dp)*u_r    + &
    &   (-205.0_dp/072.0_dp)*u      + &
    &   (+008.0_dp/005.0_dp)*u_l    + &
    &   (-001.0_dp/005.0_dp)*u_ll   + &
    &   (+008.0_dp/315.0_dp)*u_lll  + &
    &   (-001.0_dp/560.0_dp)*u_llll )

end function FD2_C8

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided 𝑑(𝒘𝑑𝒖/𝑑𝜉)/𝑑𝜉 approximation.
!! ----------------------------------------------------------------- !!
elemental function WFD2_C2(w_l, u_l, w, u, w_r, u_r)
  real(dp), intent(in) :: w_l, w, w_r
  real(dp), intent(in) :: u_l, u, u_r
  real(dp) :: WFD2_C2

  WFD2_C2 = 0.5_dp*( (w_r + w)*(u_r - u) - (w + w_l)*(u - u_l) )

end function WFD2_C2

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided 𝑑(𝒘𝑑𝒖/𝑑𝜉)/𝑑𝜉 approximation.
!! ----------------------------------------------------------------- !!
elemental function WFD2_C4(w_ll, u_ll, w_l, u_l, w, &
    &                      u, w_r, u_r, w_rr, u_rr)
  real(dp), intent(in) :: w_ll, w_l, w, w_r, w_rr
  real(dp), intent(in) :: u_ll, u_l, u, u_r, u_rr
  real(dp) :: WFD2_C4

  WFD2_C4 = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO

end function WFD2_C4

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided 𝑑(𝒘𝑑𝒖/𝑑𝜉)/𝑑𝜉 approximation.
!! ----------------------------------------------------------------- !!
elemental function WFD2_C6(w_lll, u_lll, w_ll, u_ll, w_l, u_l, w, &
    &                      u, w_r, u_r, w_rr, u_rr, w_rrr, u_rrr)
  real(dp), intent(in) :: w_lll, w_ll, w_l, w, w_r, w_rr, w_rrr
  real(dp), intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr
  real(dp) :: WFD2_C6

  WFD2_C6 = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO

end function WFD2_C6

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided 𝑑(𝒘𝑑𝒖/𝑑𝜉)/𝑑𝜉 approximation.
!! ----------------------------------------------------------------- !!
elemental function WFD2_C8(w_llll, u_llll, w_lll, u_lll, w_ll, u_ll, w_l, u_l, w, &
    &                      u, w_r, u_r, w_rr, u_rr, w_rrr, u_rrr, w_rrrr, u_rrrr)
  real(dp), intent(in) :: w_llll, w_lll, w_ll, w_l, w, w_r, w_rr, w_rrr, w_rrrr
  real(dp), intent(in) :: u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  real(dp) :: WFD2_C8

  WFD2_C8 = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO

end function WFD2_C8

end module StormRuler_FDM_Base
