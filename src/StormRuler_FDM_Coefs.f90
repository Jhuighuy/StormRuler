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
module StormRuler_FDM_Coefs

#$use 'StormRuler_Parameters.f90'

use StormRuler_Parameters, only: dp
use StormRuler_Symbolic

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

#$do order = 2, 8, 2
interface FD1_C$order
#$for T in ['R', 'S']
  module procedure FD1_C$order$T
#$end for
end interface FD1_C$order
#$end do

#$do order = 1, 8
interface FD1_F$order
#$for T in ['R', 'S']
  module procedure FD1_F$order$T
#$end for
end interface FD1_F$order
#$end do

#$do order = 2, 8, 2
interface FD2_C$order
#$for T in ['R', 'S']
  module procedure FD2_C$order$T
#$end for
end interface FD2_C$order
#$end do

#$do order = 2, 8, 2
interface WFD2_C$order
#$for T in ['R', 'S']
  module procedure WFD2_C$order$T
#$end for
end interface WFD2_C$order
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_C2$T(u_l, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_r, u_l
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C2$T = 0.5_dp*(u_r - u_l)
end function FD1_C2$T
#$end for

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_C4$T(u_ll, u_l, u_r, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_ll, u_l, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C4$T = &
    & ( (-01.0_dp/12.0_dp)*u_rr + &
    &   (+02.0_dp/03.0_dp)*u_r  + &
    &   (-02.0_dp/03.0_dp)*u_l  + &
    &   (+01.0_dp/12.0_dp)*u_ll )
end function FD1_C4$T
#$end for

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_C6$T(u_lll, u_ll, u_l, u_r, u_rr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_lll, u_ll, u_l, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C6$T = &
    & ( (+01.0_dp/60.0_dp)*u_rrr + &
    &   (-03.0_dp/20.0_dp)*u_rr  + &
    &   (+03.0_dp/04.0_dp)*u_r   + &
    &   (-03.0_dp/04.0_dp)*u_l   + &
    &   (+03.0_dp/20.0_dp)*u_ll  + &
    &   (-01.0_dp/60.0_dp)*u_lll )
end function FD1_C6$T
#$end for

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_C8$T(u_llll, u_lll, u_ll, u_l, u_r, u_rr, u_rrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_llll, u_lll, u_ll, u_l, u_r, u_rr, u_rrr, u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_C8$T = &
    & ( (-001.0_dp/280.0_dp)*u_rrrr + &
    &   (+004.0_dp/105.0_dp)*u_rrr  + &
    &   (-001.0_dp/005.0_dp)*u_rr   + &
    &   (+004.0_dp/005.0_dp)*u_r    + &
    &   (-004.0_dp/005.0_dp)*u_l    + &
    &   (+001.0_dp/005.0_dp)*u_ll   + &
    &   (-004.0_dp/105.0_dp)*u_lll  + &
    &   (+001.0_dp/280.0_dp)*u_llll )
end function FD1_C8$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! First order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_F1$T(u, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u, u_r
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F1$T = u_r - u
end function FD1_F1$T
#$end for

!! ----------------------------------------------------------------- !!
!! Second order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_F2$T(u, u_r, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F2$T = &
    & ( (-1.5_dp)*u   + &
    &   (+2.0_dp)*u_r + &
    &   (-0.5_dp)*u_rr )
end function FD1_F2$T
#$end for

!! ----------------------------------------------------------------- !!
!! Third order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_F3$T(u_l, u, u_r, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_l, u, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F3$T = &
    & ( (-1.0_dp/3.0_dp)*u_l + &
    &   (-1.0_dp/2.0_dp)*u   + &
    &                    u_r + &
    &   (-1.0_dp/6.0_dp)*u_rr )
end function FD1_F3$T
#$end for

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_F4$T(u_l, u, u_r, u_rr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_l, u, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F4$T = &
    & ( (-01.0_dp/04.0_dp)*u_l  + &
    &   (-05.0_dp/06.0_dp)*u    + &
    &   (+03.0_dp/02.0_dp)*u_r  + &
    &   (-01.0_dp/02.0_dp)*u_rr + &
    &   (+01.0_dp/12.0_dp)*u_rrr )
end function FD1_F4$T
#$end for

!! ----------------------------------------------------------------- !!
!! Fifth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_F5$T(u_ll, u_l, u, u_r, u_rr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_ll, u_l, u, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F5$T = &
    & ( (+01.0_dp/20.0_dp)*u_ll + &
    &   (-01.0_dp/02.0_dp)*u_l  + &
    &   (-01.0_dp/03.0_dp)*u    + &
    &                      u_r  + &
    &   (-01.0_dp/04.0_dp)*u_rr + &
    &   (+01.0_dp/30.0_dp)*u_rrr )
end function FD1_F5$T
#$end for

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_F6$T(u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: &
    & u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F6$T = &
    & ( (+01.0_dp/30.0_dp)*u_ll  + &
    &   (-02.0_dp/05.0_dp)*u_l   + &
    &   (-07.0_dp/12.0_dp)*u     + &
    &   (+04.0_dp/03.0_dp)*u_r   + &
    &   (-01.0_dp/02.0_dp)*u_rr  + &
    &   (+02.0_dp/15.0_dp)*u_rrr + &
    &   (-01.0_dp/60.0_dp)*u_rrrr )
end function FD1_F6$T
#$end for

!! ----------------------------------------------------------------- !!
!! Seventh order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_F7$T(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F7$T = &
    & ( (-001.0_dp/105.0_dp)*u_lll + &
    &   (+001.0_dp/010.0_dp)*u_ll  + &
    &   (-003.0_dp/005.0_dp)*u_l   + &
    &   (-001.0_dp/004.0_dp)*u     + &
    &                        u_r   + &
    &   (-003.0_dp/010.0_dp)*u_rr  + &
    &   (+001.0_dp/015.0_dp)*u_rrr + &
    &   (-001.0_dp/140.0_dp)*u_rrrr )
end function FD1_F7$T
#$end for

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy forward undivided finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD1_F8$T(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr, u_rrrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr, u_rrrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD1_F8$T = &
    & ( (-001.0_dp/168.0_dp)*u_lll  + &
    &   (+001.0_dp/014.0_dp)*u_ll   + &
    &   (-001.0_dp/002.0_dp)*u_l    + &
    &   (-009.0_dp/020.0_dp)*u      + &
    &   (+005.0_dp/004.0_dp)*u_r    + &
    &   (-001.0_dp/002.0_dp)*u_rr   + &
    &   (+001.0_dp/006.0_dp)*u_rrr  + &
    &   (-001.0_dp/028.0_dp)*u_rrrr + &
    &   (+001.0_dp/280.0_dp)*u_rrrrr )
end function FD1_F8$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD2_C2$T(u_l, u, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_l, u, u_r
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C2$T = u_r - 2.0_dp*u + u_l
end function FD2_C2$T
#$end for

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD2_C4$T(u_ll, u_l, u, u_r, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_ll, u_l, u, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C4$T = &
    & ( (-1.0_dp/12.0_dp)*u_rr + &
    &   (+4.0_dp/03.0_dp)*u_r  + &
    &   (-5.0_dp/02.0_dp)*u    + &
    &   (+4.0_dp/03.0_dp)*u_l  + &
    &   (-1.0_dp/12.0_dp)*u_ll )
end function FD2_C4$T
#$end for

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD2_C6$T(u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C6$T = &
    & ( (+01.0_dp/90.0_dp)*u_rrr + &
    &   (-03.0_dp/20.0_dp)*u_rr  + &
    &   (+03.0_dp/02.0_dp)*u_r   + &
    &   (-49.0_dp/18.0_dp)*u     + &
    &   (+03.0_dp/02.0_dp)*u_l   + &
    &   (-03.0_dp/20.0_dp)*u_ll  + &
    &   (+01.0_dp/90.0_dp)*u_lll )
end function FD2_C6$T
#$end for

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided second finite difference.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function FD2_C8$T(u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  $typename, intent(in) :: u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  FD2_C8$T = &
    & ( (-001.0_dp/560.0_dp)*u_rrrr + &
    &   (+008.0_dp/315.0_dp)*u_rrr  + &
    &   (-001.0_dp/005.0_dp)*u_rr   + &
    &   (+008.0_dp/005.0_dp)*u_r    + &
    &   (-205.0_dp/072.0_dp)*u      + &
    &   (+008.0_dp/005.0_dp)*u_l    + &
    &   (-001.0_dp/005.0_dp)*u_ll   + &
    &   (+008.0_dp/315.0_dp)*u_lll  + &
    &   (-001.0_dp/560.0_dp)*u_llll )
end function FD2_C8$T
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Second order accuracy central undivided d(w⋅du/dx)/dx approximation.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function WFD2_C2$T(w_l, u_l, w, u, w_r, u_r)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: w_l, w, w_r
  $typename, intent(in) :: u_l, u, u_r
  ! >>>>>>>>>>>>>>>>>>>>>>
  WFD2_C2$T = 0.5_dp*( (w_r+w)*(u_r-u) - (w+w_l)*(u-u_l) )
end function WFD2_C2$T
#$end for

!! ----------------------------------------------------------------- !!
!! Fourth order accuracy central undivided d(w⋅du/dx)/dx approximation.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function WFD2_C4$T(w_ll, u_ll, w_l, u_l, w, &
  &                                    u, w_r, u_r, w_rr, u_rr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: w_ll, w_l, w, w_r, w_rr
  $typename, intent(in) :: u_ll, u_l, u, u_r, u_rr
  ! >>>>>>>>>>>>>>>>>>>>>>
  WFD2_C4$T = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO
end function WFD2_C4$T
#$end for

!! ----------------------------------------------------------------- !!
!! Sixth order accuracy central undivided d(w⋅du/dx)/dx approximation.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function WFD2_C6$T(w_lll, u_lll, w_ll, u_ll, w_l, u_l, w, &
  &                                    u, w_r, u_r, w_rr, u_rr, w_rrr, u_rrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: w_lll, w_ll, w_l, w, w_r, w_rr, w_rrr
  $typename, intent(in) :: u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  WFD2_C6$T = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO
end function WFD2_C6$T
#$end for

!! ----------------------------------------------------------------- !!
!! Eighth order accuracy central undivided d(w⋅du/dx)/dx approximation.
!! ----------------------------------------------------------------- !!
#$for T, typename in [('R', 'real(dp)'), ('S', 'type(tSymbol)')]
elemental $typename function WFD2_C8$T(w_llll, u_llll, w_lll, u_lll, w_ll, u_ll, w_l, u_l, w, &
  &                                    u, w_r, u_r, w_rr, u_rr, w_rrr, u_rrr, w_rrrr, u_rrrr)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: w_llll, w_lll, w_ll, w_l, w, w_r, w_rr, w_rrr, w_rrrr
  $typename, intent(in) :: u_llll, u_lll, u_ll, u_l, u, u_r, u_rr, u_rrr, u_rrrr
  ! >>>>>>>>>>>>>>>>>>>>>>
  WFD2_C8$T = WFD2_C2(w_l, u_l, w, u, w_r, u_r) ! TODO
end function WFD2_C8$T
#$end for

end module StormRuler_FDM_Coefs
