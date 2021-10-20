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
module StormRuler_Sockets

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: c_int8_t, c_short, &
  & c_int, c_long, c_long_long, c_float, c_double, &
  & c_size_t, c_ptr, c_funptr, c_sizeof, c_loc, c_funloc  

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Server callback subroutine.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
abstract interface
  subroutine tServerCallback(socket)
    import :: ip
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer(ip), intent(in) :: socket
    ! >>>>>>>>>>>>>>>>>>>>>>
  end subroutine tServerCallback
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Read from the socket.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SocketRead(socket, data)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: socket
  class(*), intent(inout), contiguous :: data(..)
  ! >>>>>>>>>>>>>>>>>>>>>>

  interface
    subroutine cSocketRead(socket, data, length) bind(C, name='SR_SocketRead')
      import :: c_int, c_size_t, c_ptr
      ! <<<<<<<<<<<<<<<<<<<<<<
      integer(c_int), intent(in), value :: socket
      type(c_ptr), intent(in), value :: data
      integer(c_size_t), intent(in), value :: length
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cSocketRead
  end interface

end subroutine SocketRead

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Write to the socket.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine SocketWrite(socket, data)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: socket
  class(*), intent(in), target :: data
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  interface
    subroutine cSocketWrite(socket, data, length) bind(C, name='SR_SocketWrite')
      import :: c_int, c_size_t, c_ptr
      ! <<<<<<<<<<<<<<<<<<<<<<
      integer(c_int), intent(in), value :: socket
      type(c_ptr), intent(in), value :: data
      integer(c_size_t), intent(in), value :: length
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cSocketWrite
  end interface

  select type(data)
    type is (logical)
      ;
    type is (integer(ip))
      ;
    type is (real(dp))
      ;
    type is (complex(dp))
      ;
  end select

end subroutine SocketWrite

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Run a server.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine RunServer(port, callback)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer(ip), intent(in) :: port
  procedure(tServerCallback) :: callback
  ! >>>>>>>>>>>>>>>>>>>>>>

  interface
    subroutine cRunServer(port, callback) bind(C, name='SR_RunServer')
      import :: c_int, c_funptr
      ! <<<<<<<<<<<<<<<<<<<<<<
      integer(c_int), intent(in), value :: port
      type(c_funptr), intent(in), value :: callback
      ! >>>>>>>>>>>>>>>>>>>>>>
    end subroutine cRunServer
  end interface

  call cRunServer(port, c_funloc(cRunServer))

end subroutine RunServer

end module StormRuler_Sockets
