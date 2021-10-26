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
module StormRuler_Threads

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: ip

use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: c_int, c_long, &
  & c_size_t, c_ptr, c_funptr, c_null_ptr, c_funloc

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

type, bind(C) :: pthread_t
  integer(c_long) :: private = 0_c_long
end type pthread_t

interface
  function pthread_create(thread, pAttr, pStartRoutine, pArg) &
      & bind(C, name='pthread_create')
    import :: c_int, c_ptr, c_funptr, pthread_t
    type(pthread_t), intent(inout) :: thread
    type(c_ptr), intent(in), value :: pAttr
    type(c_funptr), intent(in), value :: pStartRoutine
    type(c_ptr), intent(in), value :: pArg
    integer(c_int) :: pthread_create
  end function pthread_create
end interface

interface
  function pthread_join(thread, pRetVal) bind(C, name='pthread_join')
    import :: c_int, c_ptr, pthread_t
    type(pthread_t), intent(in), value :: thread
    type(c_ptr), intent(in), value :: pRetVal
    integer(c_int) :: pthread_join
  end function pthread_join
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Posix-like thread.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tThread
  type(pthread_t), private :: mObj
contains
  procedure, non_overridable :: Join => tThread_Join
end type tThread

interface tThread
  module procedure :: tThread_New
end interface tThread

abstract interface
  subroutine tStartRoutine()
  end subroutine tStartRoutine
end interface

type, bind(C) :: sem_t
  integer(c_size_t) :: private(4) = 0_c_size_t
end type sem_t

interface
  function sem_init(sem, shared, value) bind(C, name='sem_init')
    import :: c_int, sem_t
    type(sem_t), intent(inout) :: sem
    integer(c_int), intent(in), value :: shared, value
    integer(c_int) :: sem_init
  end function sem_init
end interface

interface
  function sem_destroy(sem) bind(C, name='sem_destroy')
    import :: c_int, sem_t
    type(sem_t), intent(inout) :: sem
    integer(c_int) :: sem_destroy
  end function sem_destroy
end interface

interface
  function sem_post(sem) bind(C, name='sem_post')
    import :: c_int, sem_t
    type(sem_t), intent(inout) :: sem
    integer(c_int) :: sem_post
  end function sem_post
end interface

interface
  function sem_wait(sem) bind(C, name='sem_wait')
    import :: c_int, sem_t
    type(sem_t), intent(inout) :: sem
    integer(c_int) :: sem_wait
  end function sem_wait
end interface

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Posix-like semaphore.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tSem
  type(sem_t), private :: mObj
contains
  final :: tSem_Destroy 
  procedure, non_overridable :: Post => tSem_Post
  procedure, non_overridable :: Wait => tSem_Wait
end type tSem

interface tSem
  module procedure tSem_New
end interface tSem

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Create a thread.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
function tThread_New(StartRoutine) result(thread)
  procedure(tStartRoutine) :: StartRoutine
  type(tThread) :: thread

  procedure(tStartRoutine), pointer, save :: sStartRoutine
  integer(c_int) :: cError

  sStartRoutine => StartRoutine

  cError = pthread_create(thread%mObj, &
    & c_null_ptr, c_funloc(cStartRoutine), c_null_ptr)
  if (cError /= 0_c_int) then
    write(error_unit, *) '`pthread_create` failed, error=', cError
    error stop
  end if
  
contains
  recursive function cStartRoutine(arg) result(exitCode) bind(C)
    type(c_ptr), intent(in), value :: arg
    type(c_ptr) :: exitCode

    print *, 'Entering PThread..'

    call sStartRoutine()
    exitCode = c_null_ptr

    print *, 'Leaving PThread..'

  end function cStartRoutine
end function tThread_New

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Join the thread.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine tThread_Join(thread)
  class(tThread), intent(inout) :: thread

  integer(c_int) :: cError

  cError = pthread_join(thread%mObj, c_null_ptr)
  if (cError /= 0_c_int) then
    write(error_unit, *) '`pthread_join` failed, error=', cError
    error stop
  end if

end subroutine tThread_Join

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Create a semaphore.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
function tSem_New(value, shared) result(sem)
  integer(ip), intent(in) :: value
  logical, intent(in), optional :: shared
  type(tSem) :: sem

  integer(c_int) :: cShared
  integer(c_int) :: cError

  cShared = 0_c_int
  if (present(shared)) then
    cShared = merge(1_c_int, 0_c_int, shared)
  end if

  cError = sem_init(sem%mObj, cShared, int(value, kind=c_int))
  if (cError /= 0_c_int) then
    write(error_unit, *) '`sem_init` failed, error=', cError
    error stop
  end if

end function tSem_New

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Destroy the semaphore.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine tSem_Destroy(sem)
  type(tSem), intent(inout) :: sem

  integer(c_int) :: cError

  cError = sem_destroy(sem%mObj)
  if (cError /= 0_c_int) then
    write(error_unit, *) '`sem_destroy` failed, error=', cError
    error stop
  end if

end subroutine tSem_Destroy

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Post the semaphore.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine tSem_Post(sem)
  class(tSem), intent(inout) :: sem

  integer(c_int) :: cError

  cError = sem_post(sem%mObj)
  if (cError /= 0_c_int) then
    write(error_unit, *) '`sem_post` failed, error=', cError
    error stop
  end if

end subroutine tSem_Post

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Wait for the semaphore.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine tSem_Wait(sem)
  class(tSem), intent(inout) :: sem

  integer(c_int) :: cError

  cError = sem_wait(sem%mObj)
  if (cError /= 0_c_int) then
    write(error_unit, *) '`sem_wait` failed, error=', cError
    error stop
  end if

end subroutine tSem_Wait

end module StormRuler_Threads
