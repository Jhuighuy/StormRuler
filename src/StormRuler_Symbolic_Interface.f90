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
module StormRuler_Symbolic_Interface

use StormRuler_Parameters, only: dp

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Symbolic expressions evaluation driver interface.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type, abstract :: iSymbolicDriver
contains
    ! Handle symbol-symbol algebraic operation.
    procedure(iSymbolicDriver_Symbol_Symbol), &
      & deferred :: On_Symbol_Symbol
    ! Handle symbol-constant algebraic operation.
    procedure(iSymbolicDriver_Symbol_Constant), &
      & deferred :: On_Symbol_Constant
    ! Handle constant-symbol algebraic operation.
    procedure(iSymbolicDriver_Constant_Symbol), &
      & deferred :: On_Constant_Symbol

    ! Handle symbolic algebraic function operation.
    procedure(iSymbolicDriver_Symbolic_Func), &
      & deferred :: On_Symbolic_Func
    ! Handle symbolic algebraic binary function operation.
    procedure(iSymbolicDriver_Symbolic_BinaryFunc), &
      & deferred :: On_Symbolic_BinaryFunc
end type iSymbolicDriver

abstract interface
end interface

abstract interface
  pure function iSymbolicDriver_Symbol_Symbol(self, operator, left, right)
    import :: iSymbolicDriver
    class(iSymbolicDriver), intent(in) :: self
    character(len=*), intent(in) :: operator
    class(*), intent(in) :: left, right
    class(*), allocatable :: iSymbolicDriver_Symbol_Symbol
  end function iSymbolicDriver_Symbol_Symbol

  pure function iSymbolicDriver_Symbol_Constant(self, operator, left, right)
    import :: dp, iSymbolicDriver
    class(iSymbolicDriver), intent(in) :: self
    character(len=*), intent(in) :: operator
    class(*), intent(in) :: left 
    real(dp), intent(in) :: right
    class(*), allocatable :: iSymbolicDriver_Symbol_Constant
  end function iSymbolicDriver_Symbol_Constant
  
  pure function iSymbolicDriver_Constant_Symbol(self, operator, left, right)
    import :: dp, iSymbolicDriver
    class(iSymbolicDriver), intent(in) :: self
    character(len=*), intent(in) :: operator
    real(dp), intent(in) :: left
    class(*), intent(in) :: right 
    class(*), allocatable :: iSymbolicDriver_Constant_Symbol
  end function iSymbolicDriver_Constant_Symbol
end interface

abstract interface
  pure function iSymbolicDriver_Symbolic_Func(self, operator, argument)
    import :: dp, iSymbolicDriver
    class(iSymbolicDriver), intent(in) :: self
    character(len=*), intent(in) :: operator
    class(*), intent(in) :: argument
    class(*), allocatable :: iSymbolicDriver_Symbolic_Func
  end function iSymbolicDriver_Symbolic_Func

  pure function iSymbolicDriver_Symbolic_BinaryFunc(self, operator, first, second)
    import :: dp, iSymbolicDriver
    class(iSymbolicDriver), intent(in) :: self
    character(len=*), intent(in) :: operator
    class(*), intent(in) :: first, second
    class(*), allocatable :: iSymbolicDriver_Symbolic_BinaryFunc
  end function iSymbolicDriver_Symbolic_BinaryFunc
end interface

class(iSymbolicDriver), allocatable :: gSymbolicDriver

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Symbolic expression container class.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tSymbol
  class(*), allocatable :: inner
end type tSymbol

interface operator(+)
  module procedure Symbol_Plus_Symbol
  module procedure Symbol_Plus_Constant
  module procedure Constant_Plus_Symbol
end interface

private :: Symbol_Plus_Symbol, &
  & Symbol_Plus_Constant, Constant_Plus_Symbol

interface operator(-)
  module procedure Symbol_Minus_Symbol
  module procedure Symbol_Minus_Constant
  module procedure Constant_Minus_Symbol
end interface

private :: Symbol_Minus_Symbol, &
  & Symbol_Minus_Constant, Constant_Minus_Symbol

interface operator(*)
  module procedure Symbol_Times_Symbol
  module procedure Symbol_Times_Constant
  module procedure Constant_Times_Symbol
end interface

private :: Symbol_Times_Symbol, &
  & Symbol_Times_Constant, Constant_Times_Symbol

interface operator(/)
  module procedure Symbol_Over_Symbol
  module procedure Symbol_Over_Constant
  module procedure Constant_Over_Symbol
end interface

private :: Symbol_Over_Symbol, &
  & Symbol_Over_Constant, Constant_Over_Symbol

interface operator(**)
  module procedure Symbol_Power_Symbol
  module procedure Symbol_Power_Constant
  module procedure Constant_Power_Symbol
end interface

private :: Symbol_Power_Symbol, &
  & Symbol_Power_Constant, Constant_Power_Symbol

#$let UNARY_FUNCTIONS = [ &
  & 'acos', 'acosd', 'acosh', 'asin', 'asind', 'asinh', &
  & 'atan', 'atand', 'atanh', 'bessel_j0', 'bessel_j1', &
  & 'bessel_y0', 'bessel_y1', 'cos', 'cosd', 'cosh', &
  & 'contan', 'cotand', 'erf', 'erfc', 'exp', 'log', &
  & 'log_gamma', 'sign', 'sin', 'sind', 'sinh', 'sqrt' ]

#$let BINARY_FUNCTIONS = [ &
  & 'atan2', 'atan2d', 'hypot', 'max', 'min' ]

#$for Func in (UNARY_FUNCTIONS + BINARY_FUNCTIONS)
interface $Func
  module procedure Symbolic_$Func
end interface $Func
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! 'symbol' + 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Plus_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left, right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Plus_Symbol = tSymbol( &
    & gSymbolicDriver%On_Symbol_Symbol('+', left%inner, right%inner) )
end function Symbol_Plus_Symbol
!! ----------------------------------------------------------------- !!
!! 'symbol' + 'constant' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Plus_Constant(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left 
  real(dp), intent(in) :: right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Plus_Constant = tSymbol( &
    & gSymbolicDriver%On_Symbol_Constant('+', left%inner, right) )
end function Symbol_Plus_Constant
!! ----------------------------------------------------------------- !!
!! 'constant' + 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Constant_Plus_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: left
  type(tSymbol), intent(in) :: right 
  ! >>>>>>>>>>>>>>>>>>>>>>
  Constant_Plus_Symbol = tSymbol( &
    & gSymbolicDriver%On_Constant_Symbol('+', left, right%inner) )
end function Constant_Plus_Symbol

!! ----------------------------------------------------------------- !!
!! 'symbol' - 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Minus_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left, right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Minus_Symbol = tSymbol( &
    & gSymbolicDriver%On_Symbol_Symbol('-', left%inner, right%inner) )
end function Symbol_Minus_Symbol
!! ----------------------------------------------------------------- !!
!! 'symbol' - 'constant' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Minus_Constant(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left 
  real(dp), intent(in) :: right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Minus_Constant = tSymbol( &
    & gSymbolicDriver%On_Symbol_Constant('-', left%inner, right) )
end function Symbol_Minus_Constant
!! ----------------------------------------------------------------- !!
!! 'constant' - 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Constant_Minus_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: left
  type(tSymbol), intent(in) :: right 
  ! >>>>>>>>>>>>>>>>>>>>>>
  Constant_Minus_Symbol = tSymbol( &
    & gSymbolicDriver%On_Constant_Symbol('-', left, right%inner) )
end function Constant_Minus_Symbol

!! ----------------------------------------------------------------- !!
!! 'symbol' * 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Times_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left, right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Times_Symbol = tSymbol( &
    & gSymbolicDriver%On_Symbol_Symbol('*', left%inner, right%inner) )
end function Symbol_Times_Symbol
!! ----------------------------------------------------------------- !!
!! 'symbol' * 'constant' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Times_Constant(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left 
  real(dp), intent(in) :: right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Times_Constant = tSymbol( &
    & gSymbolicDriver%On_Symbol_Constant('*', left%inner, right) )
end function Symbol_Times_Constant
!! ----------------------------------------------------------------- !!
!! 'constant' * 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Constant_Times_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: left
  type(tSymbol), intent(in) :: right 
  ! >>>>>>>>>>>>>>>>>>>>>>
  Constant_Times_Symbol = tSymbol( &
    & gSymbolicDriver%On_Constant_Symbol('*', left, right%inner) )
end function Constant_Times_Symbol

!! ----------------------------------------------------------------- !!
!! 'symbol' / 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Over_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left, right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Over_Symbol = tSymbol( &
    & gSymbolicDriver%On_Symbol_Symbol('/', left%inner, right%inner) )
end function Symbol_Over_Symbol
!! ----------------------------------------------------------------- !!
!! 'symbol' / 'constant' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Over_Constant(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left 
  real(dp), intent(in) :: right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Over_Constant = tSymbol( &
    & gSymbolicDriver%On_Symbol_Constant('/', left%inner, right) )
end function Symbol_Over_Constant
!! ----------------------------------------------------------------- !!
!! 'constant' / 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Constant_Over_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: left
  type(tSymbol), intent(in) :: right 
  ! >>>>>>>>>>>>>>>>>>>>>>
  Constant_Over_Symbol = tSymbol( &
    & gSymbolicDriver%On_Constant_Symbol('/', left, right%inner) )
end function Constant_Over_Symbol

!! ----------------------------------------------------------------- !!
!! 'symbol' ** 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Power_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left, right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Power_Symbol = tSymbol( &
    & gSymbolicDriver%On_Symbol_Symbol('**', left%inner, right%inner) )
end function Symbol_Power_Symbol
!! ----------------------------------------------------------------- !!
!! 'symbol' ** 'constant' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_Power_Constant(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left 
  real(dp), intent(in) :: right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_Power_Constant = tSymbol( &
    & gSymbolicDriver%On_Symbol_Constant('**', left%inner, right) )
end function Symbol_Power_Constant
!! ----------------------------------------------------------------- !!
!! 'constant' ** 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Constant_Power_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: left
  type(tSymbol), intent(in) :: right 
  ! >>>>>>>>>>>>>>>>>>>>>>
  Constant_Power_Symbol = tSymbol( &
    & gSymbolicDriver%On_Constant_Symbol('**', left, right%inner) )
end function Constant_Power_Symbol

#$for Func in UNARY_FUNCTIONS
!! ----------------------------------------------------------------- !!
!! $Func('symbol') intrinsic overload.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbolic_$Func(argument)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: argument 
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbolic_$Func = tSymbol( &
    & gSymbolicDriver%On_Symbolic_Func('$Func', argument%inner) )
end function Symbolic_$Func
#$end for

#$for BinaryFunc in BINARY_FUNCTIONS
!! ----------------------------------------------------------------- !!
!! $BinaryFunc('symbol', 'symbol') intrinsic overload.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbolic_$BinaryFunc(first, second)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: first, second
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbolic_$BinaryFunc = tSymbol( &
    & gSymbolicDriver%On_Symbolic_BinaryFunc('$BinaryFunc', first%inner, second%inner) )
end function Symbolic_$BinaryFunc
#$end for

end module StormRuler_Symbolic_Interface
