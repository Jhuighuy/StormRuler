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

class(iSymbolicDriver), pointer :: gSymbolicDriver

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Symbolic expression container class.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tSymbol
  class(*), allocatable :: inner
end type tSymbol

interface assignment(=)
  module procedure Symbol_Assign_Symbol
  module procedure Symbol_Assign_Constant
end interface

#$let OPERATORS = [ &
  & ('+', 'Plus'), ('-', 'Minus'), &
  & ('*', 'Times'), ('/', 'Over'), ('**', 'Power'), ]

#$for operator, operatorName in OPERATORS
interface operator($operator)
  module procedure Symbol_${operatorName}$_Symbol
  module procedure Symbol_${operatorName}$_Constant
  module procedure Constant_${operatorName}$_Symbol
end interface
#$end for

private :: Symbol_Power_Symbol, &
  & Symbol_Power_Constant, Constant_Power_Symbol

#$let UNARY_FUNCTIONS = [ &
  & 'acos', 'acosd', 'acosh', 'asin', 'asind', 'asinh', &
  & 'atan', 'atand', 'atanh', 'bessel_j0', 'bessel_j1', &
  & 'bessel_y0', 'bessel_y1', 'cos', 'cosd', 'cosh',    &
  & 'contan', 'cotand', 'erf', 'erfc', 'exp', 'log',    &
  & 'log_gamma', 'sign', 'sin', 'sind', 'sinh', 'sqrt', ]

#$let BINARY_FUNCTIONS = [ &
  & 'atan2', 'atan2d', 'hypot', 'max', 'min', ]

#$for func in (UNARY_FUNCTIONS + BINARY_FUNCTIONS)
interface X$func
  module procedure Symbolic_$func
end interface X$func
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

elemental subroutine Symbol_Assign_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(out) :: left
  type(tSymbol), intent(in) :: right
  ! >>>>>>>>>>>>>>>>>>>>>>
end subroutine Symbol_Assign_Symbol

elemental subroutine Symbol_Assign_Constant(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(out) :: left
  real(dp), intent(in) :: right
  ! >>>>>>>>>>>>>>>>>>>>>>
end subroutine Symbol_Assign_Constant

#$for operator, operatorName in OPERATORS
!! ----------------------------------------------------------------- !!
!! 'symbol' $operator 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_${operatorName}$_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left, right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_${operatorName}$_Symbol = tSymbol( &
    & gSymbolicDriver%On_Symbol_Symbol('$operator', left%inner, right%inner) )
end function Symbol_${operatorName}$_Symbol
!! ----------------------------------------------------------------- !!
!! 'symbol' $operator 'constant' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbol_${operatorName}$_Constant(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: left 
  real(dp), intent(in) :: right
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbol_${operatorName}$_Constant = tSymbol( &
    & gSymbolicDriver%On_Symbol_Constant('$operator', left%inner, right) )
end function Symbol_${operatorName}$_Constant
!! ----------------------------------------------------------------- !!
!! 'constant' $operator 'symbol' operator.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Constant_${operatorName}$_Symbol(left, right)
  ! <<<<<<<<<<<<<<<<<<<<<<
  real(dp), intent(in) :: left
  type(tSymbol), intent(in) :: right 
  ! >>>>>>>>>>>>>>>>>>>>>>
  Constant_${operatorName}$_Symbol = tSymbol( &
    & gSymbolicDriver%On_Constant_Symbol('$operator', left, right%inner) )
end function Constant_${operatorName}$_Symbol
#$end for

#$for func in UNARY_FUNCTIONS
!! ----------------------------------------------------------------- !!
!! $func('symbol') intrinsic overload.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbolic_$func(argument)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: argument 
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbolic_$func = tSymbol( &
    & gSymbolicDriver%On_Symbolic_Func('$func', argument%inner) )
end function Symbolic_$func
#$end for

#$for binaryFunc in BINARY_FUNCTIONS
!! ----------------------------------------------------------------- !!
!! $binaryFunc('symbol', 'symbol') intrinsic overload.
!! ----------------------------------------------------------------- !!
elemental type(tSymbol) function Symbolic_$binaryFunc(first, second)
  ! <<<<<<<<<<<<<<<<<<<<<<
  type(tSymbol), intent(in) :: first, second
  ! >>>>>>>>>>>>>>>>>>>>>>
  Symbolic_$binaryFunc = tSymbol( &
    & gSymbolicDriver%On_Symbolic_BinaryFunc('$binaryFunc', first%inner, second%inner) )
end function Symbolic_$binaryFunc
#$end for

end module StormRuler_Symbolic_Interface
