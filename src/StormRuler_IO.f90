!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module StormRuler_IO

use StormRuler_Helpers
use StormRuler_Mesh

implicit none

integer, parameter :: IO_NAMELEN=10

type :: IOListNode
  character(len=IO_NAMELEN) :: Name
  real(dp), pointer :: Data(:), DataVector(:,:), DataTensor(:,:,:)
  type(IOListNode), pointer :: NextNode => null()
end type IOListNode

type :: IODataSet
  type(IOListNode), pointer :: FirstNode=>null()&
                            ,CurrentNode=>null()&
                               ,LastNode=>null() 
contains
  generic, public :: AddNode=>AddNodeScalar&
                             ,AddNodeVector&
                             ,AddNodeTensor
  procedure, private :: AllocateNode=>IODataSet_AllocateNode
  procedure, private :: AddNodeScalar=>IODataSet_AddNodeScalar&
                       ,AddNodeVector=>IODataSet_AddNodeVector&
                       ,AddNodeTensor=>IODataSet_AddNodeTensor
end type IODataSet

contains

!! -----------------------------------------------------------------  
subroutine IODataSet_AllocateNode(dataSet)
  class(IODataSet), intent(inout) :: dataSet
  if (associated(dataSet%LastNode)) then
    allocate(dataSet%LastNode%NextNode)
    dataSet%LastNode => dataSet%LastNode%NextNode
  else
    allocate(dataSet%LastNode)
    dataSet%FirstNode => dataSet%LastNode
  end if
end subroutine IODataSet_AllocateNode

!! -----------------------------------------------------------------  
!! Append a new scalar field into the IO dataset.
subroutine IODataSet_AddNodeScalar(dataSet,nodeName,nodeData)
  class(IODataSet), intent(inout) :: dataSet
  character(len=*), intent(in) :: nodeName
  real(dp), target, intent(in) :: nodeData(:)
  call dataSet%AllocateNode()
  associate(node=>dataSet%LastNode)
    node%Name = nodeName
    node%Data => nodeData
    node%DataVector => null()
    node%DataTensor => null()
  end associate
end subroutine IODataSet_AddNodeScalar
!! Append a new vector field into the IO dataset.
subroutine IODataSet_AddNodeVector(dataSet,nodeName,nodeData)
  class(IODataSet), intent(inout) :: dataSet
  character(len=*), intent(in) :: nodeName
  real(dp), target, intent(in) :: nodeData(:,:)
  call dataSet%AllocateNode()
  associate(node=>dataSet%LastNode)
    node%Name = nodeName
    node%Data => null()
    node%DataVector => nodeData
    node%DataTensor => null()
  end associate
end subroutine IODataSet_AddNodeVector
!! Append a new tensor field into the IO dataset.
subroutine IODataSet_AddNodeTensor(dataSet,nodeName,nodeData)
  class(IODataSet), intent(inout) :: dataSet
  character(len=*), intent(in) :: nodeName
  real(dp), target, intent(in) :: nodeData(:,:,:)
  call dataSet%AllocateNode()
  associate(node=>dataSet%LastNode)
    node%Name = nodeName
    node%Data => null()
    node%DataVector => null()
    node%DataTensor => nodeData
  end associate
end subroutine IODataSet_AddNodeTensor

end module StormRuler_IO