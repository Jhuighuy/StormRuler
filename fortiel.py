#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# Copyright (C) 2021 Oleg Butakov
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights  to use,
# copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

import sys
import re
import json
from json import JSONEncoder
from typing import \
  Any, List, Dict, Union, \
  Callable, Optional, Pattern, Match

sys.dont_write_bytecode = True

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


class TielError(Exception):
  '''Abstract syntax tree parse error.
  '''
  def __init__(self, message: str,
               filePath: str, line: str, lineNumber: int) -> None:
    super().__init__()
    self.message: str = message
    self.filePath: str = filePath
    self.line: str = line
    self.lineNumber: int = lineNumber

  def __str__(self) -> str:
    message = f'{self.filePath}:{self.lineNumber}:1:\n\n' \
            + f'{self.line}\n1\n' \
            + f'Fatal Error: {self.message}'
    return ''.join(message)


class TielEndError(TielError):
  '''Unexpected end of file.
  '''
  def __init__(self, filePath: str, lineNumber: int):
    super().__init__('unexpected end of file', filePath, '', lineNumber)


class TielDirError(TielError):
  '''Unexpected directive error.
  '''
  def __init__(self, directive: str,
               filePath: str, line: str, lineNumber: int):
    super().__init__(
      f'unexpected directive `{directive}`',
      filePath, line, lineNumber)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


def _regExp(pattern: str) -> Pattern[str]:
  return re.compile(pattern, re.IGNORECASE)


_DIR = _regExp(
  r'#\s*fpp\s+(?P<dir>(?P<head>\w+)(\s+.+)?)(\s+\!.*)?$')

_USE = _regExp(
  r'^(?P<dir>use|include)\s+(?P<path>(\".+\")|(\'.+\')|(\<.+\>))$')

_IF = _regExp(
  r'^if\s*\((?P<cond>.+)\)\s*then$')
_ELSE_IF = _regExp(
  r'^else\s*if\s*\((?P<cond>.+)\)\s*then$')
_ELSE = _regExp(r'^else$')
_END_IF = _regExp(r'^end\s*if$')

_DO = _regExp(
  r'^do\s+(?P<index>[a-zA-Z]\w*)\s*=\s*(?P<bounds>.*)$')
_END_DO = _regExp(r'^end\s*do$')


class TielTree:
  '''An abstract syntax tree.
  '''
  def __init__(self, filePath: str) -> None:
    self.filePath: str = filePath
    self.rootNodes: List[TielTreeNode] = []

  def __str__(self) -> str:
    class Encoder(JSONEncoder):
      def default(self, obj):
        return obj.__dict__

    string = json.dumps(self, indent=2, cls=Encoder)
    return string


class TielTreeNode:
  '''An abstract syntax tree node.
  '''
  def __init__(self, filePath: str, lineNumber: int) -> None:
    self.filePath: str = filePath
    self.lineNumber: int = lineNumber


class TielTreeNodeLineBlock(TielTreeNode):
  '''The block of regular lines syntax tree node.
  '''
  def __init__(self, filePath: str, lineNumber: int) -> None:
    super().__init__(filePath, lineNumber)
    self.lines: List[str] = []


class TielTreeNodeUse(TielTreeNode):
  '''The USE/INCLUDE directive syntax tree node.
  '''
  def __init__(self, filePath: str, lineNumber: int) -> None:
    super().__init__(filePath, lineNumber)
    self.pathToInclude: str = ''
    self.emitLineBlocks: bool = False


class TielTreeNodeIfElseEnd(TielTreeNode):
  '''The IF/ELSE IF/ELSE/END IF directive syntax tree node.
  '''
  def __init__(self, filePath: str, lineNumber: int) -> None:
    super().__init__(filePath, lineNumber)
    self.condition: str = ''
    self.thenBranch: List[TielTreeNode] = []
    self.elseIfBranches: List[TielTreeNodeElseIf] = []
    self.elseBranch: List[TielTreeNode] = []


class TielTreeNodeElseIf(TielTreeNode):
  '''The ELSE IF directive syntax tree node.'''
  def __init__(self, filePath: str, lineNumber: int) -> None:
    super().__init__(filePath, lineNumber)
    self.condition: str = ''
    self.branchBody: List[TielTreeNode] = []


class TielTreeNodeDoEnd(TielTreeNode):
  '''The DO/END DO directive syntax tree node.
  '''
  def __init__(self, filePath: str, lineNumber: int) -> None:
    super().__init__(filePath, lineNumber)
    self.indexName: str = ''
    self.bounds: str = ''
    self.loopBody: List[TielTreeNode] = []


class TielParser:
  '''Abstract syntax tree parser.
  '''
  def __init__(self, filePath: str, lines: List[str]) -> None:
    self._filePath: str = filePath
    self._lines: List[str] = lines
    self._curLineIndex: int = 0

  def _curLine(self) -> str:
    return self._lines[self._curLineIndex]

  def _curLineNumber(self) -> int:
    return self._curLineIndex + 1

  def _advanceLine(self) -> None:
    self._curLineIndex += 1

  def _matchesEnd(self) -> bool:
    return self._curLineIndex >= len(self._lines)

  def _matchLine(self, regExp: Pattern[str]) -> Match[str]:
    match = self._matchesLine(regExp)
    if match is None:
      raise RuntimeError('expected match')
    self._advanceLine()
    return match

  def _matchesLine(self, *regExpList: Pattern[str]) -> Optional[Match[str]]:
    if self._matchesEnd():
      raise TielEndError(self._filePath, self._curLineNumber())
    for regExp in regExpList:
      match = regExp.match(self._curLine())
      if match is not None:
        return match
    return None

  def _matchDirective(self, regExp: Pattern[str]) -> Match[str]:
    directive = self._matchLine(_DIR).group('dir')
    match = regExp.match(directive)
    if match is None:
      raise TielDirError(directive, self._filePath,
                         self._curLine(), self._curLineNumber())
    return match

  def _matchesDirective(self, *regExpList: Pattern[str]) -> Optional[Match[str]]:
    dirMatch = self._matchesLine(_DIR)
    if dirMatch is not None:
      directive = dirMatch.group('dir')
      for regExp in regExpList:
        match = regExp.match(directive)
        if match is not None:
          return match
    return None

  def parse(self) -> TielTree:
    '''Parse the source lines.'''
    tree = TielTree(self._filePath)
    while not self._matchesEnd():
      tree.rootNodes.append(self._parseSingle())
    return tree

  def _parseSingle(self) -> TielTreeNode:
    '''Parse a directive or a line block.'''
    if self._matchesLine(_DIR):
      return self._parseDirective()
    return self._parseLineBlock()

  def _parseLineBlock(self) -> TielTreeNodeLineBlock:
    '''Parse a line block.'''
    node = TielTreeNodeLineBlock(self._filePath,
                                 self._curLineNumber())
    while True:
      node.lines.append(self._curLine())
      self._advanceLine()
      if self._matchesEnd() or self._matchesLine(_DIR):
        break
    return node

  def _parseDirective(self) -> TielTreeNode:
    '''Parse a directive.'''
    if self._matchesDirective(_USE):
      return self._parseDirectiveUse()
    elif self._matchesDirective(_IF):
      return self._parseDirectiveIfElseEnd()
    elif self._matchesDirective(_DO):
      return self._parseDirectiveDoEnd()
    directive = self._matchesLine(_DIR)['dir']
    raise TielDirError(directive, self._filePath,
                       self._curLine(), self._curLineNumber())

  def _parseDirectiveUse(self) -> TielTreeNodeUse:
    '''Parse USE/INCLUDE directives.'''
    node = TielTreeNodeUse(self._filePath,
                           self._curLineNumber())
    directive, node.pathToInclude \
      = self._matchDirective(_USE).group('dir', 'path')
    node.pathToInclude = node.pathToInclude[1:-1]
    if directive.lower() == 'include':
      node.emitLineBlocks = True
    return node

  def _parseDirectiveIfElseEnd(self) -> TielTreeNodeIfElseEnd:
    '''Parse IF/ELSE IF/ELSE/END IF directives.'''
    node = TielTreeNodeIfElseEnd(self._filePath,
                                 self._curLineNumber())
    node.condition = self._matchDirective(_IF)['cond']
    while not self._matchesDirective(_ELSE_IF, _ELSE, _END_IF):
      node.thenBranch.append(self._parseSingle())
    while not self._matchesDirective(_ELSE, _END_IF):
      elseIfNode = TielTreeNodeElseIf(self._filePath,
                                      self._curLineNumber())
      elseIfNode.condition = self._matchDirective(_ELSE_IF)['cond']
      while not self._matchesDirective(_ELSE_IF, _ELSE, _END_IF):
        elseIfNode.branchBody.append(self._parseSingle())
      node.elseIfBranches.append(elseIfNode)
    if self._matchesDirective(_ELSE):
      self._advanceLine()
      while not self._matchesDirective(_END_IF):
        node.elseBranch.append(self._parseSingle())
    self._matchDirective(_END_IF)
    return node

  def _parseDirectiveDoEnd(self) -> TielTreeNodeDoEnd:
    '''Parse DO/END DO directives.'''
    node = TielTreeNodeDoEnd(self._filePath,
                             self._curLineNumber())
    node.indexName, node.bounds \
      = self._matchDirective(_DO).group('index', 'bounds')
    while not self._matchesDirective(_END_DO):
      node.loopBody.append(self._parseSingle())
    self._matchDirective(_END_DO)
    return node


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


class TielEvaluator:
  '''Abstract syntax tree evaluator.
  '''
  def __init__(self, scope=None):
    if scope is None:
      scope = {}
    self._scope: Dict[str, Any] = scope

  def _evalExpr(self, expression: str):
    return eval(expression, dict(self._scope))

  def eval(self,
           nodeOrTree: Union[TielTree, TielTreeNode],
           callback: Callable[[str], None]) -> None:
    '''Evaluate the syntax tree or the syntax tree node.'''
    if isinstance(nodeOrTree, TielTree):
      tree: TielTree = nodeOrTree
      self._evalNodeList(tree.rootNodes, callback)
    else:
      node: TielTreeNode = nodeOrTree
      self._evalNodeList(node, callback)

  def _evalNodeList(self,
                    nodes: Union[TielTreeNode, List[TielTreeNode]],
                    callback: Callable[[str], None]) -> None:
    '''Evaluate the syntax tree node or a list of nodes.'''
    if not isinstance(nodes, list):
      nodes = [nodes]
    for node in nodes:
      if isinstance(node, TielTreeNodeLineBlock):
        self._evalLineBlock(node, callback)
      elif isinstance(node, TielTreeNodeUse):
        self._evalUse(node, callback)
      elif isinstance(node, TielTreeNodeIfElseEnd):
        self._evalIfElseEnd(node, callback)
      elif isinstance(node, TielTreeNodeDoEnd):
        self._evalDoEnd(node, callback)
      else:
        raise RuntimeError(node.__class__.__name__)

  def _evalLineBlock(self,
                     node: TielTreeNodeLineBlock,
                     callback: Callable[[str], None]) -> None:
    '''Evaluate line block.'''
    callback(f'# {node.lineNumber} "{node.filePath}"')
    for lineNumber, line in enumerate(node.lines, start=node.lineNumber):
      self._evalLine(line, lineNumber, node.filePath, callback)

  def _evalLine(self,
                line: str, lineNumber: int, filePath: str,
                callback: Callable[[str], None]) -> None:
    '''Evaluate in-line substitutions.'''
    line = re.sub(r'{(?P<expr>.+)}',
                  lambda x: str(self._evalExpr(x['expr'])), line)
    line = re.sub(r'@(?P<expr>:(\s*,)?)',
                  lambda x: str(self._scope['__INDEX__']*x['expr']), line)
    callback(line)

  def _evalUse(self,
               node: TielTreeNodeUse,
               callback: Callable[[str], None]) -> None:
    '''Evaluate USE/INCLUDE directive.'''
    print(node.__class__.__name__)

  def _evalIfElseEnd(self,
                     node: TielTreeNodeIfElseEnd,
                     callback: Callable[[str], None]) -> None:
    '''Evaluate IF/ELSE IF/ELSE/END IF node.'''
    condition: bool = self._evalExpr(node.condition)
    if condition:
      self._evalNodeList(node.thenBranch, callback)
    else:
      for elseIfNode in node.elseIfBranches:
        condition: bool = self._evalExpr(elseIfNode.condition)
        if condition:
          self._evalNodeList(node.thenBranch, callback)
          break
      else:
        self._evalNodeList(node.elseBranch, callback)

  def _evalDoEnd(self,
                 node: TielTreeNodeDoEnd,
                 callback: Callable[[str], None]) -> None:
    '''Evaluate DO/END DO node.'''
    bounds = self._evalExpr(node.bounds)
    start, stop = bounds
    for index in range(start, stop+1):
      self._scope[node.indexName] = index
      self._scope['__INDEX__'] = index
      self._evalNodeList(node.loopBody, callback)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


def tielPreprocess(filePath: str,
                   outputFilePath: str) -> None:
  '''Preprocess the source file.'''
  with open(filePath, 'r') as fp:
    lines = fp.read().splitlines()
  tree = TielParser(filePath, lines).parse()
  with open(outputFilePath, 'w') as fp:
    scope = {'NumRanks': 1}
    TielEvaluator(dict(scope)).eval(tree,
                                    lambda x: print(x, file=fp))


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

def _tielMain() -> None:
  filePath = sys.argv[1]
  outputFilePath = sys.argv[2]
  tielPreprocess(filePath, outputFilePath)


if __name__ == '__main__':
  _tielMain()
