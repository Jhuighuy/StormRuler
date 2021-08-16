@echo off

rmdir /S /Q bin
rmdir /S /Q src_gen

rmdir /S /Q CMakeFiles
del CMakeCache.txt
del cmake_install.cmake

del build.ninja
del .ninja*

del Makefile
