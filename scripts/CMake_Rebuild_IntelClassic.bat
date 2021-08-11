@echo off

if NOT EXIST "CMakeLists.txt" (
  echo ERROR: This script should be execured from the project root directory.
  pause
  exit 1
)

set "BUILD_TYPE=Debug"
if "%1" == "Release" (
  set "BUILD_TYPE=%1"
)
echo build type is "%BUILD_TYPE%"

if EXIST "CMakeCache.txt" (
  del "CMakeCache.txt"
)
set "FC=ifort"
set "CXX=cl"
cmake -GNinja -DCMAKE_BUILD_TYPE=%BUILD_TYPE% .
if %ErrorLevel% == 0 (
  ninja clean
  ninja
)

pause
