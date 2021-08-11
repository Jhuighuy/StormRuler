@echo off

if NOT EXIST "CMakeLists.txt" (     
  echo ERROR: This script should be executed from the project root directory.
  pause
  exit 1
)

call scripts/CMake_Rebuild_IntelClassic.bat Release 
