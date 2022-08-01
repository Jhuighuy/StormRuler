
project(SpherePack 
  VERSION 1.0 
  DESCRIPTION "Sphere Packing generator"
  LANGUAGES CXX)

file(GLOB CXX_SOURCES "./*.cpp")

add_executable(SpherePack ${CXX_SOURCES})

find_package(OpenGL REQUIRED)
target_link_libraries(SpherePack PRIVATE OpenGL::GL)
target_link_libraries(SpherePack PRIVATE OpenGL::GLU)

find_package(GLUT REQUIRED)
target_link_libraries(SpherePack PRIVATE GLUT::GLUT)

# Require C++20.
set_target_properties(
  SpherePack PROPERTIES 
  CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES)

# Link with the main library.
target_link_libraries(SpherePack PRIVATE StormRuler)

# Specify binary output directory.
set_target_properties(
  SpherePack PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
