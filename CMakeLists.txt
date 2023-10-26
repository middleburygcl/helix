# setup the project
cmake_minimum_required(VERSION 3.5)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
project(helix)

# options
option(WITH_EGADS "include EGADS" TRUE)
option(WITH_TETGEN "include TetGen" TRUE)

# setup libraries and build type
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "-std=c++17")

include(cmake/dependencies.cmake)
include(cmake/setup.cmake)
include(cmake/flags.cmake)
include(cmake/debug.cmake)
include(cmake/coverage.cmake)

# setup include directories
include_directories(${HELIX_INCLUDE_DIRS})
if(PROJECT_IS_TOP_LEVEL)
  set(HELIX_INCLUDE ${HELIX_INCLUDE_DIRS})
else()
  set(HELIX_INCLUDE ${HELIX_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# build the libraries, tests and executables
add_subdirectory(helix)
