# add all the source files
set(SOURCES
  array2d.cpp
  elements.cpp
  geometry.cpp
  io.cpp
  library.cpp
  log.cpp
  mesh.cpp
  mesher.cpp
  util.cpp)

# create an object library for all the sources
add_library(helix_src OBJECT ${SOURCES})
add_definitions(-DHELIX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_features(helix_src PRIVATE cxx_std_17)
#target_compile_definitions(helix_src PUBLIC -DHELIX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
set_target_properties(helix_src PROPERTIES
            CXX_STANDARD 17
            CXX_EXTENSIONS OFF)

# create both shared and static libraries with same name
set(HELIX_LIBS helix_static helix_shared )
add_library(helix_shared SHARED $<TARGET_OBJECTS:helix_src> $<TARGET_OBJECTS:helix_external>)
set_target_properties(helix_shared PROPERTIES OUTPUT_NAME helix CLEAN_DIRECT_OUTPUT 1)
target_link_libraries(helix_shared dl)
if (UNIX AND NOT APPLE)
  target_link_libraries(helix_shared rt dl)
endif()
add_library(helix_static STATIC $<TARGET_OBJECTS:helix_src> $<TARGET_OBJECTS:helix_external>)
set_target_properties(helix_static PROPERTIES OUTPUT_NAME helix CLEAN_DIRECT_OUTPUT 1)

add_custom_target(helix_lib DEPENDS ${HELIX_LIBS})
target_link_libraries(helix_shared ${HELIX_EXTERNAL_LIBRARIES})
target_link_libraries(helix_static ${HELIX_EXTERNAL_LIBRARIES})

# create executables
add_executable(helix helix.cpp)
add_dependencies(helix helix_lib)
target_link_libraries(helix helix_static dl)

# function to add unit test files with a suffix in directories
function(add_test_files var type)
  set(sources)
  foreach(dir ${ARGN})
      file(GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${dir}/*_${type}.cpp")
      list(APPEND sources ${_sources})
  endforeach()
  set(${var} ${${var}} ${sources} PARENT_SCOPE)
endfunction()

add_test_files(TEST_FILES "test"
  ${CMAKE_CURRENT_SOURCE_DIR} 
  ${CMAKE_CURRENT_SOURCE_DIR}/math)

# set the directory for the test executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)

# create a target for each test source
foreach(test_src ${TEST_FILES})

  # replace the '/' with '_' and "../" with ""
  string(REPLACE "../" "" test_src_underscored ${test_src})
  string(REPLACE "/" "_" test_src_underscored ${test_src_underscored})

  # remove the extension
  get_filename_component(test_bin ${test_src_underscored} NAME_WE)

  # assign the name of the executable
  set(test_exe ${test_bin}_exe)

  # the actual executable: target for simply building
  add_executable(${test_exe} ${test_src})
  target_link_libraries(${test_exe} helix_shared ${HELIX_EXTERNAL_LIBRARIES})
  target_compile_definitions(${test_exe} PUBLIC -DSTANDALONE -DHELIX_FULL_UNIT_TEST=false)

  # create the target
  add_custom_target(${test_bin} command $<TARGET_FILE:${test_exe}> 1 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data)

  # targets for memory & optionally stack checking & gdb
  ADD_DEBUG_TARGETS(${test_bin} ${CMAKE_SOURCE_DIR}/data)

endforeach()

# create the full unit executable
add_executable(unit_exe tester.cpp ${TEST_FILES})
target_link_libraries(unit_exe helix_shared)
target_compile_definitions(unit_exe PUBLIC -DSTDOUT_REDIRECT="${CMAKE_BINARY_DIR}/unit_tests_output.txt" -DHELIX_FULL_UNIT_TEST=true)

# target for running all unit tests from a single executable
add_custom_target(unit command $<TARGET_FILE:unit_exe> WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data)

ADD_DEBUG_TARGETS(unit ${CMAKE_SOURCE_DIR})
ADD_COVERAGE_UT(unit_coverage unit)

# reset the directory for the executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
