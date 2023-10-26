set(extern_repositories)
function(add_extern_repository name)
	set(options "")
	set(one_value_args GIT_REPOSITORY FULL_HISTORY SKIP_CONFIG)
	set(multi_value_args "")
	cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})
	if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/${name})	
		if(ARG_GIT_REPOSITORY)
			set(clone_opts --recursive --depth=1 --single-branch)
			if(ARG_FULL_HISTORY)
				set(clone_opts --recursive)
			endif()
			set(fetch_get git)
			set(fetch_url ${ARG_GIT_REPOSITORY})
			set(fetch_arg clone ${clone_opts} ${ARG_GIT_REPOSITORY} ${CMAKE_CURRENT_SOURCE_DIR}/extern/${name})
		else()
			message(FATAL_ERROR "unknown repository type")
		endif()
		message(STATUS "fetching ${name} from ${fetch_url}")
		execute_process(COMMAND ${fetch_get} ${fetch_arg} RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
	endif()
	if(NOT ARG_SKIP_CONFIG)
		add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/${name} EXCLUDE_FROM_ALL)
	endif()
	list(APPEND extern_repositories "${CMAKE_CURRENT_SOURCE_DIR}/extern/${name}")
	set(extern_repositories ${extern_repositories} PARENT_SCOPE)
endfunction()

set(WITH_GMF_AIO TRUE)
set(WITH_GMF_FORTRAN FALSE)

add_extern_repository(fmt GIT_REPOSITORY "https://github.com/fmtlib/fmt")
add_extern_repository(libmeshb GIT_REPOSITORY "https://github.com/LoicMarechal/libMeshb" SKIP_CONFIG TRUE)
add_extern_repository(tinyobjloader GIT_REPOSITORY "https://github.com/tinyobjloader/tinyobjloader")
add_extern_repository(argparse GIT_REPOSITORY "https://github.com/p-ranav/argparse")
add_extern_repository(morton GIT_REPOSITORY "https://github.com/morton-nd/morton-nd")
add_extern_repository(stlext GIT_REPOSITORY "https://github.com/middpolymer/stlext" SKIP_CONFIG TRUE)
add_extern_repository(OpenNL GIT_REPOSITORY "https://github.com/middpolymer/geogram.psm.OpenNL" SKIP_CONFIG TRUE)
add_extern_repository(PCK GIT_REPOSITORY "https://github.com/middpolymer/geogram.psm.Predicates" SKIP_CONFIG TRUE)
add_extern_repository(tetgen GIT_REPOSITORY "https://github.com/libigl/tetgen.git")

if(WITH_EGADS)
	add_extern_repository(egads GIT_REPOSITORY "https://github.com/middpolymer/egads")
endif()

# utilities to clean up and update repositories
add_custom_target(clean_extern COMMAND rm -rf ${extern_repositories})
foreach(repo ${extern_repositories})
	string(REPLACE ${CMAKE_CURRENT_SOURCE_DIR}/extern/ "" name ${repo})
	add_custom_target(update_${name} COMMAND cd ${repo} && git pull)
endforeach()

# external repositories
set(external_libraries fmt argparse)
if(WITH_EGADS)
	set(external_libraries ${external_libraries} egads)
endif()

# TBB
find_package(TBB)
if(TBB_FOUND)
	message(STATUS "found TBB")
	set(external_libraries ${external_libraries} TBB::tbb)
	add_definitions(-DHELIX_WITH_TBB=1)
else()
	message(WARNING "did not find TBB")
	add_definitions(-DHELIX_WITH_TBB=0)
endif()

# OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
	message(STATUS "found OpenMP")
	set(external_libraries ${external_libraries} OpenMP::OpenMP_CXX)
	add_definitions(-DHELIX_WITH_OMP=1)
else()
	message(STATUS "did not find OpenMP")
	add_definitions(-DHELIX_WITH_OMP=0)
endif()

set(HELIX_EXTERNAL_LIBRARIES ${external_libraries} ${GL_LIBRARIES})

# set all include directories
set(HELIX_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/extern
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/libmeshb/sources
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/PCK
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/fmt/include
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/morton/include
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/tetgen
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/stlext
)

if(WITH_EGADS)
	set(HELIX_INCLUDE_DIRS ${HELIX_INCLUDE_DIRS}
		${CMAKE_CURRENT_SOURCE_DIR}/extern/egads/include
		${CMAKE_CURRENT_SOURCE_DIR}/extern/egads/src
		${EGADS_OCC_INCLUDE}
	)
	set(HELIX_EXTERNAL_LIBRARIES ${HELIX_EXTERNAL_LIBRARIES} ${EGADS_OCC_LIBS})
	add_definitions(-DHELIX_WITH_EGADS=1)
else()
	message(STATUS "not including EGADS")
	add_definitions(-DHELIX_WITH_EGADS=0)
endif()

set(EXTERN_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/extern/OpenNL/OpenNL_psm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/PCK/Predicates_psm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/extern/tinyobjloader/tiny_obj_loader.cc
	${CMAKE_CURRENT_SOURCE_DIR}/extern/libmeshb/sources/libmeshb7.c
)
add_library(helix_external OBJECT ${EXTERN_SOURCES})

if(WITH_TETGEN)
	set(HELIX_INCLUDE_DIRS ${HELIX_INCLUDE_DIRS}
		${CMAKE_CURRENT_SOURCE_DIR}/extern/tetgen
	)
	add_definitions(-DHELIX_WITH_TETGEN=1)
	target_compile_definitions(tetgen PRIVATE -DTETLIBRARY)
	target_compile_options(tetgen PRIVATE "-w")
	set(HELIX_EXTERNAL_LIBRARIES ${HELIX_EXTERNAL_LIBRARIES} tetgen)
else()
	message(STATUS "not including TetGen")
	add_definitions(-DHELIX_WITH_TETGEN=0)
endif()

set_target_properties(helix_external PROPERTIES COMPILE_FLAGS "-w")