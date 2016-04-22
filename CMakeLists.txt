cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

project(ParBrain)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-Wall -std=c99 -g")
set(CMAKE_LD_FLAGS ${CMAKE_LD_FLAGS} "-lm")

# message(STATUS ${CMAKE_C_FLAGS})

# Add local CMake directory to CMake modules path to load project-specific modules. 
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

find_package(CXSparse REQUIRED)
include_directories(${CXSPARSE_INCLUDE_DIRS})

find_package (VTK REQUIRED)
include (${VTK_USE_FILE})

set(HDRs src/solver.h src/brain.h src/matops.h src/nvu.h)
set(SRCs src/run_parbrain.c src/solver.c src/brain.c src/matops.c src/nvu.c src/adjacency.c)

add_executable(parBrainSim ${SRCs} ${HDRs})
add_executable (ConvertBinToVtu util/binToVtu/src/bin_to_vtu.cxx)

target_link_libraries(parBrainSim 
	${MPI_LIBRARIES}
	${CXSPARSE_LIBRARIES}
	m)

target_link_libraries(ConvertBinToVtu 
	${VTK_LIBRARIES})
