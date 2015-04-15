cmake_minimum_required (VERSION 2.6)

project (binToVtu)

find_package (VTK REQUIRED)
include (${VTK_USE_FILE})

add_executable (ConvertBinToVtu src/binToVtu.cxx)

target_link_libraries (ConvertBinToVtu ${VTK_LIBRARIES})

