cmake_minimum_required(VERSION 3.10)
project(src C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(src
        adjacency.c
        brain.c
        brain.h
        constants.h
        diffusion.c
        diffusion.h
        matops.c
        matops.h
        nvu.c
        nvu.h
        run_parbrain.c
        solver.c
        solver.h)
