cmake_minimum_required(VERSION 3.23)
project(optimization)

set(CMAKE_CXX_STANDARD 14)

add_executable(optimization
        DRIVER_CODE.cpp
        classes/matrix.cpp
        libs/matrix.h
        classes/ode_solver.cpp
        libs/ode_solver.h
        classes/opt_alg.cpp
        libs/opt_alg.h
        classes/solution.cpp
        libs/solution.h
        classes/user_funs.cpp
        libs/user_funs.h
        )
