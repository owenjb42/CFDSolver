cmake_minimum_required(VERSION 3.12)
project(CFDSolver VERSION 0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ======== Raylib ========
add_subdirectory(dependencies/raylib)

# ======== CFD Solver ========
add_executable(CFDSolver CFDSolver/main.cpp 
			   CFDSolver/Interface.cpp)
target_link_libraries(CFDSolver PUBLIC raylib)