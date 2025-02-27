cmake_minimum_required(VERSION 3.12)
project(CFDSolver VERSION 0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ======== Raylib ========
add_subdirectory(dependencies/raylib)

# ======== Raygui ========
add_library(raygui CFDSolver/raygui.cpp)
target_include_directories(raygui PUBLIC dependencies/raygui/src)
target_link_libraries(raygui PUBLIC raylib)

# ======== CFD Solver ========
add_executable(CFDSolver CFDSolver/main.cpp 
			   CFDSolver/Interface.cpp 
			   CFDSolver/Helper.cpp 
			   CFDSolver/BoundaryControl.cpp)
target_link_libraries(CFDSolver PUBLIC raylib raygui)