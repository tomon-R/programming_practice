cmake_minimum_required(VERSION 3.10)
project(FEMProject)

set(CMAKE_CXX_STANDARD 17)

include_directories(include model solver)

file(GLOB MODEL_SOURCES "model/*.cpp")
file(GLOB SOLVER_SOURCES "solver/*.cpp")
file(GLOB SRC_SOURCES "src/*.cpp")

# Define the output directory for generated files
set(OUTPUT_DIR ${CMAKE_BINARY_DIR}/output)

# Ensure the output directory exists at configure time
file(MAKE_DIRECTORY ${OUTPUT_DIR})

# Pass the output directory to the compiler as a definition
add_definitions(-DOUTPUT_DIR="${OUTPUT_DIR}")

# Add main executable
add_executable(FEMProject main/main.cpp ${MODEL_SOURCES} ${SOLVER_SOURCES} ${SRC_SOURCES})
add_executable(Iso_parametric main/iso_parametric.cpp ${MODEL_SOURCES} ${SOLVER_SOURCES} ${SRC_SOURCES})

# Add test executable(s)
add_executable(Test1 test/test1.cpp ${MODEL_SOURCES} ${SOLVER_SOURCES} ${SRC_SOURCES})
add_executable(Test2 test/test2.cpp ${MODEL_SOURCES} ${SOLVER_SOURCES} ${SRC_SOURCES})

# Pass OUTPUT_DIR definition to all executables
target_compile_definitions(FEMProject PRIVATE OUTPUT_DIR="${OUTPUT_DIR}")
target_compile_definitions(Iso_parametric PRIVATE OUTPUT_DIR="${OUTPUT_DIR}")
target_compile_definitions(Test1 PRIVATE OUTPUT_DIR="${OUTPUT_DIR}")
target_compile_definitions(Test2 PRIVATE OUTPUT_DIR="${OUTPUT_DIR}")