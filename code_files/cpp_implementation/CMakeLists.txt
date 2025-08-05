cmake_minimum_required(VERSION 3.10)
project(option_solver)

# Compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")

# Directories
include_directories(include)

# Find all source files in src/
file(GLOB SOURCES "src/*.cpp")

# Create the executable
add_executable(option_solver ${SOURCES})

# Optional: custom commands like "make run"
add_custom_target(run
    COMMAND option_solver
    DEPENDS option_solver
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Optional: clean CSV data files like "make clean-data"
add_custom_target(clean-data
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/data/*.csv
)
