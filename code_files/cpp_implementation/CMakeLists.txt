cmake_minimum_required(VERSION 3.14)
project(option_pricing)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_executable(option_pricing ${SOURCES})
