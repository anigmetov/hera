cmake_minimum_required(VERSION 3.12)
project(hera VERSION 1.1.1)

# C++ 14 required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)

# Hera is header-only library
# adding this target just makes it possible to automatically add
# correct include directories by target_link_libraries(user_prog PRIVATE hera)
add_library(hera INTERFACE)
target_include_directories(hera INTERFACE include/ extern/)

# TODO: detect system Catch2?
add_subdirectory(extern/Catch2)
enable_testing()

add_subdirectory(bottleneck)
add_subdirectory(wasserstein)
add_subdirectory(matching)
add_subdirectory(bindings/python)
