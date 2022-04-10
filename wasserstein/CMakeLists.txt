project(wasserstein)
cmake_minimum_required(VERSION 3.5.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)

# Boost
find_package(Boost)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include SYSTEM ${Boost_INCLUDE_DIR})

find_package (Threads)
set (libraries ${libraries}     ${CMAKE_THREAD_LIBS_INIT})

add_subdirectory(example)
#add_subdirectory(tests)

# Tests
add_executable(wasserstein_test  ${CMAKE_CURRENT_SOURCE_DIR}/tests/tests_main.cpp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hera_wasserstein.cpp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hera_wasserstein_pure_geom.cpp  include/hera_infinity.h tests/tests_reader.h)
target_link_libraries(wasserstein_test PUBLIC ${libraries})
