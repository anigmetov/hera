project(matching_distance)
cmake_minimum_required(VERSION 3.5.1)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GenerateExportHeader)
find_package(Boost REQUIRED)

# Default to Release

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif (NOT CMAKE_BUILD_TYPE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
        SYSTEM ${BOOST_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../bottleneck/include)

set(CMAKE_CXX_STANDARD 14)


if (NOT WIN32)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${OpenMP_CXX_FLAGS} -Wall -Wextra ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -D_GLIBCXX_DEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g -ggdb")
endif (NOT WIN32)

file(GLOB BT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../bottleneck/include/*.h ${CMAKE_CURRENT_SOURCE_DIR}/../bottleneck/include/*.hpp)
file(GLOB MD_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)

file(GLOB SRC_TEST_FILES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
 
find_package(Threads)
set(libraries ${libraries} "stdc++fs" ${CMAKE_THREAD_LIBS_INIT})


find_package(OpenMP)
if (OPENMP_FOUND)
    set(libraries ${libraries} ${OpenMP_CXX_LIBRARIES})
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_executable(matching_dist "example/matching_dist.cpp" ${MD_HEADERS} ${BT_HEADERS} )
target_link_libraries(matching_dist PUBLIC ${libraries})

add_executable(module_example "example/module_example.cpp" ${MD_HEADERS} ${BT_HEADERS} )
target_link_libraries(module_example PUBLIC ${libraries})


add_executable(matching_distance_test ${SRC_TEST_FILES} ${BT_HEADERS} ${MD_HEADERS})
target_link_libraries(matching_distance_test PUBLIC ${libraries})
