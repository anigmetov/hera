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
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${OpenMP_CXX_FLAGS} -Wall pedantic -Wextra ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${OpenMP_CXX_FLAGS} -Wall -Wextra ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -D_GLIBCXX_DEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g -ggdb")
endif (NOT WIN32)

file(GLOB BT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../bottleneck/include/*.h ${CMAKE_CURRENT_SOURCE_DIR}/../bottleneck/include/*.hpp)
file(GLOB MD_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

file(GLOB SRC_TEST_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/*.cpp)

find_package(Threads)
set(libraries ${libraries} "stdc++fs" ${CMAKE_THREAD_LIBS_INIT} ${OpenMP_CXX_LIBRARIES})


find_package(OpenMP)
if (OPENMP_FOUND)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


#add_executable(matching_distance ${SRC_FILES} ${BT_HEADERS} ${MD_HEADERS})
add_executable(matching_distance "src/main.cpp"
        "src/box.cpp"
        "src/common_util.cpp"
        "src/persistence_module.cpp"
        "src/simplex.cpp"
        "src/bifiltration.cpp"
        "src/matching_distance.cpp"
        "src/dual_box.cpp"
        "src/dual_point.cpp"
        "include/box.h"
        "include/common_util.h"
        "include/persistence_module.h"
        "include/simplex.h"
        "include/bifiltration.h"
        "include/matching_distance.h"
        "include/dual_box.h"
        "include/dual_point.h"
        ${BT_HEADERS} include/cell_with_value.h src/cell_with_value.cpp)
target_link_libraries(matching_distance PUBLIC ${libraries})

#add_executable(matching_distance_test ${SRC_TEST_FILES} ${BT_HEADERS} ${MD_HEADERS})
add_executable(matching_distance_test ${SRC_TEST_FILES}
        "src/box.cpp"
        "src/common_util.cpp"
        "src/persistence_module.cpp"
        "src/simplex.cpp"
        "src/bifiltration.cpp"
        "src/matching_distance.cpp"
        "src/dual_box.cpp"
        "src/dual_point.cpp"
        "include/box.h"
        "include/common_util.h"
        "include/persistence_module.h"
        "include/simplex.h"
        "include/bifiltration.h"
        "include/matching_distance.h"
        "include/dual_box.h"
        "include/dual_point.h"
        ${BT_HEADERS} src/tests/test_common.cpp  src/common_util.cpp src/tests/test_matching_distance.cpp src/cell_with_value.cpp)
target_link_libraries(matching_distance_test PUBLIC ${libraries})

add_executable(test_generator "src/test_generator.cpp"
        "src/box.cpp"
        "src/common_util.cpp"
        "src/persistence_module.cpp"
        "src/simplex.cpp"
        "src/bifiltration.cpp"
        "src/matching_distance.cpp"
        "src/dual_box.cpp"
        "src/dual_point.cpp"
        "include/box.h"
        "include/common_util.h"
        "include/persistence_module.h"
        "include/simplex.h"
        "include/bifiltration.h"
        "include/matching_distance.h"
        "include/dual_box.h"
        "include/dual_point.h"
        ${BT_HEADERS} src/cell_with_value.cpp)
target_link_libraries(test_generator PUBLIC ${libraries})

#add_executable(matching_distance "src/main.cpp" "src/box.cpp" "src/common_util.cpp" "src/line.cpp" "src/persistence_module.cpp" ${BT_HEADERS} ${MD_HEADERS})
