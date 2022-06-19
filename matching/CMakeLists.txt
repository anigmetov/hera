include(GenerateExportHeader)

# Boost
find_package(Boost REQUIRED)

find_package(Threads)

set(libraries ${libraries} "stdc++fs" ${CMAKE_THREAD_LIBS_INIT} hera)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(libraries ${libraries} ${OpenMP_CXX_LIBRARIES})
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

foreach(exec matching_dist module_example)
    add_executable(${exec} "${exec}.cpp")
    target_include_directories(${exec} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../extern/opts/include)
    target_include_directories(${exec} SYSTEM PRIVATE ${Boost_INCLUDE_DIR})
    target_link_libraries(${exec} PUBLIC ${libraries})

    if(MSVC)
        target_compile_options(${exec} PRIVATE /W4 /WX)
    else()
        target_compile_options(${exec} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endforeach()

#if(HERA_BUILD_TESTS)
#    add_subdirectory(tests)
#endif()
