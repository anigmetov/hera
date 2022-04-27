# Boost
find_package(Boost REQUIRED)

find_package(Threads)
set (libraries ${libraries} ${CMAKE_THREAD_LIBS_INIT})

foreach(wd wasserstein_dist wasserstein_dist_dipha wasserstein_dist_point_cloud)
    add_executable(${wd} "${wd}.cpp")

    target_include_directories(${wd} SYSTEM PRIVATE ${Boost_INCLUDE_DIR})
    target_include_directories(${wd} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                             ${CMAKE_CURRENT_SOURCE_DIR}/../extern
                                             ${CMAKE_CURRENT_SOURCE_DIR}/../extern/opts/include)

target_link_libraries(${wd} PUBLIC ${libraries})

    if(MSVC)
        target_compile_options(${wd} PRIVATE /W4 /WX)
    else()
        target_compile_options(${wd} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endforeach()

add_subdirectory(tests)
