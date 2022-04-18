find_package (Threads)
set (libraries ${libraries}  ${CMAKE_THREAD_LIBS_INIT})

add_executable(bottleneck_dist bottleneck_dist.cpp)
target_include_directories(bottleneck_dist PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(bottleneck_dist PUBLIC ${libraries})

add_subdirectory(tests)
