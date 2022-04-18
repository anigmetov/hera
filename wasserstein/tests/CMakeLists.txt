add_executable(wasserstein_tests test_hera_wasserstein.cpp test_hera_wasserstein_pure_geom.cpp)
target_link_libraries(wasserstein_tests PRIVATE ${libraries})

target_include_directories(wasserstein_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_include_directories(wasserstein_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../extern/Catch2/src)
target_include_directories(wasserstein_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../extern)
target_include_directories(wasserstein_tests PRIVATE SYSTEM ${Boost_INCLUDE_DIR})

target_link_libraries(wasserstein_tests PRIVATE Catch2::Catch2WithMain)

# copy test data to the directory with bottleneck_tests binary after build
add_custom_command(TARGET wasserstein_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
                                               ${CMAKE_CURRENT_BINARY_DIR}/data)


add_test(NAME wasserstein COMMAND wasserstein_tests)
