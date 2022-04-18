add_executable(bottleneck_tests test_hera_bottleneck.cpp)

target_include_directories(bottleneck_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/../../extern/Catch2/src)

target_link_libraries(bottleneck_tests PRIVATE Catch2::Catch2WithMain)

# copy test data to the directory with bottleneck_tests binary after build
add_custom_command(TARGET bottleneck_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
                                               ${CMAKE_CURRENT_BINARY_DIR}/data)

add_test(NAME bottleneck COMMAND bottleneck_tests)
