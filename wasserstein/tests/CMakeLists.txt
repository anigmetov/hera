add_executable(wasserstein_test  tests_main.cpp test_hera_wasserstein.cpp test_hera_wasserstein_pure_geom.cpp)
target_link_libraries(wasserstein_test PRIVATE ${libraries})

target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                         ${CMAKE_CURRENT_SOURCE_DIR}../extern/Catch2/src)

target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

add_test(NAME all COMMAND wasserstein_test)
