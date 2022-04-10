add_executable(wasserstein_dist wasserstein_dist.cpp)
target_link_libraries(wasserstein_dist PUBLIC ${libraries})

add_executable(wasserstein_dist_dipha wasserstein_dist_dipha.cpp)
target_link_libraries(wasserstein_dist_dipha PUBLIC ${libraries})

# pure geometric version, arbitrary dimension
add_executable(wasserstein_dist_point_cloud wasserstein_dist_point_cloud.cpp)
target_link_libraries(wasserstein_dist_point_cloud PUBLIC ${libraries})

if(MSVC)
    target_compile_options(wasserstein_dist             PRIVATE /W4 /WX)
    target_compile_options(wasserstein_dist_dipha       PRIVATE /W4 /WX)
    target_compile_options(wasserstein_dist_point_cloud PRIVATE /W4 /WX)
else()
    target_compile_options(wasserstein_dist             PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(wasserstein_dist_dipha       PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(wasserstein_dist_point_cloud PRIVATE -Wall -Wextra -Wpedantic)
endif()


