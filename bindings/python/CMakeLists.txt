add_subdirectory(pybind11)

if (CMAKE_LIBRARY_OUTPUT_DIRECTORY)        # to work with pip
    set (MODULE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
else()
    set (MODULE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

#file (MAKE_DIRECTORY "${MODULE_OUTPUT_DIRECTORY}/hera")
#file (GLOB HERA_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/hera/*.py")
#add_custom_target           (hera ALL ${CMAKE_COMMAND} -E copy_directory
#                                       ${CMAKE_CURRENT_SOURCE_DIR}/hera ${MODULE_OUTPUT_DIRECTORY}/hera
#                                       DEPENDS ${HERA_PYTHON})

pybind11_add_module         (hera hera.cpp)
#target_link_libraries       (_hera PRIVATE ${libraries})
target_include_directories  (hera PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../wasserstein/include"
                                          "${CMAKE_CURRENT_SOURCE_DIR}/../../bottleneck/include")
#set_target_properties       (hera PROPERTIES OUTPUT_NAME hera/_hera)
