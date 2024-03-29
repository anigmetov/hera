add_subdirectory(pybind11)

set(CMAKE_CXX_STANDARD 14)

if (CMAKE_LIBRARY_OUTPUT_DIRECTORY)        # to work with pip
    set (MODULE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
else()
    set (MODULE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

file (MAKE_DIRECTORY "${MODULE_OUTPUT_DIRECTORY}/hera")
file (GLOB HERA_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/hera/*.py")

add_custom_target  (pyhera ALL ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_CURRENT_SOURCE_DIR}/hera ${MODULE_OUTPUT_DIRECTORY}/hera
                    DEPENDS ${HERA_PYTHON})

pybind11_add_module         (_hera hera.cpp hera_ws_geom.cpp)
target_link_libraries       (_hera PRIVATE hera)

set_target_properties       (_hera PROPERTIES OUTPUT_NAME hera/_hera)
