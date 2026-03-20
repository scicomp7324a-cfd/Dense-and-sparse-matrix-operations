cmake_minimum_required(VERSION 3.16)

project(PressureVelocityCoupling LANGUAGES CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# including the folders as libraries
add_library(mesh INTERFACE)
target_include_directories(mesh INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Mesh
    )

add_library(sparse_operations INTERFACE)
target_include_directories(sparse_operations INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/SparseOperations
    )
target_link_libraries(sparse_operations INTERFACE mesh)

add_library(dense_operations INTERFACE)
target_include_directories(dense_operations INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/DenseOperations
)

add_library(linear_vector INTERFACE)
target_include_directories(linear_vector INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearVector
)

add_library(test INTERFACE)
target_include_directories(test INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Test
)

add_library(viz_library INTERFACE)
target_include_directories(viz_library INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/VisualizationLibrary
)



add_executable(dsm_app main.cpp)
find_package(BLAS REQUIRED)
target_include_directories(dsm_app PRIVATE /opt/homebrew/opt/openblas/include)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(dsm_app PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(dsm_app PRIVATE
        /W4
    )
endif()

target_link_libraries(dsm_app PRIVATE
    mesh
    sparse_operations
    dense_operations
    linear_vector
    test
    BLAS::BLAS
)

add_custom_command(TARGET dsm_app POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/Files
            $<TARGET_FILE_DIR:dsm_app>/Files
)

add_custom_command(TARGET dsm_app POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:dsm_app>/Output
)
