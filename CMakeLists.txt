cmake_minimum_required(VERSION 3.2)

project(adaptSpotEnergies)
set (TARGET_NAME adaptSpotEnergies)



set (PROJECT_ROOT_DIR ..)
set (CMAKE_INSTALL_PREFIX ".")
set (EXECUTABLE_OUTPUT_PATH "${PROJECT_ROOT_DIR}")
set (SEPARATE_COMPILATION ON)
message ( STATUS "PROJECT_ROOT_DIR:       " ${PROJECT_ROOT_DIR} )
message ( STATUS "CMAKE_INSTALL_PREFIX:   " ${CMAKE_INSTALL_PREFIX} )
message ( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

#-C++ compiler-------------------------------------------------------
set (CMAKE_C_COMPILER "/usr/bin/clang")
set (CMAKE_CXX_COMPILER "/usr/bin/clang++")
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Wno-variadic-macros -Wno-language-extension-token")

#-Cuda compiler------------------------------------------------------
find_package(CUDA REQUIRED)
include_directories(src/cuda_helper)
set (CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
set (CUDA_PROPAGATE_HOST_FLAGS ON)
set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_20,code=sm_20 -D_FORCE_INLINES)
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    if (DEVICE_DEBUGGING)
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -G)
    endif()
endif()
if (SEPARATE_COMPILATION)
    set(CUDA_SEPARABLE_COMPILATION ON)
endif()
# set (CUDA_SDK_ROOT_DIR "/home/gmoc/NVIDIA_GPU_Computing_SDK")
# include_directories ("${CUDA_SDK_ROOT_DIR}")

# Definitions------------------------------------------------
option(DEBUG_SCORER "Use debug mode to export a CT volume with the ray paths" OFF)
if(DEBUG_SCORER)
  add_definitions(-D__DEBUG_MODE__)
endif()


# Files------------------------------------------------------
file (GLOB cpp_sources  src/*.cpp)
file (GLOB cpp_headers  src/*.hpp)
file (GLOB cuda_sources src/*.cu)
file (GLOB cuda_headers src/*.cuh)
# message ( STATUS "SRC_FILES:              " ${SRC_FILES} )

# Create objects--------------------------------------------------
# cuda_compile (cuda_obj ${cuda_sources})
# compile (cpp_obj OBJECT ${cpp_sources})

cuda_add_executable (${TARGET_NAME} ${cuda_sources} ${cpp_sources})
# add_executable (${TARGET_NAME} ${cpp_sources} ${cpp_headers} ${cuda_obj})

# Boost------------------------------------------------------
find_package (Boost COMPONENTS program_options REQUIRED)
include_directories ("${Boost_INCLUDE_DIRS}")
target_link_libraries (${TARGET_NAME} ${Boost_LIBRARIES})

# Math library ------------------------------------------------------
target_link_libraries(${TARGET_NAME} m)
