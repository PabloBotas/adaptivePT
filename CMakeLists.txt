cmake_minimum_required(VERSION 3.2)

project(adaptive)
set (TARGET_NAME adaptive)

set (PROJECT_ROOT_DIR "../adaptive")
set (CMAKE_INSTALL_PREFIX ".")
set (EXECUTABLE_OUTPUT_PATH "${PROJECT_ROOT_DIR}")
set (SEPARATE_COMPILATION ON)
message ( STATUS "PROJECT_ROOT_DIR:       " ${PROJECT_ROOT_DIR} )
message ( STATUS "CMAKE_INSTALL_PREFIX:   " ${CMAKE_INSTALL_PREFIX} )
message ( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )
include_directories(src
                    src/cuda src/cuda/cuda_helper
                    src/patient_data
                    src/phys_data
                    src/utils
                    src/volumes_representations)

#-C++ compiler-------------------------------------------------------
set (CMAKE_C_COMPILER "/usr/bin/gcc")
set (CMAKE_CXX_COMPILER "/usr/bin/g++")
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-variadic-macros -Wno-language-extension-token")

#-Cuda compiler------------------------------------------------------
find_package(CUDA REQUIRED)
include_directories(src/cuda_helper)
set (CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
set (CUDA_PROPAGATE_HOST_FLAGS ON)
set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++11 -gencode arch=compute_35,code=sm_35 -D_FORCE_INLINES)
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    if (DEVICE_DEBUGGING)
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -G)
    endif()
endif()
if (SEPARATE_COMPILATION)
    set(CUDA_SEPARABLE_COMPILATION ON)
endif()
option(REGISTERS_INFO "Outputs register and memory information per kernel" OFF)
if (REGISTERS_INFO)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --ptxas-options=-v)
endif()
option(COMPILE_SM20 "Compiles for compute capability of 2.0" OFF)
if (COMPILE_SM20)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_20,code=sm_20)
endif()
# set (CUDA_SDK_ROOT_DIR "/home/gmoc/NVIDIA_GPU_Computing_SDK")
# include_directories ("${CUDA_SDK_ROOT_DIR}")

 

# Definitions------------------------------------------------
option(DEBUG_READ_CT "Output a 32bit-float binary file with raw ct values" OFF)
if(DEBUG_READ_CT)
  add_definitions(-D__DEBUG_OUTPUT_READ_CT__)
endif()
get_filename_component(abs_bin_path ${CMAKE_INSTALL_PREFIX} ABSOLUTE)
add_definitions(-DBIN_PATH=${abs_bin_path} )
# option(DEBUG_RAYS "Use mode to export a CT volume with the ray paths" OFF)
# if(DEBUG_RAYS)
#   add_definitions(-D__OUTPUT_SCORER_VOLUME__)
# endif()


# Files------------------------------------------------------
file (GLOB cpp_sources  src/*.cpp src/**/*.cpp)
file (GLOB cpp_headers  src/*.hpp src/**/*.hpp)
file (GLOB cuda_sources src/*.cu  src/**/*.cu)
file (GLOB cuda_headers src/*.cuh src/**/*.cuh)
# message ( STATUS "cpp_sources:            " ${cpp_sources} )
# message ( STATUS "cpp_headers:            " ${cpp_headers} )
# message ( STATUS "cuda_sources:           " ${cuda_sources} )
# message ( STATUS "cuda_headers:           " ${cuda_headers} )

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
