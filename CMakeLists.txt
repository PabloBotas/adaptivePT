cmake_minimum_required(VERSION 3.2)

project(adaptive)

# set (PROJECT_ROOT_DIR "../adaptive")
set (CMAKE_INSTALL_PREFIX ".")
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set (SEPARATE_COMPILATION ON)
# message ( STATUS "PROJECT_ROOT_DIR:       " ${PROJECT_ROOT_DIR} )
message ( STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX} )
message ( STATUS "OUTPUT_DIRECTORY:     " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )

#-C++ compiler-------------------------------------------------------
set (CMAKE_C_COMPILER "/usr/bin/gcc")
set (CMAKE_CXX_COMPILER "/usr/bin/g++")
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-variadic-macros")

#-Cuda compiler------------------------------------------------------
find_package(CUDA REQUIRED)
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

### Definitions------------------------------------------------
# Debug outputs
option(DEBUG_READ_CT "Output a 32bit-float binary file with raw ct values" OFF)
if(DEBUG_READ_CT)
    add_definitions(-D__DEBUG_OUTPUT_READ_CT__)
endif()
get_filename_component(abs_bin_path ${CMAKE_INSTALL_PREFIX} ABSOLUTE)
add_definitions(-DBIN_PATH=${abs_bin_path} )

option(DEBUG_INFLUENCE_MATRICES "Output the 64bit-float binary influence matrices between spots" OFF)
if(DEBUG_INFLUENCE_MATRICES)
    add_definitions(-D__DEBUG_INFLUENCE_MATRICES__)
endif()

# Stepping definitions
option(STEP_DEFAULT "Step using average at 1 sigma" ON)
if(STEP_DEFAULT)
    add_definitions(-D__STEP_DEFAULT__)
    add_definitions(-D__STEP_AVERAGE__)
endif()

option(STEP_CENTRAL_AXIS "Step using only the central axis" OFF)
if(STEP_CENTRAL_AXIS)
    add_definitions(-D__STEP_CENTRAL_AXIS__)
endif()

option(STEP_Q50 "Step using Q50 statistics" OFF)
if(STEP_Q50)
    add_definitions(-D__STEP_Q50__)
endif()

option(STEP_AVERAGE "Step using average at 1 sigma" OFF)
if(STEP_AVERAGE)
    add_definitions(-D__STEP_AVERAGE__)
endif()

option(STEP_AVERAGE2 "Step using average at 0.5 sigma" OFF)
if(STEP_AVERAGE2)
    add_definitions(-D__STEP_AVERAGE2__)
endif()

option(STEP_AVERAGE3 "Step using average at 0.25 sigma" OFF)
if(STEP_AVERAGE3)
    add_definitions(-D__STEP_AVERAGE3__)
endif()
# option(DEBUG_RAYS "Use mode to export a CT volume with the ray paths" OFF)
# if(DEBUG_RAYS)
#   add_definitions(-D__OUTPUT_SCORER_VOLUME__)
# endif()

# Range shifter definitions
set(RANGE_SHIFTER_THICKNESS "0.01" CACHE STRING "Thickness in cm of created range shifters (default
                                                 is 0.01). The WEPL will be adjusted with the
                                                 density.")
string(TOLOWER ${RANGE_SHIFTER_THICKNESS} temp)
if ("default" STREQUAL ${temp})
    set(RANGE_SHIFTER_THICKNESS "0.01")
elseif((NOT RANGE_SHIFTER_THICKNESS MATCHES "^[+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$") OR
    (RANGE_SHIFTER_THICKNESS LESS 0) OR
    (RANGE_SHIFTER_THICKNESS EQUAL 0))
    message(FATAL_ERROR "RANGE_SHIFTER_THICKNESS must be a number bigger than zero.")
endif()
add_definitions(-D__RANGE_SHIFTER_THICKNESS__=${RANGE_SHIFTER_THICKNESS})

set(RANGE_SHIFTER_DENSITY "1.15" CACHE STRING "Nominal density of the range shifter (default is
                                               1.15).")
string(TOLOWER ${RANGE_SHIFTER_DENSITY} temp)
if ("default" STREQUAL ${temp})
    set(RANGE_SHIFTER_DENSITY "1.15")
elseif((NOT RANGE_SHIFTER_DENSITY MATCHES "^[+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$") OR
    (RANGE_SHIFTER_DENSITY LESS 0) OR
    (RANGE_SHIFTER_DENSITY EQUAL 0))
    message(FATAL_ERROR "RANGE_SHIFTER_DENSITY must be a number bigger than zero.")
endif()
add_definitions(-D__RANGE_SHIFTER_DENSITY__=${RANGE_SHIFTER_DENSITY})

option(RANGE_SHIFTER_STEPS "If the range shifter should be discretized and how." ON)
add_definitions(-D__RANGE_SHIFTER_DISCRETE_MGH__)

add_subdirectory(src)
# add_subdirectory(src/extra)

# Opt4D ------------------------------------------------------
find_program (OPT4D_PATH_BINARY opt4D "/opt/opt4D-build")
