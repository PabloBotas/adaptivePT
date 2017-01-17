cmake_minimum_required(VERSION 3.2)

project(adaptSpotEnergies)
set (TARGET_NAME adaptSpotEnergies)
set (CMAKE_C_COMPILER "/usr/bin/gcc")
set (CMAKE_CXX_COMPILER "/usr/bin/clang++")
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")

set (PROJECT_ROOT_DIR ..)
set (CMAKE_INSTALL_PREFIX ".")
set (EXECUTABLE_OUTPUT_PATH "${PROJECT_ROOT_DIR}")
message ( STATUS "PROJECT_ROOT_DIR:       " ${PROJECT_ROOT_DIR} )
message ( STATUS "CMAKE_INSTALL_PREFIX:   " ${CMAKE_INSTALL_PREFIX} )
message ( STATUS "EXECUTABLE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

find_package (Boost COMPONENTS program_options REQUIRED)
include_directories ("${Boost_INCLUDE_DIRS}")

file(GLOB SRC_FILES src/*.cpp)
# message ( STATUS "SRC_FILES:              " ${SRC_FILES} )

add_executable (${TARGET_NAME} ${SRC_FILES})
target_link_libraries (${TARGET_NAME} ${Boost_LIBRARIES})
target_link_libraries(${TARGET_NAME} m)
