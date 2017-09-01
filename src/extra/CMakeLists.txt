# CREATE DIJ ----------------------------------------------
set (TARGET_NAME dij_creator)

include_directories(../cuda/cuda_helper
                    ../patient_data
                    ../utils
                    ../volumes_representations)
set (CREATE_DIJ_SRCS calc_dij.cu ../patient_data/tramp.cpp ../patient_data/spot.cpp
                                 ../utils/energy_range.cpp ../utils/data_table.cpp ../utils/spline.cpp)
cuda_add_executable (${TARGET_NAME} ${CREATE_DIJ_SRCS})
CUDA_ADD_CUBLAS_TO_TARGET( ${TARGET_NAME} )
# Boost
find_package (Boost COMPONENTS program_options REQUIRED)
include_directories ("${Boost_INCLUDE_DIRS}")
target_link_libraries (${TARGET_NAME} ${Boost_LIBRARIES})
# Math library
target_link_libraries(${TARGET_NAME} m)
