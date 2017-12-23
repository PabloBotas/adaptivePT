set (TARGET_NAME adaptive)

include_directories(. cuda cuda/cuda_helper
                    patient_data
                    phys_data
                    utils
                    volumes_representations)

# Files------------------------------------------------------
file (GLOB cpp_sources  ./*.cpp cuda/*.cpp patient_data/*.cpp phys_data/*.cpp utils/*.cpp volumes_representations/*.cpp)
file (GLOB cpp_headers  ./*.hpp cuda/*.hpp patient_data/*.hpp phys_data/*.hpp utils/*.hpp volumes_representations/*.hpp)
file (GLOB cuda_sources ./*.cu  cuda/*.cu  patient_data/*.cu  phys_data/*.cu  utils/*.cu  volumes_representations/*.cu )
file (GLOB cuda_headers ./*.cuh cuda/*.cuh patient_data/*.cuh phys_data/*.cuh utils/*.cuh volumes_representations/*.cuh)
set (EXCLUDE_DIR "extra")
foreach (TMP_PATH ${cpp_sources})
    string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM ${cpp_sources} ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)
foreach (TMP_PATH ${cpp_headers})
    string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM ${cpp_headers} ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)
foreach (TMP_PATH ${cuda_sources})
    string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM ${cuda_sources} ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)
foreach (TMP_PATH ${cuda_headers})
    string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM ${cuda_headers} ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)
# message ( STATUS "cpp_sources:    " ${cpp_sources} )
# message ( STATUS "cpp_headers:    " ${cpp_headers} )
# message ( STATUS "cuda_sources:   " ${cuda_sources} )
# message ( STATUS "cuda_headers:   " ${cuda_headers} )

# Create objects--------------------------------------------------
# cuda_compile (cuda_obj ${cuda_sources})
# compile (cpp_obj OBJECT ${cpp_sources})

cuda_add_executable (${TARGET_NAME} ${cuda_sources} ${cpp_sources})
# add_executable (${TARGET_NAME} ${cpp_sources} ${cpp_headers} ${cuda_obj})

# Boost------------------------------------------------------
find_package (Boost COMPONENTS program_options system REQUIRED)
include_directories ("${Boost_INCLUDE_DIRS}")
target_link_libraries (${TARGET_NAME} ${Boost_LIBRARIES})

# Math library ------------------------------------------------------
target_link_libraries(${TARGET_NAME} m)
