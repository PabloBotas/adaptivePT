#ifndef __GPU_RAY_POSITIONING_CUH__
#define __GPU_RAY_POSITIONING_CUH__

#include "patient_parameters.hpp"
#include "special_types.hpp"
#include <vector>

__host__ void virtual_src_to_treatment_plane(const unsigned int& num,
                                             const std::vector<BeamAngles_t>& angles,
                                             const double3& ct_offsets);

__global__ void virtual_src_to_treatment_plane_kernel(const int num,
                                                      const double2 *angles,
                                                      const double3 ct_offsets);
__host__ void treatment_plane_to_virtual_src(Array4<double>& pos,
                                             const Array4<double>& pos2,
                                             const Patient_Parameters_t& pat);
__global__ void treatment_plane_to_virtual_src_kernel(const int num,
                                                      const int nbeams,
                                                      double4* pos_,
                                                      const double4* dir_,
                                                      const double2* angles,
                                                      const double4* plane_dir,
                                                      const double4* plane_pos,
                                                      const double3 ct_offsets);
void correct_offsets(const unsigned int& num,
                     const double3& offsets,
                     const double3& original_offsets);
Array4<double> offset_endpoints(const Array4<double>& orig_endpoints,
                               const double3& offsets,
                               const double3& original_offsets);
__global__ void correct_offsets_kernel(const int num,
                                       const double3 offsets,
                                       const double3 original_offsets);
__global__ void correct_offsets_kernel(const int num,
                                       double4* dev_orig_endpoints,
                                       const double3 offsets,
                                       const double3 original_offsets);

__device__ double4 ray_trace_to_CT_volume(const double4& p,
                                         const double4& v);
__device__ __host__ double4 ray_trace_to_CT_volume(const double4& p,
                                                  const double4& v,
                                                  const int3 nvox,
                                                  const double3 dvox);
__device__ __host__ double3 ray_trace_to_CT_volume(const double3& p,
                                                  const double3& v,
                                                  const int3 nvox,
                                                  const double3 dvox);

// __host__ Planes_t
// get_treatment_planes (const Patient_Parameters_t& pat,
//                       const std::vector<BeamAngles_t>& angles);

__device__ __host__ double3 ext_to_int_coordinates(double3 a);
__device__ __host__ double4 ext_to_int_coordinates(double4 a);
__device__ __host__ double3 int_to_ext_coordinates(double3 a);
__device__ __host__ double4 int_to_ext_coordinates(double4 a);
__device__ double4 rotate(const double4& p, const double& gantry, const double& couch);

#endif
