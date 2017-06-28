#ifndef __WARPER_HPP__
#define __WARPER_HPP__

#include "special_types.hpp"
#include "vector3.hpp"
#include "vector4.hpp"
#include <string>
#include <vector>

void 
apply_vf (Array4<float>& p,
          const Array3<float>& vf);

void
export_vf(const Array3<float>& vf,
          const Array4<float>& p,
          const std::string& file,
          const std::vector<short>& spots_per_field);
void
flip_positions_X (Array4<float>& vec,
                  const CT_Dims_t dims);

void
flip_direction_X (Array4<float>& vec);

void 
get_unitary_vector (Array3<float>& r,
                    const Array4<float>& p,
                    const Array4<float>& p2);

Array3<float>
get_vf_from_stdout (std::string str);

void
probe_vf (Array3<float>& vf,
          const Array4<float>& p,
          std::string vf_file);

Array3<float>
probe_vf (const Array4<float>& p,
          std::string vf_file);

void 
project_vector_on_plane (Array3<float>& p,
                         const Array4<float>& n,
                         const std::vector<short>& spots_per_field);

std::string
to_location_str (const Vector3_t<float>& p, const bool last);

void
warp_data (Array4<float>& endpoints,
           Array4<float>& init_pos,
           const std::string vf_file,
           const std::string output_vf,
           const CT_Dims_t& ct,
           Array4<float> treatment_plane,
           const std::vector<short>& spots_per_field);

#endif
