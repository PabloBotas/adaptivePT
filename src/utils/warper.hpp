#ifndef __WARPER_HPP__
#define __WARPER_HPP__

#include "special_types.hpp"
#include "vector3.hpp"
#include "vector4.hpp"
#include <string>
#include <vector>

void 
apply_vf (std::vector< Vector4_t<float> >& p,
          const std::vector< Vector3_t<float> >& vf);

void
flip_positions_X (std::vector< Vector4_t<float> >& vec,
                  const CT_Dims_t dims);

void
flip_direction_X (std::vector< Vector4_t<float> >& vec);

void 
get_unitary_vector (std::vector< Vector3_t<float> >& r,
                    const std::vector< Vector4_t<float> >& p,
                    const std::vector< Vector4_t<float> >& p2);

std::vector< Vector3_t<float> >
get_vf_from_stdout (std::string str);

void
probe_vf (std::vector< Vector3_t<float> >& vf,
          const std::vector< Vector4_t<float> >& p,
          std::string vf_file);

std::vector< Vector3_t<float> >
probe_vf (const std::vector< Vector4_t<float> >& p,
          std::string vf_file);

void 
project_vector_on_plane (std::vector< Vector3_t<float> >& p,
                         const std::vector< Vector4_t<float> >& n);

std::string
to_location_str (const Vector3_t<float>& p, const bool last);

void
warp_data(std::vector< Vector4_t<float> >& endpoints,
          std::vector< Vector4_t<float> >& init_pos,
          const std::string vf_file,
          const CT_Dims_t& ct,
          std::vector< Vector4_t<float> > treatment_plane);

#endif