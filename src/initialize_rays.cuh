#ifndef __INITIALIZE_RAYS_CUH__
#define __INITIALIZE_RAYS_CUH__

#include "patient_parameters.hpp"
#include "helper_math.h"

#include <vector>

void create_virtual_source_buffers(const Patient_Parameters_t& pat,
                              std::vector<double4>& xbuffer,
                              std::vector<double4>& vxbuffer,
                              std::vector<short2>& ixbuffer);

void create_treatment_plane_buffers(const Patient_Parameters_t& pat,
                                    const Array4<double>& endpoints,
                                    const Array4<double>& init_pos,
                                    std::vector<double4>& xbuffer,
                                    std::vector<double4>& vxbuffer,
                                    std::vector<short2>& ixbuffer);

double3 iso_to_virtual_src_pos(double z, double2 SAD, double2 spot);

double2 virtual_src_to_iso_pos(double3 p, double2 SAD);
double2 virtual_src_to_iso_pos(double3 pos, double3 cos);
void virtual_src_to_iso_pos(Array4<double>& pos, SAD_t SAD);

double3 getDirection(double3 pos, double2 spot);

short2 get_beam_spot_id (size_t num, const std::vector<short>& spots_per_field);
#endif
