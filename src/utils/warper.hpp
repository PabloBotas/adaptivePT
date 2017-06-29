#ifndef __WARPER_HPP__
#define __WARPER_HPP__

#include "special_types.hpp"
#include "vector3.hpp"
#include "vector4.hpp"
#include <string>
#include <vector>

class Warper_t
{
public:
    Warper_t(const std::string vf_file,
             const std::string output_vf);
    void apply_to (Array4<float>& endpoints,
                   Array4<float>& init_pos,
                   const CT_Dims_t& ct,
                   Array4<float> treatment_plane,
                   const std::vector<short>& spots_per_field);
    
private:
    // Constructor helpers
    void set_vf_origins ();
    // Main functions
    void probe (const Array4<float>& p);
    void write_to_file (const Array4<float>& p,
                        const std::vector<short>& spots_per_field);
    void warp_points (Array4<float>& p);
    void project_vf_on_plane (const Array4<float>& n,
                              const std::vector<short>& spots_per_field);

    // Utils
    void flip_positions_X (Array4<float>& vec, const CT_Dims_t dims);
    void flip_direction_X (Array4<float>& vec);
    std::string to_location_str (const Vector3_t<float>& p, const bool last);

    bool exp_file;
    std::string file;
    std::string output;
    Vector3_t<float> origins;
    Array3<float> vf;
};



#endif
