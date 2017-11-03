#ifndef __WARPER_HPP__
#define __WARPER_HPP__

#include "program_options.hpp"
#include "special_types.hpp"
#include "vector3.hpp"
#include "vector4.hpp"
#include <string>
#include <vector>

class Warper_t
{
public:
    Warper_t();
    Warper_t(const std::string vf_file,
             const std::string output_vf);
    void set(const std::string vf_file, const std::string output_vf);
    void apply_to_plan (Array4<double>& endpoints,
                        Array4<double>& init_pos,
                        const CT_Dims_t& ct,
                        Planes_t treatment_plane,
                        const std::vector<BeamAngles_t>& angles,
                        const std::vector<short>& spots_per_field,
                        const Warp_opts_t opts);
    Array4<double> apply_to_points (const Array4<double>& pos,
                                    const CT_Dims_t& ct);
    void print_vf (unsigned int n = 10);
    Vector3_t<double> vf_ave;
    Array3<double> vf_ave_planes;
    
private:
    // Constructor helpers
    void set_vf_origins ();
    // Main functions
    void apply_position_options (Warp_opts_t options, const std::vector<short>& spots_per_field);
    void apply_rigid_positions ();
    void apply_rigid_positions_per_beam (const std::vector<short>& spots_per_field);
    void probe (const Array4<double>& p, const CT_Dims_t& ct);
    void write_to_file (const Array4<double>& p,
                        const std::vector<short>& spots_per_field);
    void warp_points (Array4<double>& p);
    void warp_init_points (Array4<double>& init_pos, const Planes_t& pln,
                           const std::vector<short>& spots_per_field,
                           const std::vector<BeamAngles_t>& angles);
    void project_vf_on_plane (const Planes_t& pln,
                              const std::vector<short>& spots_per_field);

    // Utils
    void flip_positions_X (Array4<double>& vec, const CT_Dims_t dims);
    void flip_direction_X (Array4<double>& vec);
    std::string to_location_str (const Vector3_t<double>& p, const bool last);
    void set_average ();
    void set_probes (const Array4<double>& p);

    bool exp_file;
    std::string file;
    std::string output;
    Vector3_t<double> origins;
    Array3<double> vf;
    Array3<double> vf_planes;
    Array3<double> probes;
};



#endif
