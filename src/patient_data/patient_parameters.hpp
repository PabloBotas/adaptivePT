#ifndef __PATIENT_PARAMETERS_HPP__
#define __PATIENT_PARAMETERS_HPP__

#include "patient_parameters_parser.hpp"
#include "volume.hpp"
#include "special_types.hpp"

#include <string>
#include <vector>

class Patient_Parameters_t
{
public:
    Patient_Parameters_t(std::string directory);
    void set_treatment_planes();
    void print();
    
    std::string patient_dir;
    std::string input_dir;
    std::string tramp_dir;
    std::string machine;
    SAD_t virtualSAD;

    size_t nbeams;
    std::vector<std::string> beam_dirs;
    std::vector<std::string> run_dir;
    std::vector<std::string> beam_names;
    std::vector<std::string> tramp_files;
    std::vector<std::string> adapted_tramp_names;
    std::vector<std::string> geometric_tramp_names;
    std::vector<std::string> topas_files;

    ushort total_spots;
    std::vector<short> spots_per_field;
    std::vector<short> accu_spots_per_field;
    std::vector<std::vector<short>> energy_layers;
    std::vector<float> source_energies;
    std::vector<float> source_weights;
    std::vector<uint> n_fractions;

    CT_Dims_t ct;
    CT_Dims_t original_ct;
    std::string planning_ct_file;
    std::vector<Aperture_Dims_t>     apertures;
    std::vector<RangeShifter_Dims_t> range_shifters;

    std::vector<BeamAngles_t> angles;
    std::vector<float>        isocenter_to_beam_distance;
    Planes_t                  treatment_planes;

    void ct_to_int_coordinates();
    void ct_to_ext_coordinates();
    void update_geometry(const Volume_t& vol);

private:
    void exploreBeamDirectories();
    void parseTopasFiles();
    void getTopasGlobalParameters();
    void getTopasBeamParameters();
    void set_spots_data();
    void set_planning_CT_file();
    void set_new_tramps();

    void consolidate_originals();
    std::vector<BeamAngles_t> original_angles;
    std::vector<float>       original_isocenter_to_beam_distance;

//     std::string rel_beam_location  = "/input";
//     std::string rel_topas_location = "/run";
};

#endif
