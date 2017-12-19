#ifndef __INFLUENCE_MANAGER_HPP__
#define __INFLUENCE_MANAGER_HPP__

#include "command_line_parser.hpp"
#include "patient_parameters.hpp"
#include "program_options.hpp"
#include "vector3.hpp"
#include "vector4.hpp"
#include "volume.hpp"
#include "warper.hpp"

#include <string>
#include <vector>

class Influence_manager
{
public:
    Influence_manager(const Parser& parser_, const Patient_Parameters_t& pat,
                      Warper_t warper_, const Volume_metadata_t ct_metadata_);
    ~Influence_manager();

    void calculate_at_plan();
    void calculate_at_fraction(std::vector<float> vec);
    void calculate(Array4<float>& influence, Array3<float>& positions,
                   std::vector<float> new_energies = std::vector<float>());

    Influence_engines_t engine;
    Array4<float> matrix_at_plan;
    Array4<float> matrix_at_fraction;

private:
    Influence_manager();
    void influence_from_beam_model(Array4<float>& infe, const std::vector<float>& ne = std::vector<float>());
    void influence_from_gpmc_dij();
    void influence_from_gpmc_dose();

    // beam_model and dose
    void structure_sampler();

    Array3<float> pos_at_plan;
    Array3<float> pos_at_fraction;
    Array3<float> boundary_at_plan;

    uint matrix_elements = 0;
    uint n_probing_positions = 0;
    uint nspots = 0;
    float mask_sampling_percentage = 25;
    bool get_boundary = false;

    const Patient_Parameters_t* patient_parameters;
    std::string outputdir;

    // Needed for structure sampling
    std::string ct_mask_file;
    const CT_Dims_t* ctdims;
    Warper_t warper;

    // Needed for beam model calculation
    const Volume_metadata_t ct_metadata;

    bool save_memory = true;
};

#endif