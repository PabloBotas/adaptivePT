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

    void get_dij_at_plan();
    void get_dij_at_frac(const std::vector<float>& new_energies = std::vector<float>());
    void get_dose_at_plan();
    void get_dose_at_frac(const std::vector<float>& new_energies = std::vector<float>());

    Influence_engines_t engine;
    std::vector<float> dose_at_plan;
    std::vector<float> dose_at_frac;
    Array4<float> matrix_at_plan;
    Array4<float> matrix_at_frac;
    uint n_voxels = 0;
    uint n_spots = 0;
private:
    Influence_manager();
    void influence_from_beam_model(std::string outfile, Array4<float>& infe,
                                   const std::vector<float>& ne = std::vector<float>());
    void influence_from_gpmc_dij();
    void influence_from_gpmc_dose();

    void get_dij(std::string outfile,
                 Array4<float>& influence, Array3<float>& positions,
                 const std::vector<float>& new_energies = std::vector<float>());
    void get_dose(std::string dose_file,
                  Array4<float>& matrix, std::vector<float>& dose,
                  const std::vector<float>& new_energies = std::vector<float>());
    void read_dose_file (std::string file, std::vector<float>& dose);

    void structure_sampler();

    std::string dose_plan_file;
    std::string dose_frac_file;
    std::string gpmc_dij_plan_file;
    std::string gpmc_dij_frac_file;
    std::string beam_model_dij_plan_file;
    std::string beam_model_dij_frac_file;

    Array3<float> pos_at_plan;
    Array3<float> pos_at_frac;
    Array3<float> boundary_at_plan;

    const Patient_Parameters_t* patient_parameters;
    const CT_Dims_t* ctdims;
    const Volume_metadata_t ct_metadata;
    std::string outputdir;

    std::string ct_mask_file;
    Warper_t warper;
    uint matrix_elements = 0;
    float mask_sampling_percentage = 25;
    bool get_boundary = false;
    bool save_memory = true;


};

#endif