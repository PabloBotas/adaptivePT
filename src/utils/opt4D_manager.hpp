#ifndef __OPT4D_MANAGER_HPP__
#define __OPT4D_MANAGER_HPP__

#include "enviroment.hpp"
#include "vector4.hpp"
#include "volume.hpp"

#include <string>
#include <valarray>
#include <vector>

class Opt4D_manager
{
public:
    Opt4D_manager(std::string outdir);
    ~Opt4D_manager();
    // void populate_directory(const uint& n_spots_, const uint& n_voxels_,
    //                         const Array4<float>& influence1,
    //                         const Array4<float>& influence2);
    void populate_directory(const uint& n_spots_, const uint& n_voxels_,
                            const uint& target_nvox_, const uint& rim_nvox_,
                            const uint& oars_nvox_, const std::vector<float>& mask,
                            const Volume_t& target_mask,
                            const Volume_t& rim_mask, const Volume_t& oars_mask,
                            const std::valarray<float>& target_dose,
                            const std::vector<std::valarray<float>>& adapt_field_dij);
    void launch_optimization();
    std::vector<float> get_weight_scaling();
private:
    Opt4D_manager();
    void read_bwf_file();
    void set_write_reference_dose(const std::vector<float>& dose);
    void set_write_reference_dose(const std::valarray<float>& dose);
    void set_write_reference_influence(const Array4<float>& influence);
    void write_templates();
    void write_dif();
    void write_vv(const std::vector<float>& mask, const Volume_t& target_mask,
                  const Volume_t& rim_mask, const Volume_t& oars_mask);
    void write_dij(const Array4<float>& influence);
    void write_dij(const std::vector<float>& data);
    void write_dij(const std::vector<std::valarray<float>>& data);

    uint n_spots;
    uint n_voxels;
    uint target_nvox;
    uint rim_nvox;
    uint oars_nvox;
    std::string out_directory;
    float min_average_constrain;
    std::vector<float> weight_scaling;

    // default file names
    std::string launcher_file;
    std::string launcher_file_base = "opt4D_launcher.sh";
    std::string template_launcher_file = std::string(INSTALLATION_PATH) +
                                         "/src/extra/opt4D_launcher_template.sh";
    std::string plan_file;
    std::string plan_file_base = "opt4D_planfile.pln";
    std::string template_plan_file = std::string(INSTALLATION_PATH) +
                                     "/src/extra/opt4D_planfile_template.pln";
    std::string dif_file;
    std::string dif_file_base = "dimensions.dif";
    std::string vv_file;
    std::string vv_file_base = "structures.vv";
    std::string dij_file;
    std::string dij_file_base = "beam_1.dij";
    std::string reference_file;
    std::string reference_file_base = "missing_target_dose.dat";
    std::string bwf_file;
    std::string bwf_file_base = "beam_1.bwf";

};

#endif