#include "opt4D_manager.hpp"

#include "utils.hpp"
#include "enviroment.hpp"
#include "vector4.hpp"

#include <fstream>
#include <map>
#include <limits>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <valarray>
#include <vector>


Opt4D_manager::Opt4D_manager(std::string outdir) : out_directory(outdir),
                                                   min_average_constrain(0)
{
    launcher_file  = outdir + '/' + launcher_file_base;
    plan_file      = outdir + '/' + plan_file_base;
    dif_file       = outdir + '/' + dif_file_base;
    vv_file        = outdir + '/' + vv_file_base;
    dij_file       = outdir + '/' + dij_file_base;
    reference_file = outdir + '/' + reference_file_base;
    bwf_file       = outdir + '/' + bwf_file_base;
}


Opt4D_manager::~Opt4D_manager()
{
    
}


void Opt4D_manager::launch_optimization()
{
    std::string cmd = "bash " + launcher_file +" "+ plan_file_base +" "+ out_directory;
    std::string std_out = utils::run_command(cmd);
    size_t pos = std_out.find_last_of('\n');
    pos = std_out.find_last_of('\n', pos-1);
    pos = std_out.find_last_of('\n', pos-1);
    pos = std_out.find_last_of('\n', pos-1);

    std::cout << std_out.substr(pos+1, std_out.size()-1) << std::endl;
}


// void Opt4D_manager::populate_directory(const uint& n_spots_, const uint& n_voxels_,
//                                        const Array4<float>& influence_ct,
//                                        const Array4<float>& influence_cbct)
// {
//     n_spots = n_spots_;
//     n_voxels = n_voxels_;

//     std::cout << "Writting Opt4D files:" << std::endl;
//     mkdir(out_directory.c_str(), 0774);
//     std::cout << "\t- " << dif_file_base << std::endl;
//     write_dif();
//     std::cout << "\t- " << vv_file_base << std::endl;
//     write_vv();
//     std::cout << "\t- " << dij_file_base << std::endl;
//     write_dij(influence_cbct);
//     std::cout << "\t- " << reference_file_base << std::endl;
//     set_write_reference_influence(influence_ct);
//     std::cout << "\t- templates" << std::endl;
//     write_templates();
// }


void Opt4D_manager::populate_directory(const uint& n_spots_, const uint& n_voxels_,
                                       const uint& target_nvox_, const uint& rim_nvox_,
                                       const uint& oars_nvox_, const std::vector<float>& mask,
                                       const Volume_t& target_mask,
                                       const Volume_t& rim_mask, const Volume_t& oars_mask,
                                       const std::valarray<float>& target_dose,
                                       const std::vector<std::valarray<float>>& adapt_field_dij)
{
    // TODO: The number of voxels must be computed without overlaps!!
    n_spots = n_spots_;
    n_voxels = n_voxels_;
    target_nvox = target_nvox_;
    rim_nvox = rim_nvox_;
    oars_nvox = oars_nvox_;

    std::cout << "Writting Opt4D files:" << std::endl;
    mkdir(out_directory.c_str(), 0774);
    std::cout << "\t- " << dif_file_base << std::endl;
    write_dif();
    std::cout << "\t- " << vv_file_base << std::endl;
    write_vv(mask, target_mask, rim_mask, oars_mask);
    std::cout << "\t- " << dij_file_base << std::endl;
    write_dij(adapt_field_dij);
    std::cout << "\t- " << reference_file_base << std::endl;
    set_write_reference_dose(target_dose);
    std::cout << "\t- templates" << std::endl;
    write_templates();
}


void Opt4D_manager::write_templates()
{
    if (min_average_constrain > 0) {
        std::map<std::string, std::string> replace_map {
            {"MINAVERAGECONS", std::to_string(min_average_constrain)},
        };
        utils::copy_replace_in_file(template_plan_file, plan_file, replace_map);
    } else {
        utils::copy_file(template_plan_file, plan_file);
    }
    // Bash script
    utils::copy_file(template_launcher_file, launcher_file);
}


void Opt4D_manager::write_dif() {
    std::ofstream stream (dif_file);
    if (stream.is_open()) {
        stream << "Delta-X 1" << std::endl;
        stream << "Delta-Y 1" << std::endl;
        stream << "Delta-Z 1" << std::endl;
        stream << "Dimension-CT-X 1" << std::endl;
        stream << "Dimension-CT-Y 1" << std::endl;
        stream << "Dimension-CT-Z 1" << std::endl;
        stream << "Dimension-Dose-X " << n_voxels << std::endl;
        stream << "Dimension-Dose-Y 1" << std::endl;
        stream << "Dimension-Dose-Z 1" << std::endl;
        stream << "ISO-Index-CT-X 0" << std::endl;
        stream << "ISO-Index-CT-Y 0" << std::endl;
        stream << "ISO-Index-CT-Z 0" << std::endl;
        stream << "ISO-Index-Dose-X 0" << std::endl;
        stream << "ISO-Index-Dose-Y 0" << std::endl;
        stream << "ISO-Index-Dose-Z 0" << std::endl;
        stream.close();
    } else {
        std::cerr << "Can't open " + dif_file + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Opt4D_manager::write_vv(const std::vector<float>& mask, const Volume_t& target_mask,
                             const Volume_t& rim_mask, const Volume_t& oars_mask)
{
    struct vvheader {
        char file_type_name[8];     // "vv-opt4D"
        short endian_test_pattern;  // 25185 = "ab" or \x6261 on little
                                    // endian systems
        char newline_character;     // "\n"
        char nVois_text[6];         // "#VOIs:"
    };
    vvheader header = {
        {'v','v','-','o','p','t','4','D'}, 25185, '\n',
        {'#','V','O','I','s',':'}};

    std::string strvv = vv_file;
    // Write header
    std::ofstream stream1 (strvv);
    if (stream1.is_open()) {
        stream1.write((char*)(&header), 17);
        char nVois_text1[6]  = {'0','0','0','0','3','\n'};
        char nVois_text2[13] = {'0','0','0','0','0',' ','t','a','r','g','e','t','\n'};
        char nVois_text3[14] = {'0','0','0','0','1',' ','f','a','l','l','o','f','f','\n'};
        char nVois_text4[11] = {'0','0','0','0','2',' ','o','a','r','s','\n'};
        stream1.write((char*) &nVois_text1, 6);
        stream1.write((char*) &nVois_text2, 13);
        stream1.write((char*) &nVois_text3, 14);
        stream1.write((char*) &nVois_text4, 11);
        stream1.close();
    } else {
        std::cerr << "Can't open " + strvv + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Write patient structures
    std::ofstream stream2 (strvv, std::ios::binary | std::ios::app);
    if (stream2.is_open()) {
        unsigned char dummy = 0;
        stream2.write((char*) &target_nvox, sizeof(uint));
        for (size_t i = 0, imask = 0, istruct = 0;
             i < mask.size() && istruct < target_nvox;
             ++i) {
            if (mask.at(i) != 0) {
                if (target_mask.data.at(i) != 0) {
                    stream2.write((char*) &imask, sizeof(uint));
                    stream2.write((char*) &dummy, sizeof(unsigned char));
                    istruct++;
                }
                imask++;
            }
        }
        stream2.write((char*) &rim_nvox, sizeof(uint));
        for (size_t i = 0, imask = 0, istruct = 0;
             i < mask.size() && istruct < rim_nvox;
             ++i) {
            if (mask.at(i) != 0) {
                if (rim_mask.data.at(i) != 0) {
                    stream2.write((char*) &imask, sizeof(uint));
                    stream2.write((char*) &dummy, sizeof(unsigned char));
                    istruct++;
                }
                imask++;
            }
        }
        stream2.write((char*) &oars_nvox, sizeof(uint));
        for (size_t i = 0, imask = 0, istruct = 0;
             i < mask.size() && istruct < oars_nvox;
             ++i) {
            if (mask.at(i) != 0) {
                if (oars_mask.data.at(i) != 0) {
                    stream2.write((char*) &imask, sizeof(uint));
                    stream2.write((char*) &dummy, sizeof(unsigned char));
                    istruct++;
                }
                imask++;
            }
        }
        stream2.close();
    } else {
        std::cerr << "Can't open " + strvv + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Opt4D_manager::write_dij(const Array4<float>& influence_cbct)
{
    std::vector<float> data(influence_cbct.size());
    for (size_t i = 0; i < influence_cbct.size(); ++i) {
        data.at(i) = influence_cbct.at(i).w;
    }
    write_dij(data);
}


void Opt4D_manager::write_dij(const std::vector<std::valarray<float>>& data)
{
    uint total_size = 0;
    for (const auto& field_dij : data)
        total_size += field_dij.size();
    std::vector<float> full_dij;
    full_dij.reserve(total_size);
    for (const auto& field_dij : data) {
#ifdef __GOOD_VALARRAYS
        full_dij.insert(full_dij.end(), field_dij.begin(), field_dij.end());
#else
        full_dij.insert(full_dij.end(), std::begin(field_dij), std::end(field_dij));
#endif
    }
    write_dij(full_dij);
}


void Opt4D_manager::write_dij(const std::vector<float>& data)
{
    // Get normalizing factor
    float m = 0;
    for (size_t i = 0; i < data.size(); ++i)
        m = data.at(i) < m ? m : data.at(i);
    float factor = m / std::numeric_limits<short>::max();

    // Get non-zero values and indexes
    std::vector<int> non_zeros(n_spots, 0);
    std::vector< std::vector<int>> indexes(n_spots);
    std::vector< std::vector<short>> values(n_spots);
    for (size_t i = 0; i<n_spots; i++) {
        for (size_t j = 0; j<n_voxels; j++) {
            int index = n_voxels*i+j;
            short value = (short)( data.at(index)/factor + 0.5 );
            if (value > 0) {
                indexes.at(i).push_back( (int)j );
                values.at(i).push_back( value );
                non_zeros.at(i) += 1;
            }
        }
    }

    // Write data
    float dummy_float_0 = 0;
    float dummy_float_1 = 1;
    float dummy_float = 0;
    int n_vox_int = n_voxels;
    int n_spots_int = n_spots;
    int dummy_int_1 = 1;
    std::ofstream stream (dij_file, std::ios::binary);
    if (stream.is_open()) {
        // gantry, table, collimator angles
        stream.write((char*) &dummy_float_0, sizeof(float));
        stream.write((char*) &dummy_float_0, sizeof(float));
        stream.write((char*) &dummy_float_0, sizeof(float));
        // bixel spacing
        stream.write((char*) &dummy_float_1, sizeof(float));
        stream.write((char*) &dummy_float_1, sizeof(float));
        // voxel spacing
        stream.write((char*) &dummy_float_1, sizeof(float));
        stream.write((char*) &dummy_float_1, sizeof(float));
        stream.write((char*) &dummy_float_1, sizeof(float));
        // voxel cube dimensions
        stream.write((char*) &n_vox_int, sizeof(int));
        stream.write((char*) &dummy_int_1, sizeof(int));
        stream.write((char*) &dummy_int_1, sizeof(int));
        // number of pencil beams
        stream.write((char*) &n_spots_int, sizeof(int));
        // dose scale factor
        stream.write((char*) &factor, sizeof(float));

        int totalNonZero=0;
        for (size_t i=0; i < n_spots; i++, dummy_float++) {
            stream.write((char*) &dummy_float, sizeof(float)); // Energy
            stream.write((char*) &dummy_float, sizeof(float)); // Spot-x
            stream.write((char*) &dummy_float, sizeof(float)); // Spot-y
            stream.write((char*) &non_zeros.at(i), sizeof(int)); // Non zeros (elements to read)
            for (int j=0; j < non_zeros.at(i); j++) {
                stream.write((char*) &indexes.at(i).at(j), sizeof(int));
                stream.write((char*) &values.at(i).at(j), sizeof(short));
            }
            totalNonZero += non_zeros.at(i);
        }
    } else {
        std::cerr << "Can't open " + dij_file + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Opt4D_manager::set_write_reference_dose(const std::valarray<float>& dose)
{
    float min_ref = dose.min();
    float max_ref = dose.max();
    float ave_ref = dose.sum()/dose.size();
    
    std::cout << "Reference in target for optimization (min, max, ave): ";
    std::cout << min_ref << ", " << max_ref << ", " << ave_ref << std::endl;
    min_average_constrain = ave_ref;

    // Write data
    std::ofstream stream (reference_file, std::ios::binary);
    if (stream.is_open()) {
        stream.write((char*) &dose[0], dose.size()*sizeof(float));
    } else {
        std::cerr << "Can't open " + reference_file + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Opt4D_manager::set_write_reference_dose(const std::vector<float>& dose)
{
    float min_ref = 1000000000;
    float max_ref = 0;
    uint nbatches = 100;
    std::vector<float> ave_ref_batches(nbatches, 0);
    float ave_ref = 0;
    for(size_t i = 0; i<n_voxels; i++) {
        min_ref = std::min(min_ref, dose.at(i));
        max_ref = std::max(max_ref, dose.at(i));
        ave_ref_batches.at(int(i%nbatches)) += dose.at(i)/n_spots;
    }
    for (uint i = 0; i < nbatches; ++i) {
        ave_ref += ave_ref_batches.at(i);
    }

    std::cout << "Reference in target for optimization (min, max, ave): ";
    std::cout << min_ref << ", " << max_ref << ", " << ave_ref << std::endl;
    min_average_constrain = ave_ref;

    // Write data
    std::ofstream stream (reference_file, std::ios::binary);
    if (stream.is_open()) {
        stream.write((char*) dose.data(),  dose.size()*sizeof(float));
    } else {
        std::cerr << "Can't open " + reference_file + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Opt4D_manager::set_write_reference_influence(const Array4<float>& influence)
{
    std::vector<float> dose(n_voxels, 0);
    // Accumulate influence on voxel j by all spots i(0 -> n_spots)
    for (size_t j = 0; j<n_voxels; j++)
        for(size_t i = 0; i<n_spots; i++)
            dose.at(j) += influence.at(n_voxels*i+j).w;
    set_write_reference_dose(dose);
}


void Opt4D_manager::read_bwf_file()
{
    std::ifstream infile( bwf_file );
    if (!infile.is_open()) {
        std::cout << "Can't open file: " << bwf_file << std::endl;
        exit(EXIT_FAILURE);
    }

    weight_scaling.resize(n_spots);

    std::string line;
    size_t i = 0;
    while (std::getline(infile, line))
    {
        if (line.at(0) == '#')
            continue;
        float w_;
        float dummy;
        std::istringstream iss(line);
        iss >> dummy >> w_ >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
        weight_scaling.at(i) = w_;
        i++;
    }
}


std::vector<float> Opt4D_manager::get_weight_scaling()
{
    if (weight_scaling.empty())
        read_bwf_file();
    return weight_scaling;
}

