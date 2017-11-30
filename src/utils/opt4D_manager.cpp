#include "opt4D_manager.hpp"

#include "utils.hpp"
#include "enviroment.hpp"
#include "vector4.hpp"

#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <sys/stat.h>


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

    std::cout << std_out.substr(pos+1, std_out.size()-1) << std::endl;
}


void Opt4D_manager::populate_directory(const Array4<float>& influence_ct,
                                       const Array4<float>& influence_cbct)
{
    std::cout << "Writting Opt4D files:" << std::endl;
    mkdir(out_directory.c_str(), 0774);
    n = (unsigned int)(sqrt(influence_ct.size() + 0.5));
    std::cout << "\t- " << dif_file_base << std::endl;
    write_dif();
    std::cout << "\t- " << vv_file_base << std::endl;
    write_vv();
    std::cout << "\t- " << dij_file_base << std::endl;
    write_dij(influence_cbct);
    std::cout << "\t- " << reference_file_base << std::endl;
    set_write_reference_influence(influence_ct);
    std::cout << "\t- templates" << std::endl;
    write_templates();
}


void Opt4D_manager::write_templates()
{
    std::ifstream src1(template_plan_file);
    if (!src1) {
        std::cerr << "Can't open " + template_plan_file + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ofstream dst1(plan_file);
    if (min_average_constrain > 0) {
        std::string planfile_contents;
        for (char ch; src1.get(ch); planfile_contents.push_back(ch)) {
        }
        std::string to_replace = "MINAVERAGECONS";

        size_t pos = planfile_contents.find(to_replace);
        while (pos != std::string::npos) {
            planfile_contents.replace(pos, to_replace.length(),
                                  std::to_string(min_average_constrain));
            pos = planfile_contents.find(to_replace, pos);
        }
        
        dst1 << planfile_contents;
    } else {
        dst1 << src1.rdbuf();
    }


    // Bash script
    std::ifstream src2(template_launcher_file);
    if (!src2) {
        std::cerr << "Can't open " + template_launcher_file +
                     " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ofstream dst2(launcher_file);
    dst2 << src2.rdbuf();
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
        stream << "Dimension-Dose-X " << n << std::endl;;
        stream << "Dimension-Dose-Y 1" << std::endl;;
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


void Opt4D_manager::write_vv()
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
        char nVois_text1[6]  = {'0','0','0','0','1','\n'};
        char nVois_text2[14] = {'0','0','0','0','0',' ','p','a','t','i','e','n','t','\n'};
        stream1.write((char*) &nVois_text1, 6);
        stream1.write((char*) &nVois_text2, 14);
        stream1.close();
    } else {
        std::cerr << "Can't open " + strvv + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Write first dummy patient structure
    std::ofstream stream2 (strvv, std::ios::binary | std::ios::app);
    if (stream2.is_open()) {
        unsigned char dummy = 0;
        stream2.write((char*) &n, sizeof(unsigned int));
        for (unsigned int i = 0; i < n; ++i) {
            stream2.write((char*) &i, sizeof(unsigned int));
            stream2.write((char*) &dummy, sizeof(unsigned char));
        }
        stream2.close();
    } else {
        std::cerr << "Can't open " + strvv + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//-o /opt/utils/adaptive/test/water_patient/adaptation_x10_y20_z30/opt4D_reoptimization --initial_fluence_noise 0.5 --lbfgs --lbfgs_m 20 --linesearch_steps 20 --dont_project_on_bounds --unit_initial_fluence --random_seed 0 --add_step_event "UPDATE_LAGRANGE_MULTIPLIERS,20,20" --add_step_event "UPDATE_PENALTY,20,20" --constraint_penalty_multiplier 2 --max_steps 2000 --min_steps 600 --max_time 1200 --write_dose --write_beam_dose opt4D_planfile.pln
void Opt4D_manager::write_dij(const Array4<float>& influence_cbct)
{
    // Get normalizing factor
    float m = 0;
    for (size_t i = 0; i < influence_cbct.size(); ++i)
        m = influence_cbct.at(i).w < m ? m : influence_cbct.at(i).w;
    float factor = m / std::numeric_limits<short>::max();

    // Get non-zero values and indexes
    std::vector<int> non_zeros(n, 0);
    std::vector< std::vector<int> > indexes(n);
    std::vector< std::vector<short> > values(n);
    for (size_t i = 0; i<n; i++) {
        for (size_t j = 0; j<n; j++) {
            int index = n*i+j;
            short value = (short)( influence_cbct.at(index).w/factor + 0.5 );
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
    int n_int = n;
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
        stream.write((char*) &n_int, sizeof(int));
        stream.write((char*) &dummy_int_1, sizeof(int));
        stream.write((char*) &dummy_int_1, sizeof(int));
        // number of pencil beams
        stream.write((char*) &n_int, sizeof(int));
        // dose scale factor
        stream.write((char*) &factor, sizeof(float));

        int totalNonZero=0;
        for (size_t i=0; i < n; i++) {
            float i_float = i;
            stream.write((char*) &i_float, sizeof(float));
            stream.write((char*) &i_float, sizeof(float));
            stream.write((char*) &dummy_float_0, sizeof(float));
            stream.write((char*) &non_zeros.at(i), sizeof(int));

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


void Opt4D_manager::set_write_reference_influence(const Array4<float>& influence)
{
    // Accumulate influence on voxel j by all spots i(0 -> n)
    std::vector<float> reference(n, 0);
    float min_ref = 1000000000;
    float max_ref = 0;
    float ave_ref = 0;
    for (size_t j = 0; j<n; j++) {
        for(size_t i = 0; i<n; i++)
            reference.at(j) += influence.at(n*i+j).w;
        if (min_ref > reference.at(j))
            min_ref = reference.at(j);
        if (max_ref < reference.at(j))
            max_ref = reference.at(j);
        ave_ref += reference.at(j)/n;
    }
    std::cout << "Reference influence for optimization (min, max, ave): ";
    std::cout << min_ref << ", " << max_ref << ", " << ave_ref << std::endl;
    min_average_constrain = ave_ref;

    // Write data
    std::ofstream stream (reference_file, std::ios::binary);
    if (stream.is_open()) {
        stream.write((char*) reference.data(), n*sizeof(float));
    } else {
        std::cerr << "Can't open " + reference_file + " file" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Opt4D_manager::read_bwf_file()
{
    std::ifstream infile( bwf_file );
    if (!infile.is_open()) {
        std::cout << "Can't open file: " << bwf_file << std::endl;
        exit(EXIT_FAILURE);
    }

    weight_scaling.resize(n);

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

