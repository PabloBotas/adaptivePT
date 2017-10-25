#include "opt4D_manager.hpp"

#include "utils.hpp"
#include "enviroment.hpp"
#include "vector4.hpp"

#include <fstream>
#include <limits>
#include <string>
#include <sys/stat.h>


Opt4D_manager::Opt4D_manager(std::string outdir) : out_directory(outdir)
{
}


Opt4D_manager::~Opt4D_manager()
{
    
}


void Opt4D_manager::launch_optimization()
{
    std::string cmd = "bash " + out_directory + "/opt4D_launcher.sh " +
                      out_directory + "/opt4d_planfile.pln" + out_directory;
    utils::run_command(cmd);
}


void Opt4D_manager::populate_directory(const Array4<double>& influence_ct,
                                       const Array4<double>& influence_cbct)
{
    std::cout << "Writting Opt4D files ..." << std::endl;
    n = (unsigned int)(sqrt(influence_ct.size() + 0.5));
    write_templates();
    write_dif();
    write_vv();
    write_dij(influence_cbct);
    write_reference_influence(influence_ct);
}


void Opt4D_manager::write_templates()
{
    mkdir(out_directory.c_str(), 0774);

    std::ifstream src1(std::string(INSTALLATION_PATH) +
                       "/src/extra/opt4D_planfile_template.pln");
    if(!src1)
    {
        std::cerr << "Can't open " +
                     std::string(INSTALLATION_PATH) +
                     "/src/extra/opt4D_planfile_template.pln" +
                     " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ofstream dst1(out_directory+"/opt4D_planfile.pln");
    dst1 << src1.rdbuf();
    std::ifstream src2(std::string(INSTALLATION_PATH) +
                       "/src/extra/opt4D_launcher_template.sh");
    if(!src2)
    {
        std::cerr << "Can't open " +
                     std::string(INSTALLATION_PATH) +
                     "/src/extra/opt4D_planfile_template.pln" +
                     " file" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ofstream dst2(out_directory+"/opt4D_launcher.sh");
    dst2 << src2.rdbuf();
}


void Opt4D_manager::write_dif()
{
    std::ofstream stream (out_directory+"/dimensions.dif");
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

    std::string strvv = out_directory+"/structures.vv";
    // Write header
    std::ofstream stream1 (strvv);
    if (stream1.is_open()) {
        stream1.write((char*)(&header), 17);
        char nVois_text1[6]  = {'0','0','0','0','1','\n'};
        char nVois_text2[14] = {'0','0','0','0','0',' ','p','a','t','i','e','n','t','\n'};
        stream1.write((char*) &nVois_text1, 6);
        stream1.write((char*) &nVois_text2, 14);
        stream1.close();
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
    }
}

//-o /opt/utils/adaptive/test/water_patient/adaptation_x10_y20_z30/opt4D_reoptimization --lbfgs --lbfgs_m 20 --linesearch_steps 20 --dont_project_on_bounds --unit_initial_fluence --random_seed 0 --add_step_event "UPDATE_LAGRANGE_MULTIPLIERS,20,20" --add_step_event "UPDATE_PENALTY,20,20" --constraint_penalty_multiplier 2 --max_steps 2000 --min_steps 600 --max_time 1200 --write_dose --write_beam_dose opt4D_planfile.pln
void Opt4D_manager::write_dij(const Array4<double>& influence_cbct)
{
    // Get normalizing factor
    double m = 0;
    for (size_t i = 0; i < influence_cbct.size(); ++i)
        m = influence_cbct.at(i).w < m ? m : influence_cbct.at(i).w;
    float factor = m / std::numeric_limits<short>::max();

    // Get non-zero values and indexes
    std::vector<int> non_zeros(n, 0);
    std::vector< std::vector<int> > indexes(n);
    std::vector< std::vector<short> > values(n);
    for(size_t i = 0; i<n; i++) {
        for (size_t j = 0; j<n; j++) {
            int index = n*i+j;
            short value = (short)( influence_cbct.at(index).w/factor + 0.5 );
            if(value > 0) {
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
    std::ofstream stream (out_directory+"/beam_1.dij", std::ios::binary);
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
        for(size_t i=0; i < n; i++)
        {
            float energy = influence_cbct.at(i*n).x;
            float i_float = i;
            stream.write((char*) &energy, sizeof(float));
            stream.write((char*) &i_float, sizeof(float));
            stream.write((char*) &dummy_float_0, sizeof(float));
            stream.write((char*) &non_zeros.at(i), sizeof(int));

            for(int j=0; j < non_zeros.at(i); j++)
            {
                stream.write((char*) &indexes.at(i).at(j), sizeof(int));
                stream.write((char*) &values.at(i).at(j), sizeof(short));
            }
            totalNonZero += non_zeros.at(i);
        }
    }
}


void Opt4D_manager::write_reference_influence(const Array4<double>& influence)
{
    // Accumulate influence on voxel j by all spots i(0 -> n)
    std::vector<float> reference(n, 0);
    for(size_t j = 0; j<n; j++)
        for(size_t i = 0; i<n; i++)
            reference.at(j) += influence.at(n*i+j).w;

    // Write data
    std::ofstream stream (out_directory+"/reference_influence.dat", std::ios::binary);
    if (stream.is_open()) {
        stream.write((char*) reference.data(), n*sizeof(float));
    }
}



