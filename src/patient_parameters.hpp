#ifndef __PATIENT_PARAMETERS_HPP__
#define __PATIENT_PARAMETERS_HPP__

#include "patient_parameters_parser.hpp"
#include "special_types.hpp"

#include <string>
#include <vector>

/*
Patient dir:
    - input:
        * tramps
        * ct
        * beam1:
            - run:
                * MCAUTO_DICOM.txt
                * ...
            - ...
        * ...
    -...
*/



class Patient_Parameters_t
{
public:
    Patient_Parameters_t(std::string directory);
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
    std::vector<std::string> topas_files;

    CT_Dims_t ct;
    std::vector<Aperture_Dims_t>     apertures;
    std::vector<RangeShifter_Dims_t> range_shifters;

    std::vector<BeamAngles_t> angles;
    std::vector<float>        isocenter_to_beam_distance;

    std::string results_dir;

    void add_results_directory(std::string s);

private:
    void exploreBeamDirectories();
    void parseTopasFiles();
    void getTopasGlobalParameters();
    void getTopasBeamParameters();

//     std::string rel_beam_location  = "/input";
//     std::string rel_topas_location = "/run";
};

#endif