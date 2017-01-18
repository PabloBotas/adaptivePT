#ifndef __PLAN_PARAMETERS_HPP__
#define __PLAN_PARAMETERS_HPP__

#include "plan_parameters_parser.hpp"

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

template<class T=double, class D=T>
struct Vector_t
{
    T x;
    T y;
    D z;
};

struct BeamAngles_t
{
    double gantry;
    double couch;
};

struct CT_Dims_t
{
    Vector_t<> offset;
    Vector_t<double, std::vector<double> > d;
    Vector_t<unsigned int, std::vector<unsigned int> > n;
};

struct Aperture_Dims_t
{
    bool   exists;
    double thick;
    double zdown;
    Aperture_Dims_t() : exists(false), thick(0), zdown(0) {}
};

struct RangeShifter_Dims_t
{
    bool   exists;
    double thick;
    double zdown;
    double zup;
    RangeShifter_Dims_t() : exists(false), thick(0), zdown(0), zup(0) {}
};

std::vector<std::string> getFilesWithSuffix(std::string folderpath,
                                            std::string suffix,
                                            std::string contains="");
std::vector<std::string> getFoldersWithFile(std::string folderpath, std::string name);

class Plan_Parameters_t
{
public:
    Plan_Parameters_t(std::string directory);
    void print();
    
    std::string patient_dir;
    std::string input_dir;
    std::string tramp_dir;

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
    std::vector<double>       isocenter_to_beam_distance;

private:
    void exploreBeamDirectories();
    void parseTopasFiles();
    void getTopasGlobalParameters();
    void getTopasBeamParameters();

//     std::string rel_beam_location  = "/input";
//     std::string rel_topas_location = "/run";
};

#endif