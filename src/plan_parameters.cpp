#include "plan_parameters.hpp"
#include "plan_parameters_parser.hpp"

#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <sstream>
#include <math.h>

Plan_Parameters_t::Plan_Parameters_t(std::string dir) : patient_dir(dir),
                                                        input_dir(dir+"/input"),
                                                        tramp_dir(input_dir+"/tramps")
{
    exploreBeamDirectories();
    parseTopasFiles();
}

void Plan_Parameters_t::exploreBeamDirectories()
{
    beam_dirs = getFoldersWithFile(input_dir, "run/MCAUTO_DICOM.txt");
    nbeams    = beam_dirs.size();
    
    // Get beam names
    beam_names.resize(nbeams);
    for(unsigned int i=0; i<nbeams; i++)
    {
        size_t pos = beam_dirs.at(i).rfind("/");
        if (pos != std::string::npos)
        {
            if (pos != beam_dirs.at(i).size() - 1)
                beam_names.at(i) = beam_dirs.at(i).substr(pos + 1);
            else
                beam_names.at(i) = beam_dirs.at(i);
        }
        else
        {
            beam_names.push_back(beam_dirs.at(i));
        }
    }

    // Set run directories
    run_dir.resize(nbeams);
    for (unsigned int i = 0; i < nbeams; i++)
    {
        run_dir[i] = beam_dirs[i];
        run_dir[i].append("/run");
    }

    // Get tramp files
    tramp_files.resize(nbeams);
    tramp_files = getFilesWithSuffix(tramp_dir, ".tramp");

    // Get topas parameter files
    topas_files.resize(nbeams);
    for(unsigned int i = 0; i < nbeams; i++)
    {
        topas_files.at(i) = getFilesWithSuffix(beam_dirs[i], "MCAUTO_DICOM.txt")[0];
    }
}

void Plan_Parameters_t::parseTopasFiles()
{
    getTopasGlobalParameters();
    getTopasBeamParameters();
}

void Plan_Parameters_t::getTopasGlobalParameters()
{
    Plan_Parameters_Parser_t pars(topas_files.front());

    // CT grid shift
    ct.offset.x = pars.readReal<double>("Rt/CT/ImgCenterX") - pars.readReal<double>("Rt/beam/IsoCenter0");
    ct.offset.y = pars.readReal<double>("Rt/CT/ImgCenterY") - pars.readReal<double>("Rt/beam/IsoCenter1");
    ct.offset.z = pars.readReal<double>("Rt/CT/ImgCenterZ") - pars.readReal<double>("Rt/beam/IsoCenter2");

    // CT grid resolution
    ct.d.x = pars.readReal<double>("Rt/CT/PixelSpacing0");
    ct.d.y = pars.readReal<double>("Rt/CT/PixelSpacing1");
    ct.d.z = pars.readVector<double>("Rt/CT/SliceThicknessSpacing", true);

    // CT grid number of voxels
    ct.n.x = pars.readReal<double>("Rt/CT/Columns");
    ct.n.y = pars.readReal<double>("Rt/CT/Rows");
    ct.n.z = pars.readVectorInts<unsigned int>("Rt/CT/SliceThicknessSections", true);
}

void Plan_Parameters_t::getTopasBeamParameters()
{
    // Resize containers
    apertures.resize(nbeams);
    range_shifters.resize(nbeams);
    angles.resize(nbeams);
    isocenter_to_beam_distance.resize(nbeams);

    for(size_t i = 0; i < nbeams; i++)
    {
        // Beam data from MCAUTO_DICOM.txt
        Plan_Parameters_Parser_t pars(topas_files.at(i));

        // Create shortcuts
        Aperture_Dims_t&     ap        = apertures.at(i);
        RangeShifter_Dims_t& rs        = range_shifters.at(i);
        BeamAngles_t&        angle     = angles.at(i);
        double&              isoToBeam = isocenter_to_beam_distance.at(i);

        // Aperture
        ap.exists = pars.readBool("Rt/beam/IncludeAperture", false);
        if(ap.exists)
        {
            ap.thick =  pars.readReal<double>("Rt/beam/BlockThickness", 0);
            ap.zdown = -pars.readReal<double>("Rt/beam/IsocenterToBlockTrayDistance", 0);
        }

        // Range shifter
        rs.exists = pars.readBool("Rt/beam/IncludeRangeShifter", false);
        if(rs.exists)
        {
            rs.thick =  pars.readReal<double>("Rt/beam/RangeShifterThickness", 0);
            rs.zdown = -pars.readReal<double>("Rt/beam/IsocenterToRangeShifterTrayDistance", 0);
            rs.zup   =  rs.zdown - rs.thick;
        }

        // Angles
        angle.gantry = (270.0 * (M_PI/180.0)) - pars.readReal<double>("Rt/beam/Gantry");
        angle.couch  = pars.readReal<double>("Rt/beam/PatientSupportAngle");

        // Distance from isocenter to phase space plane
        // Downstream edge of range shifter - gets moved upstream later
        if (rs.zdown != 0)
            isoToBeam = rs.zdown;
        else
            isoToBeam = - pars.readReal<double>("Rt/beam/IsocenterToBeamDistance");
    }
}


void Plan_Parameters_t::print()
{
    std::cout << "Directories:" << std::endl;
    std::cout << "    - Patient:     " << patient_dir << std::endl;
    std::cout << "    - Beam parent: " << input_dir << std::endl;
    std::cout << "    - Beams:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
        std::cout << "        " << beam_dirs.at(i) << std::endl;
    
    std::cout << "Files:" << std::endl;
    std::cout << "    - Tramp files:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
        std::cout << "        " << tramp_files.at(i) << std::endl;
    std::cout << "    - Topas files:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
        std::cout << "        " << topas_files.at(i) << std::endl;

    std::cout << "Beam names:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
        std::cout << "        " << beam_names.at(i) << std::endl;

    std::cout << "CT data:" << std::endl;
    std::cout << "    - Voxels:     " << ct.n.x << ", " << ct.n.y << ", (";
    for (size_t i = 0; i < ct.n.z.size(); i++)
    {
        std::cout << ct.n.z.at(i);
        if (i != ct.n.z.size()-1)
            std::cout << ", ";
        else
            std::cout << ")" << std::endl;
    }
    std::cout << "    - Voxel size: " << ct.d.x << ", " << ct.d.y << ", (";
    for (size_t i = 0; i < ct.d.z.size(); i++)
    {
        std::cout << ct.d.z.at(i);
        if (i != ct.d.z.size()-1)
            std::cout << ", ";
        else
            std::cout << ") cm" << std::endl;
    }
    std::cout << "    - Offset:     " << ct.offset.x << ", " << ct.offset.y << ", " << ct.offset.z << " cm" << std::endl;

    std::cout << "Beam angles:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
    {
        std::cout << "    - " << beam_names.at(i) << ":" << std::endl;
        std::cout << "        - gantry: " << angles.at(i).gantry*180.0/M_PI << " deg" << std::endl;
        std::cout << "        - couch:  " << angles.at(i).couch*180.0/M_PI  << " deg" << std::endl;
    }
    std::cout << "Isocenter to distance:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
    {
        std::cout << "    - " << beam_names.at(i) << ": " << isocenter_to_beam_distance.at(i) << " cm" << std::endl;
    }
}


std::vector<std::string> getFoldersWithFile(std::string folderpath, std::string name)
//      return all folder names of all folders
{
    std::vector<std::string> returnStrings; //empty returnvalue

    struct dirent *direntp = NULL;
    DIR *dirp = NULL;

    dirp = opendir(folderpath.c_str());
    if (dirp == NULL)
    {
        std::cerr << "ERROR, " << folderpath << " cannot be openned!" << std::endl;
        std::cerr << "Was searching for " << name  << std::endl;
        perror("Cannot open directory");
        exit(EXIT_FAILURE);
    }

    // For every directory entry...
    while ((direntp = readdir(dirp)) != NULL)
    {
        // Ignore special directories.
        if ((std::strcmp(direntp->d_name, ".") == 0) ||
                (std::strcmp(direntp->d_name, "..") == 0))
            continue;

        struct stat fstat;
        std::string full_name;

        full_name = std::string(folderpath);
        if(*(full_name.end()-1) != '/')
            full_name += std::string("/");
        full_name += std::string(direntp->d_name);

        // Get only if it is really directory.
        if (stat(full_name.c_str(), &fstat) < 0)
            continue;
        if (S_ISDIR(fstat.st_mode))
        {
            // Check name
            if (stat((full_name + "/" + name).c_str() , &fstat) < 0)
                continue;
            returnStrings.push_back(full_name);
        }
    }

    return returnStrings;
}

std::vector<std::string> getFilesWithSuffix(std::string folderpath, std::string suffix, std::string contains)
//  get a list of files inside of the folderpath with given suffix
{
    std::vector<std::string> returnStrings;

    struct dirent *direntp = NULL;
    DIR *dirp = NULL;

    dirp = opendir(folderpath.c_str());
    if (dirp == NULL)
    {
        std::cerr << "ERROR, " << folderpath << " cannot be openned!" << std::endl;
        std::cerr << "Was searching for " << suffix  << std::endl;
        perror("Cannot open directory");
        exit(EXIT_FAILURE);
    }

    // For every directory entry...
    while ((direntp = readdir(dirp)) != NULL)
    {
        // Ignore special directories.
        if ((std::strcmp(direntp->d_name, ".") == 0) || (std::strcmp(direntp->d_name, "..") == 0))
            continue;

        struct stat fstat;
        std::string full_name;

        //  get full name
        full_name = std::string(folderpath);
        if (*(full_name.end()-1) != '/')
            full_name += std::string("/");
        full_name += std::string(direntp->d_name);

        //  get info for this entry(folder or file)
        if (stat(full_name.c_str(), &fstat) < 0)
            continue;

        //  if a file, check it
        if(!S_ISDIR(fstat.st_mode))
        {
            std::string temp(direntp->d_name);
            int n = temp.size();
            int start = n-suffix.size();
            std::string MacFilesPref = "._";
            if(start<0) continue;  //the current file name is shorter than suffix
            if(suffix == temp.substr(n-suffix.size(),n) && temp.substr(0,2) != MacFilesPref) {
                // Check name
                if (contains.empty() || temp.find(contains) != std::string::npos) {
                    returnStrings.push_back(full_name);
                }
            }
        }
        //  if a folder, recursive
        else
        {
            std::vector<std::string> temp = getFilesWithSuffix(full_name, suffix);
            returnStrings.insert(returnStrings.end(), temp.begin(), temp.end() ); //append
        }
    }

    return returnStrings;
}
