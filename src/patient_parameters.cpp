#include "patient_parameters.hpp"
#include "patient_parameters_parser.hpp"

#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <sstream>
#include <math.h>

Patient_Parameters_t::Patient_Parameters_t(std::string dir) : patient_dir(dir),
                                                        input_dir(dir+"/input"),
                                                        tramp_dir(input_dir+"/tramps")
{
    exploreBeamDirectories();
    parseTopasFiles();
}

void Patient_Parameters_t::exploreBeamDirectories()
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

void Patient_Parameters_t::parseTopasFiles()
{
    getTopasGlobalParameters();
    getTopasBeamParameters();
}

void Patient_Parameters_t::getTopasGlobalParameters()
{
    Patient_Parameters_Parser_t pars(topas_files.front());

    // Machine
    machine = pars.readString("Rt/beam/TreatmentMachineName");
    virtualSAD.a = pars.readReal<float>("Rt/beam/VirtualSourceAxisDistances0");
    virtualSAD.b = pars.readReal<float>("Rt/beam/VirtualSourceAxisDistances1");

    // CT grid shift
    float ImgCenterX = pars.readReal<float>("Rt/CT/ImgCenterX");
    float ImgCenterY = pars.readReal<float>("Rt/CT/ImgCenterY");
    float ImgCenterZ = pars.readReal<float>("Rt/CT/ImgCenterZ");
    float IsoCenter0 = pars.readReal<float>("Rt/beam/IsoCenter0");
    float IsoCenter1 = pars.readReal<float>("Rt/beam/IsoCenter1");
    float IsoCenter2 = pars.readReal<float>("Rt/beam/IsoCenter2");
    ct.offset.x = ImgCenterX - IsoCenter0;
    ct.offset.y = ImgCenterY - IsoCenter1;
    ct.offset.z = ImgCenterZ - IsoCenter2;

    // CT grid resolution
    ct.d.x = pars.readReal<float>("Rt/CT/PixelSpacing0");
    ct.d.y = pars.readReal<float>("Rt/CT/PixelSpacing1");
    ct.d.z = pars.readVector<float>("Rt/CT/SliceThicknessSpacing", true);

    // CT grid number of voxels
    ct.n.x = pars.readReal<unsigned int>("Rt/CT/Columns");
    ct.n.y = pars.readReal<unsigned int>("Rt/CT/Rows");
    ct.n.z = pars.readVectorInts<unsigned int>("Rt/CT/SliceThicknessSections", true);

    ct.total = ct.n.x*ct.n.y*ct.n.z.front();
}

void Patient_Parameters_t::getTopasBeamParameters()
{
    // Resize containers
    apertures.resize(nbeams);
    range_shifters.resize(nbeams);
    angles.resize(nbeams);
    isocenter_to_beam_distance.resize(nbeams);

    for(size_t i = 0; i < nbeams; i++)
    {
        // Beam data from MCAUTO_DICOM.txt
        Patient_Parameters_Parser_t pars(topas_files.at(i));

        // Create shortcuts
        Aperture_Dims_t&     ap        = apertures.at(i);
        RangeShifter_Dims_t& rs        = range_shifters.at(i);
        BeamAngles_t&        angle     = angles.at(i);
        float&              isoToBeam = isocenter_to_beam_distance.at(i);

        // Aperture
        ap.exists = pars.readBool("Rt/beam/IncludeAperture", false);
        if(ap.exists)
        {
            ap.thick =  pars.readReal<float>("Rt/beam/BlockThickness", 0);
            ap.zdown = -pars.readReal<float>("Rt/beam/IsocenterToBlockTrayDistance", 0);
        }

        // Range shifter
        rs.exists = pars.readBool("Rt/beam/IncludeRangeShifter", false);
        if(rs.exists)
        {
            rs.thick =  pars.readReal<float>("Rt/beam/RangeShifterThickness", 0);
            rs.zdown = -pars.readReal<float>("Rt/beam/IsocenterToRangeShifterTrayDistance", 0);
            rs.zup   =  rs.zdown - rs.thick;
        }

        // Angles
        angle.gantry = (270.0 * (M_PI/180.0)) - pars.readReal<float>("Rt/beam/Gantry");
        angle.couch  = - pars.readReal<float>("Rt/beam/PatientSupportAngle");

        // Distance from isocenter to phase space plane
        // Downstream edge of range shifter - gets moved upstream later
        if (rs.zdown != 0)
            isoToBeam = rs.zdown;
        else
            isoToBeam = - pars.readReal<float>("Rt/beam/IsocenterToBeamDistance");
    }
}


void Patient_Parameters_t::print()
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

