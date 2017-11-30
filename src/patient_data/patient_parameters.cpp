#include "patient_parameters.hpp"
#include "patient_parameters_parser.hpp"
#include "tramp.hpp"
#include "special_types.hpp"
#include "utils.hpp"

#include <iostream>
#include <string>
#include <cmath>
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
    set_spots_data();
    set_planning_CT_file();
}

void Patient_Parameters_t::set_planning_CT_file()
{
    planning_ct_file = patient_dir + "/input/ctbinary/ctvolume.dat";
}

void Patient_Parameters_t::exploreBeamDirectories()
{
    beam_dirs = getFoldersWithFile(input_dir, "run/MCAUTO_DICOM.txt");
    nbeams    = beam_dirs.size();
    
    // Get beam names
    beam_names.resize(nbeams);
    for (unsigned int i=0; i<nbeams; i++) {
        size_t pos = beam_dirs.at(i).rfind("/");
        if (pos != std::string::npos) {
            if (pos != beam_dirs.at(i).size() - 1)
                beam_names.at(i) = beam_dirs.at(i).substr(pos + 1);
            else
                beam_names.at(i) = beam_dirs.at(i);
        } else {
            beam_names.push_back(beam_dirs.at(i));
        }
    }

    // Set run directories
    run_dir.resize(nbeams);
    for (unsigned int i = 0; i < nbeams; i++) {
        run_dir[i] = beam_dirs[i];
        run_dir[i].append("/run");
    }

    // Get tramp files
    tramp_files.resize(nbeams);
    tramp_files = getFilesWithSuffix(tramp_dir, ".tramp");

    // Get topas parameter files
    topas_files.resize(nbeams);
    for (unsigned int i = 0; i < nbeams; i++) {
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

    // CT grid resolution
    ct.d.x = pars.readReal<float>("Rt/CT/PixelSpacing0");
    ct.d.y = pars.readReal<float>("Rt/CT/PixelSpacing1");
    ct.d.z = pars.readLastRealInVector<float>("Rt/CT/SliceThicknessSpacing", true);

    // CT grid number of voxels
    ct.n.x = pars.readInteger<unsigned int>("Rt/CT/Columns");
    ct.n.y = pars.readInteger<unsigned int>("Rt/CT/Rows");
    ct.n.z = pars.readLastIntInVector<unsigned int>("Rt/CT/SliceThicknessSections", true);

    ct.total = ct.n.x*ct.n.y*ct.n.z;

    // CT grid shift
    float ImgCenterX = pars.readReal<float>("Rt/CT/ImgCenterX");
    float ImgCenterY = pars.readReal<float>("Rt/CT/ImgCenterY");
    float ImgCenterZ = pars.readReal<float>("Rt/CT/ImgCenterZ");
    ct.isocenter.x = pars.readReal<float>("Rt/beam/IsoCenter0");
    ct.isocenter.y = pars.readReal<float>("Rt/beam/IsoCenter1");
    ct.isocenter.z = pars.readReal<float>("Rt/beam/IsoCenter2");

    // Offset of first voxel corner to isocenter
    ct.offset.x = ImgCenterX - ct.isocenter.x;
    ct.offset.y = ImgCenterY - ct.isocenter.y;
    ct.offset.z = ImgCenterZ - ct.isocenter.z;
}

void Patient_Parameters_t::consolidate_originals()
{
    original_ct = ct;
    original_angles = angles;
    original_isocenter_to_beam_distance = isocenter_to_beam_distance;
}

void Patient_Parameters_t::getTopasBeamParameters()
{
    // Resize containers
    apertures.resize(nbeams);
    range_shifters.resize(nbeams);
    angles.resize(nbeams);
    isocenter_to_beam_distance.resize(nbeams);

    for (size_t i = 0; i < nbeams; i++) {
        // Beam data from MCAUTO_DICOM.txt
        Patient_Parameters_Parser_t pars(topas_files.at(i));

        // Create shortcuts
        Aperture_Dims_t& ap = apertures.at(i);
        RangeShifter_Dims_t& rs = range_shifters.at(i);
        BeamAngles_t& angle = angles.at(i);
        float& isoToBeam = isocenter_to_beam_distance.at(i);

        // Aperture
        ap.exists = pars.readBool("Rt/beam/IncludeAperture", false);
        if (ap.exists) {
            ap.thick =  pars.readReal<float>("Rt/beam/BlockThickness", 0);
            ap.zdown = -pars.readReal<float>("Rt/beam/IsocenterToBlockTrayDistance", 0);
        }

        // Range shifter
        rs.exists = pars.readBool("Rt/beam/IncludeRangeShifter", false);
        if (rs.exists) {
            rs.thick = pars.readReal<float>("Rt/beam/RangeShifterThickness", 0);
            rs.zdown = pars.readReal<float>("Rt/beam/IsocenterToRangeShifterTrayDistance", 0);
        }

        // Angles
        angle.gantry = pars.readReal<float>("Rt/beam/Gantry");
        angle.couch  = pars.readReal<float>("Rt/beam/PatientSupportAngle");

        // Distance from isocenter to phase space plane
        // Downstream edge of range shifter - gets moved upstream later
        if (rs.thick != 0)
            isoToBeam = rs.zdown;
        else
            isoToBeam = pars.readReal<float>("Rt/beam/IsocenterToBeamDistance");
    }
}


void Patient_Parameters_t::print()
{
    std::cout << "Files:" << std::endl;
    std::cout << "    - Patient: " << patient_dir << std::endl;
    std::cout << "    - Beams:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
        std::cout << "        " << beam_dirs.at(i) << std::endl;
    std::cout << "    - Topas files:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
        std::cout << "        " << topas_files.at(i) << std::endl;
    std::cout << "    - Tramp files:" << std::endl;
    for (size_t i = 0; i < nbeams; i++)
        std::cout << "        " << tramp_files.at(i) << std::endl;

    std::cout << "Beam angles:" << std::endl;
    for (size_t i = 0; i < nbeams; i++) {
        std::cout << "    - " << beam_names.at(i) << ":" << std::endl;
        std::cout << "        - gantry: " << angles.at(i).gantry*180.0/M_PI << " deg" << std::endl;
        std::cout << "        - couch:  " << angles.at(i).couch*180.0/M_PI  << " deg" << std::endl;
    }
    std::cout << "Isocenter to beam distance:" << std::endl;
    for (size_t i = 0; i < nbeams; i++) {
        std::cout << "    - " << beam_names.at(i) << ": " << isocenter_to_beam_distance.at(i) << " cm" << std::endl;
    }
    std::cout << "CT data:" << std::endl;
    std::cout << "    - Voxels:     " << ct.n.x;
    std::cout << ", " << ct.n.y;
    std::cout << ", " << ct.n.z << std::endl;
    std::cout << "    - Voxel size: " << ct.d.x;
    std::cout << ", " << ct.d.y;
    std::cout << ", " << ct.d.z << " cm" << std::endl;
    std::cout << "    - Offset:     " << ct.offset.x;
    std::cout << ", " << ct.offset.y;
    std::cout << ", " << ct.offset.z << " cm" << std::endl;

}

void Patient_Parameters_t::ct_to_int_coordinates()
{
    std::swap(ct.d.x, ct.d.z);
    std::swap(ct.n.x, ct.n.z);
    std::swap(ct.offset.x, ct.offset.z);
    ct.offset.x *= -1;

    // From the center to the corner
    ct.offset.x -= 0.5*ct.n.x*ct.d.x;
    ct.offset.y -= 0.5*ct.n.y*ct.d.y;
    ct.offset.z -= 0.5*ct.n.z*ct.d.z;

    for (size_t i = 0; i < angles.size(); i++) {
        angles[i].couch *= -1.0; //reverse couch angle
        angles[i].gantry = (270.0 * (M_PI/180.0)) - angles[i].gantry;
    }
}

void Patient_Parameters_t::ct_to_ext_coordinates()
{
    // From the center to the corner
    ct.offset.z += 0.5*ct.n.z*ct.d.z;
    ct.offset.y += 0.5*ct.n.y*ct.d.y;
    ct.offset.x += 0.5*ct.n.x*ct.d.x;

    ct.offset.x *= -1;
    std::swap(ct.d.x, ct.d.z);
    std::swap(ct.n.x, ct.n.z);
    std::swap(ct.offset.x, ct.offset.z);

    for (size_t i = 0; i < angles.size(); i++) {
        angles[i].couch *= -1.0; //reverse couch angle
        angles[i].gantry = (270.0 * (M_PI/180.0)) - angles[i].gantry;
    }
}

void Patient_Parameters_t::update_geometry_with_external(const Volume_t& v)
{
    consolidate_originals();

    ct.offset = ct.offset + 0.5*(ct.n*ct.d - v.n*v.d);
    ct.n.x = v.n.x;
    ct.n.y = v.n.y;
    ct.n.z = v.n.z;
    ct.d.x = v.d.x;
    ct.d.y = v.d.y;
    ct.d.z = v.d.z;
    ct.total = v.nElements;
}

void Patient_Parameters_t::set_spots_data()
{
    spots_per_field.reserve(nbeams);
    short acc = 0;
    for (size_t i=0; i < nbeams; i++) {
        Tramp_t tramp(tramp_files.at(i));

        spots_per_field.push_back(tramp.nspots);
        acc += tramp.nspots;
        accu_spots_per_field.push_back(acc);

        float temp = tramp.energies.front();
        std::vector<short> idxs;
        for (size_t j=1; j < tramp.nspots; j++) {
            if (temp != tramp.energies.at(j))
                idxs.push_back(j);
        }
        energy_layers.push_back(idxs);
    }
    // Set total spots
    total_spots = accu_spots_per_field.back();
}

void Patient_Parameters_t::set_treatment_planes()
{
    // External to internal coordinates would do: x -> -y; y -> -x, they are all 0 here.
    treatment_planes = Planes_t(nbeams);
    for (size_t i = 0; i < nbeams; i++) {
        // Set default
        treatment_planes.p.at(i) = Vector3_t<float>(0, 0, -isocenter_to_beam_distance.at(i));
        treatment_planes.dir.at(i) = Vector3_t<float>(0, 0, 1);
        treatment_planes.source_a.at(i) = Vector3_t<float>(0, 0, -std::abs(virtualSAD.a));
        treatment_planes.source_b.at(i) = Vector3_t<float>(0, 0, -std::abs(virtualSAD.b));
        // Rotate
        treatment_planes.p.at(i).rotate(angles.at(i).gantry, angles.at(i).couch);
        treatment_planes.dir.at(i).rotate(angles.at(i).gantry, angles.at(i).couch);
        treatment_planes.source_a.at(i).rotate(angles.at(i).gantry, angles.at(i).couch);
        treatment_planes.source_b.at(i).rotate(angles.at(i).gantry, angles.at(i).couch);
        // Add offsets
        treatment_planes.p.at(i) -= ct.offset;
        treatment_planes.source_a.at(i) -= ct.offset;
        treatment_planes.source_b.at(i) -= ct.offset;
    }
}
