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
    set_new_tramps();
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
        run_dir.at(i) = beam_dirs.at(i);
        run_dir.at(i).append("/run");
    }

    // Get tramp files
    tramp_files.resize(nbeams);
    std::vector<std::string> temp = getFilesWithSuffix(tramp_dir, ".tramp_modified");
    if (temp.size() == 0)
        temp = getFilesWithSuffix(tramp_dir, ".tramp");
    if (temp.size() != nbeams) {
        std::cerr << "ERROR! The number of tramp files found != number of beams!" << std::endl;
        std::cerr << "Found beams:" << std::endl;
        for (uint i = 0; i < nbeams; ++i)
            std::cerr << beam_dirs.at(i) << " -> " << temp.at(i) << std::endl;
        std::cerr << "Found tramps:" << std::endl;
        for (uint i = 0; i < temp.size(); ++i)
            std::cerr << temp.at(i) << std::endl;
        exit(EXIT_FAILURE);
    }
    std::copy(temp.begin(), temp.end(), tramp_files.begin());

    // Get topas parameter files
    topas_files.resize(nbeams);
    for (unsigned int i = 0; i < nbeams; i++) {
        topas_files.at(i) = getFilesWithSuffix(beam_dirs[i], "MCAUTO_DICOM.txt")[0];
    }
}

void Patient_Parameters_t::set_new_tramps()
{
    geometric_tramp_names.resize(nbeams);
    adapted_tramp_names.resize(nbeams);
    for (uint i = 0; i < nbeams; ++i) {
        adapted_tramp_names.at(i) = std::string("adapted_beam_" + std::to_string(i) + "_" +
                                                 beam_names.at(i) + ".tramp_modified");
        geometric_tramp_names.at(i) = std::string("geometric_beam_" + std::to_string(i) + "_" +
                                                  beam_names.at(i) + ".tramp_modified");
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
    machine = utils::toLower(pars.readString("Rt/beam/TreatmentMachineName"));
    virtualSAD.a = pars.readReal<float>("Rt/beam/VirtualSourceAxisDistances0");
    virtualSAD.b = pars.readReal<float>("Rt/beam/VirtualSourceAxisDistances1");

    // CT location
    ct.origin.x = pars.readReal<float>("Rt/CT/ImagePositionPatient0");
    ct.origin.y = pars.readReal<float>("Rt/CT/ImagePositionPatient1");
    ct.origin.z = pars.readReal<float>("Rt/CT/ImagePositionPatient2");

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
    beam_names.resize(nbeams);
    n_fractions.resize(nbeams);
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

        // Fractions
        n_fractions.at(i) = pars.readInteger<uint>("Rt/beam/NumberOfFractionsPlanned", 1);

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
    for (size_t i = 0; i < nbeams; i++) {
        std::cout << "        " << beam_dirs.at(i) << std::endl;
    }
    std::cout << "    - Topas files:" << std::endl;
    for (size_t i = 0; i < nbeams; i++) {
        std::cout << "        " << topas_files.at(i) << std::endl;
    }
    std::cout << "    - Source files (" << total_spots << " spots):" << std::endl;
    for (size_t i = 0; i < nbeams; i++) {
        std::cout << "        " << tramp_files.at(i) << ": ";
        std::cout << spots_per_field.at(i) << " spots" << std::endl;
    }

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
    std::swap(ct.isocenter.x, ct.isocenter.z);
    std::swap(ct.origin.x, ct.origin.z);
    ct.offset.x *= -1;
    ct.isocenter.x *= -1;
    ct.origin.x *= -1;

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
    ct.isocenter.x *= -1;
    ct.origin.x *= -1;
    std::swap(ct.d.x, ct.d.z);
    std::swap(ct.n.x, ct.n.z);
    std::swap(ct.offset.x, ct.offset.z);
    std::swap(ct.isocenter.x, ct.isocenter.z);
    std::swap(ct.origin.x, ct.origin.z);

    for (size_t i = 0; i < angles.size(); i++) {
        angles[i].couch *= -1.0; //reverse couch angle
        angles[i].gantry = angles[i].gantry - (270.0 * (M_PI/180.0));
    }
}

void Patient_Parameters_t::update_geometry(const Volume_t& v)
{
    consolidate_originals();

    ct.offset = ct.offset + 0.5*(ct.n*ct.d - v.n*v.d);
    ct.origin.x = v.origin.x;
    ct.origin.y = v.origin.y;
    ct.origin.z = v.origin.z;
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
        source_energies.reserve(source_energies.size() + tramp.energies.size());
        source_energies.insert(source_energies.end(), tramp.energies.begin(), tramp.energies.end());
        source_weights.reserve(source_weights.size() + tramp.weights.size());
        source_weights.insert(source_weights.end(), tramp.weights.begin(), tramp.weights.end());

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
