#include "influence_manager.hpp"

#include "gpmc_manager.hpp"
#include "gpu_influence_launchers.cuh"
#include "utils.hpp"

#include <algorithm>
#include <fstream>


Influence_manager::Influence_manager() : ct_metadata(Volume_metadata_t())
{
}


Influence_manager::~Influence_manager()
{
}


Influence_manager::Influence_manager(const Parser& parser_, const Patient_Parameters_t& pat,
                                     Warper_t warper_, const Volume_metadata_t ct_metadata_) :
    engine(parser_.influence_opts),
    n_spots(patient_parameters->total_spots),
    dose_plan_file(parser_.dose_plan_file),
    dose_frac_file(parser_.dose_frac_file),
    // gpmc_dij_plan_file(parser_.gpmc_dij_plan_file),
    gpmc_dij_frac_file(parser_.gpmc_dij_frac_file),
    beam_model_dij_plan_file(parser_.beam_model_dij_plan_file),
    beam_model_dij_frac_file(parser_.beam_model_dij_frac_file),
    patient_parameters(&pat),
    ctdims(&pat.ct),
    ct_metadata(ct_metadata_),
    outputdir(parser_.out_dir),
    ct_mask_file(parser_.ct_mask_file),
    // Needed for structure sampling
    warper(warper_)
{
    // engine = parser_.influence_opts;
    // outputdir = parser_.out_dir;
    // dose_plan_file = parser_.dose_plan_file;
    // dose_frac_file = parser_.dose_frac_file;
    // // gpmc_dij_plan_file = parser_.gpmc_dij_plan_file;
    // gpmc_dij_frac_file = parser_.gpmc_dij_frac_file;
    // beam_model_dij_plan_file = parser_.beam_model_dij_plan_file;
    // beam_model_dij_frac_file = parser_.beam_model_dij_frac_file;
    // gpmc_dij_file = parser_.out_dir
    // n_spots = patient_parameters->total_spots;

    // // Needed for structure sampling
    // ct_mask_file = parser_.ct_mask_file;
    // warper = warper_;
}


void Influence_manager::get_dose_at_plan()
{
    get_dose(dose_plan_file, matrix_at_plan, dose_at_plan);
}


void Influence_manager::get_dose_at_frac(const std::vector<float>& new_energies)
{
    get_dose(dose_frac_file, matrix_at_frac, dose_at_frac, new_energies);
}


void Influence_manager::get_dose(std::string dose_file,
                                 Array4<float>& matrix, std::vector<float>& dose,
                                 const std::vector<float>& new_energies)
{
    if (engine == Influence_engines_t::BEAM_MODEL) {
        if (!matrix.size())
            get_dij_at_frac(new_energies);
        // Accumulate influence on position j by all spots i(0 -> n)
        dose.resize(n_voxels);
        for (size_t j = 0; j<n_voxels; j++)
            for(size_t i = 0; i<n_spots; i++)
                dose.at(j) += matrix.at(n_voxels*i+j).w;
    } else if (engine == Influence_engines_t::GPMC_DIJ) {
        read_dose_file(dose_file, dose);
    } else if (engine == Influence_engines_t::GPMC_DOSE) {
        read_dose_file(dose_file, dose);
    }
}


void Influence_manager::get_dij_at_plan()
{
    if (engine == Influence_engines_t::BEAM_MODEL) {
        get_dij(beam_model_dij_plan_file, matrix_at_plan, pos_at_plan);
    } else {
        std::cerr << "Dij at plan is only needed in \"beam model\" mode." << std::endl;
        std::cerr << "Why are you even requesting this anyway?" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Influence_manager::get_dij_at_frac(const std::vector<float>& new_energies)
{
    get_dij(beam_model_dij_frac_file, matrix_at_frac, pos_at_frac, new_energies);
}


void Influence_manager::get_dij(std::string outfile,
                                Array4<float>& influence, Array3<float>& positions,
                                const std::vector<float>& new_energies)
{
    if (engine == Influence_engines_t::BEAM_MODEL) {
        // Fill matrix plan with starting positions
        if (n_spots*n_voxels == 0)
            structure_sampler();

        influence.resize(n_spots*n_voxels);

        Array4<float> temp(n_voxels);
        for (uint i = 0; i < n_voxels; ++i)
            temp[i] = positions[i];

        for (uint i = 0; i < n_spots; ++i) {
            std::copy(temp.begin(), temp.end(),
                      influence.begin()+i*n_voxels);
        }

        if (save_memory) {
            positions.clear();
            boundary_at_plan.clear();
        }
        
        influence_from_beam_model (outfile, influence, new_energies);
    } else if (engine == Influence_engines_t::GPMC_DIJ) {
        Gpmc_manager gpmc(*patient_parameters, outfile, "dosedij", outputdir);
        gpmc.calculate_dij(1000, false, ct_mask_file);
    } else if (engine == Influence_engines_t::GPMC_DOSE) {
        // NOT IMPLEMENTED YET!!
    }
}


void Influence_manager::influence_from_beam_model (std::string outfile, Array4<float>& influence,
                                                   const std::vector<float>& new_energies)
{
    influence_from_beam_model_launcher(outfile, influence, new_energies, ct_metadata,
        *patient_parameters, n_spots, n_voxels, outputdir);
}


void Influence_manager::read_dose_file (std::string file, std::vector<float>& dose)
{
    std::ifstream stream(file, std::ios::binary | std::ios::ate);
    utils::check_fs(stream, file, "to read plan dose.");
    n_voxels = stream.tellg()/sizeof(float);
    size_t bytes_to_read = n_voxels*sizeof(float);
    stream.seekg(0, std::ios::beg);
    dose.resize(n_voxels);
    stream.read(reinterpret_cast<char*>(&dose[0]), bytes_to_read);
}

void Influence_manager::structure_sampler ()
{
    // Read structure volume
    Volume_t vol(ct_mask_file, Volume_t::Source_type::MHA);
    vol.ext_to_int_coordinates();

    // Store active voxel positions in array
    Array3<float> positions;
    const float mask_threshold = 0.5;
    for (uint index = 0; index < vol.nElements; ++index) {
        if (vol.data.at(index) > mask_threshold) {
            Vector3_t<float> pos_temp;
            Vector3_t<int> idx_temp;
            // Get index
            idx_temp.x = index/(vol.n.y*vol.n.z) % vol.n.x;
            idx_temp.y = index/vol.n.z % vol.n.y;
            idx_temp.z = index % vol.n.z;
            // Check if it's boundary of target and store it
            bool is_boundary = false;
            if (get_boundary) {
                for (int i = -1; i <= 1; ++i) {
                    for (int j = -1; j <= 1; ++j) {
                        for (int k = -1; k <= 1; ++k) {
                            if (i != 0 && j != 0 && k != 0) {
                                if (idx_temp.x+i >= (int)vol.n.x || idx_temp.x+i < 0 ||
                                    idx_temp.y+j >= (int)vol.n.y || idx_temp.y+j < 0 ||
                                    idx_temp.z+k >= (int)vol.n.z || idx_temp.z+k < 0) {
                                    is_boundary = true;
                                    break;
                                }
                                uint index_bound = (idx_temp.z+k) + (idx_temp.y+j)*vol.n.z +
                                                   (idx_temp.x+i)*vol.n.z*vol.n.y;
                                if (vol.data.at(index_bound) < mask_threshold) {
                                    is_boundary = true;
                                    break;
                                }
                            }
                            
                        }
                        if (is_boundary)
                            break;
                    }
                    if (is_boundary)
                        break;
                }
            }
            if (idx_temp.x >= (int)vol.n.x ||
                idx_temp.y >= (int)vol.n.y ||
                idx_temp.z >= (int)vol.n.z) {
                std::cerr << "ERROR in structure sampling, out of the volume!!" << std::endl;
                exit(EXIT_FAILURE);
            }
            pos_temp.x = vol.d.x * (idx_temp.x + 0.5);
            pos_temp.y = vol.d.y * (idx_temp.y + 0.5);
            pos_temp.z = vol.d.z * (idx_temp.z + 0.5);
            positions.push_back(pos_temp);
            if (get_boundary && is_boundary)
                boundary_at_plan.push_back(pos_temp);
        }
    }

    // Calculate number of probes
    n_voxels = uint(std::ceil(positions.size()*std::min(1.0f, mask_sampling_percentage/100)));
    matrix_elements = n_spots*n_voxels;

    // Adjust arrays size
    if (pos_at_plan.size() < (size_t)n_voxels &&
        pos_at_frac.size() < (size_t)n_voxels) {
        // std::cout << "Allocating " << n_voxels << " elements" << std::endl;
        pos_at_plan.resize(n_voxels);
        pos_at_frac.resize(n_voxels);
    } else if (pos_at_plan.size() > (size_t)n_voxels &&
        pos_at_frac.size() > (size_t)n_voxels) {
        // std::cout << "Allocating " << n_voxels << " elements" << std::endl;
        pos_at_plan.resize(n_voxels);
        pos_at_plan.shrink_to_fit();
        pos_at_frac.resize(n_voxels);
        pos_at_frac.shrink_to_fit();
    }

    // Sample n_voxels voxels from array
    std::cout << "Sampling " << n_voxels;
    std::vector<uint> rand_indexes(positions.size());
    std::iota (rand_indexes.begin(), rand_indexes.end(), 0);
    if (n_voxels != positions.size()) {
        std::cout << " random";
        // std::srand(std::time(0));
        std::srand(0);
        std::random_shuffle (rand_indexes.begin(), rand_indexes.end());
        rand_indexes.resize(n_voxels);
        std::sort(rand_indexes.begin(), rand_indexes.end());
    }
    std::cout << " points (" << 100*float(n_voxels)/positions.size();
    std::cout << " %) of structure in " << ct_mask_file << std::endl;

    for (uint i=0; i < n_voxels; ++i) {
        pos_at_plan.at(i).x = positions.at(rand_indexes.at(i)).x;
        pos_at_plan.at(i).y = positions.at(rand_indexes.at(i)).y;
        pos_at_plan.at(i).z = positions.at(rand_indexes.at(i)).z;
    }

    // CBCT:
    // Fill positions at fraction
    Array4<float> temp(pos_at_plan.begin(), pos_at_plan.end());
    temp = warper.apply_to_points(temp, *ctdims);
    std::copy(temp.begin(), temp.end(), pos_at_frac.begin());    
}
