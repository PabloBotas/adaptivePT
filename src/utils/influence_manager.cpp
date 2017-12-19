#include "influence_manager.hpp"

#include "gpu_influence_launchers.cuh"

#include <algorithm>


Influence_manager::Influence_manager() : ct_metadata(Volume_metadata_t())
{
}


Influence_manager::~Influence_manager()
{
}


Influence_manager::Influence_manager(const Parser& parser_, const Patient_Parameters_t& pat,
                                     Warper_t warper_, const Volume_metadata_t ct_metadata_)
                                    : patient_parameters(&pat), ctdims(&pat.ct),
                                      ct_metadata(ct_metadata_)
{
    engine = parser_.influence_opts;
    outputdir = parser_.output_opt4D_files.empty() ? parser_.out_dir :
                parser_.output_opt4D_files;
    nspots = patient_parameters->total_spots;

    // Needed for structure sampling
    ct_mask_file = parser_.ct_mask_file;
    warper = warper_;
}


void Influence_manager::calculate_at_plan()
{
    calculate(matrix_at_plan, pos_at_plan);
}


void Influence_manager::calculate_at_fraction(std::vector<float> new_energies)
{
    calculate(matrix_at_fraction, pos_at_fraction, new_energies);
}


void Influence_manager::calculate(Array4<float>& influence, Array3<float>& positions,
                                  std::vector<float> new_energies)
{
    if (engine == Influence_engines_t::BEAM_MODEL) {
        // Fill matrix plan with starting positions
        if (matrix_elements == 0)
            structure_sampler();

        influence.resize(matrix_elements);

        Array4<float> temp(n_probing_positions);
        for (uint i = 0; i < n_probing_positions; ++i)
            temp[i] = positions[i];

        for (uint i = 0; i < nspots; ++i) {
            std::copy(temp.begin(), temp.end(),
                      influence.begin()+i*n_probing_positions);
        }

        if (save_memory) {
            positions.clear();
            boundary_at_plan.clear();
        }
        
        influence_from_beam_model (influence, new_energies);
    } else {
        // NOT IMPLEMENTED YET!!
    }
}


void Influence_manager::influence_from_beam_model (Array4<float>& influence,
                                                   const std::vector<float>& new_energies)
{
    influence_from_beam_model_launcher(influence, new_energies, ct_metadata,
        *patient_parameters, nspots, n_probing_positions, outputdir);
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
    n_probing_positions = uint(std::ceil(positions.size()*std::min(1.0f, mask_sampling_percentage/100)));
    matrix_elements = nspots*n_probing_positions;

    // Adjust arrays size
    if (pos_at_plan.size() < (size_t)n_probing_positions &&
        pos_at_fraction.size() < (size_t)n_probing_positions) {
        // std::cout << "Allocating " << n_probing_positions << " elements" << std::endl;
        pos_at_plan.resize(n_probing_positions);
        pos_at_fraction.resize(n_probing_positions);
    } else if (pos_at_plan.size() > (size_t)n_probing_positions &&
        pos_at_fraction.size() > (size_t)n_probing_positions) {
        // std::cout << "Allocating " << n_probing_positions << " elements" << std::endl;
        pos_at_plan.resize(n_probing_positions);
        pos_at_plan.shrink_to_fit();
        pos_at_fraction.resize(n_probing_positions);
        pos_at_fraction.shrink_to_fit();
    }

    // Sample n_probing_positions voxels from array
    std::cout << "Sampling " << n_probing_positions;
    std::vector<uint> rand_indexes(positions.size());
    std::iota (rand_indexes.begin(), rand_indexes.end(), 0);
    if (n_probing_positions != positions.size()) {
        std::cout << " random";
        // std::srand(std::time(0));
        std::srand(0);
        std::random_shuffle (rand_indexes.begin(), rand_indexes.end());
        rand_indexes.resize(n_probing_positions);
        std::sort(rand_indexes.begin(), rand_indexes.end());
    }
    std::cout << " points (" << 100*float(n_probing_positions)/positions.size();
    std::cout << " %) of structure in " << ct_mask_file << std::endl;

    for (uint i=0; i < n_probing_positions; ++i) {
        pos_at_plan.at(i).x = positions.at(rand_indexes.at(i)).x;
        pos_at_plan.at(i).y = positions.at(rand_indexes.at(i)).y;
        pos_at_plan.at(i).z = positions.at(rand_indexes.at(i)).z;
    }

    // CBCT:
    // Fill positions at fraction
    Array4<float> temp(pos_at_plan.begin(), pos_at_plan.end());
    temp = warper.apply_to_points(temp, *ctdims);
    std::copy(temp.begin(), temp.end(), pos_at_fraction.begin());    
}
