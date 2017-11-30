#include "structure_sampler.hpp"

#include "vector4.hpp"
#include "volume.hpp"
#include "warper.hpp"

#include <algorithm>    // std::sort, std::min
#include <fstream>
#include <limits>
#include <numeric>      // std::iota
#include <random>
#include <string>
#include <vector>


template<class T>
void structure_sampler (const std::string& file, const float& percentage,
                        const uint nspots, Warper_t warper, const CT_Dims_t& pat_ct,
                        std::vector<T>& ct_pos, std::vector<T>& cbct_pos)
{
    // Read structure volume
    Volume_t vol(file, Volume_t::Source_type::MHA);
    vol.ext_to_int_coordinates();

    // Store active voxel positions in array
    Array3<float> positions;
    // Array3<float> boundary_positions;
    for (uint index = 0; index < vol.nElements; ++index) {
        if (vol.data.at(index) > 0.5) {
            Vector3_t<float> pos_temp;
            Vector3_t<uint> idx_temp;
            // Get index
            idx_temp.x = index/(vol.n.y*vol.n.z) % vol.n.x;
            idx_temp.y = index/vol.n.z % vol.n.y;
            idx_temp.z = index % vol.n.z;
            // // Check if it's boundary of target and store it
            // bool is_boundary = false;
            // for (int i = -1; i < 2; ++i) {
            //     for (int j = -1; j < 2; ++j) {
            //         for (int k = -1; k < 2; ++k) {
            //             if (i == 0 && j == 0 && k == 0)
            //                 continue;
            //             if (idx_temp.x+i >= (int)vol.n.x || idx_temp.x+i < 0 ||
            //                 idx_temp.y+j >= (int)vol.n.y || idx_temp.y+j < 0 ||
            //                 idx_temp.z+k >= (int)vol.n.z || idx_temp.z+k < 0) {
            //                 is_boundary = true;
            //                 break;
            //             }
            //             uint index_bound = (idx_temp.z+k) + (idx_temp.y+j)*vol.n.z +
            //                                (idx_temp.x+i)*vol.n.z*vol.n.y;
            //             if (vol.data.at(index_bound) < 0.5) {
            //                 is_boundary = true;
            //                 break;
            //             }
            //         }
            //         if (is_boundary)
            //             break;
            //     }
            //     if (is_boundary)
            //         break;
            // }
            if (idx_temp.x >= vol.n.x || idx_temp.y >= vol.n.y || idx_temp.z >= vol.n.z) {
                std::cerr << "ERROR in structure sampling, out of the volume!!" << std::endl;
                exit(EXIT_FAILURE);
            }
            pos_temp.x = vol.d.x * (idx_temp.x + 0.5);
            pos_temp.y = vol.d.y * (idx_temp.y + 0.5);
            pos_temp.z = vol.d.z * (idx_temp.z + 0.5);
            positions.push_back(pos_temp);
            // if (is_boundary) {
            //     boundary_positions.push_back(pos_temp);
            // }
        }
    }

    // Calculate number of probes
    const uint nprobes = uint(std::ceil(positions.size()*std::min(1.0f, percentage/100)));

    // Adjust arrays size
    if (ct_pos.size() < (size_t)nprobes*nspots && cbct_pos.size() < (size_t)nprobes*nspots) {
        // std::cout << "Allocating " << nprobes*nspots << " x2 elements" << std::endl;
        ct_pos.resize(nprobes*nspots);
        cbct_pos.resize(nprobes*nspots);
    } else if (ct_pos.size() > (size_t)nprobes*nspots && cbct_pos.size() > (size_t)nprobes*nspots) {
        // std::cout << "Allocating " << nprobes*nspots << " x2 elements" << std::endl;
        ct_pos.resize(nprobes*nspots);
        ct_pos.shrink_to_fit();
        cbct_pos.resize(nprobes*nspots);
        cbct_pos.shrink_to_fit();
    }

    // Sample nprobes voxels from array
    std::cout << "Sampling " << nprobes;
    std::vector<uint> rand_indexes(positions.size());
    std::iota (rand_indexes.begin(), rand_indexes.end(), 0);
    if (nprobes != positions.size()) {
        std::cout << " random";
        // std::srand(std::time(0));
        std::srand(0);
        std::random_shuffle (rand_indexes.begin(), rand_indexes.end());
        rand_indexes.resize(nprobes);
        std::sort(rand_indexes.begin(), rand_indexes.end());
    }
    std::cout << " points (";
    std::cout << 100*float(nprobes)/positions.size() << " %) of structure in " << file << std::endl;

    for (uint i=0; i < nspots; ++i) {
        for (uint j=0; j < nprobes; ++j) {
            uint index = i*nspots + j;
            ct_pos.at(index).x = positions.at(rand_indexes.at(j)).x;
            ct_pos.at(index).y = positions.at(rand_indexes.at(j)).y;
            ct_pos.at(index).z = positions.at(rand_indexes.at(j)).z;
        }
    }

    // CBCT:
    // Fill influence array for CBCT
    Array4<float> temp(ct_pos.begin(), ct_pos.begin()+nprobes);
    Array4<float> temp2 = warper.apply_to_points(temp, pat_ct);
    for (uint i = 0; i < nspots; ++i)
        std::copy (temp2.begin(), temp2.end(), cbct_pos.begin()+i*nprobes);
}

template
void structure_sampler< Vector4_t<float> >(const std::string&, const float&,
                                            const uint, Warper_t, const CT_Dims_t&,
                                            std::vector< Vector4_t<float> >&,
                                            std::vector< Vector4_t<float> >&);