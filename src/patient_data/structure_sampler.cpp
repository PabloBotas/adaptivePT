#include "structure_sampler.hpp"

#include "vector4.hpp"
#include "volume.hpp"
#include "warper.hpp"

#include <algorithm>    // std::sort
#include <random>
#include <string>
#include <vector>


template<class T>
void structure_sampler (const std::string& file, const ushort nprobes,
                        const ushort nspots, Warper_t warper, const CT_Dims_t& pat_ct,
                        std::vector<T>& ct_pos, std::vector<T>& cbct_pos)
{
    std::cout << "Launching structure sampler " << file << std::endl;
    // Read structure volume
    Volume_t vol(file, Volume_t::Source_type::MHA);
    vol.ext_to_int_coordinates();

    // Store active voxel positions in array
    Array3<double> positions;
    Array3<double> boundary_positions;    
    Vector3_t<double> pos_temp;
    Vector3_t<int> idx_temp;
    for (uint index = 0; index < vol.nElements; ++index) {
        if (vol.data.at(index) > 0.5) {
            // Get index
            idx_temp.x = index/(vol.n.y*vol.n.z) % vol.n.x;
            idx_temp.y = index/vol.n.z % vol.n.y;
            idx_temp.z = index % vol.n.z;
            // Check if it's boundary of target and store it
            bool is_boundary = false;
            for (int i = -1; i < 2; ++i) {
                for (int j = -1; j < 2; ++j) {
                    for (int k = -1; k < 2; ++k) {
                        if (i == 0 && j == 0 && k == 0)
                            continue;
                        if (idx_temp.x+i >= (int)vol.n.x || idx_temp.x+i < 0 ||
                            idx_temp.y+j >= (int)vol.n.y || idx_temp.y+j < 0 ||
                            idx_temp.z+k >= (int)vol.n.z || idx_temp.z+k < 0) {
                            is_boundary = true;
                            break;
                        }
                        uint index_bound = (idx_temp.z+k) + (idx_temp.y+j)*vol.n.z +
                                           (idx_temp.x+i)*vol.n.z*vol.n.y;
                        if (vol.data.at(index_bound) < 0.5) {
                            is_boundary = true;
                            break;
                        }
                    }
                    if (is_boundary)
                        break;
                }
                if (is_boundary)
                    break;
            }
            pos_temp = vol.d * (idx_temp + 0.5);
            if (is_boundary) {
                boundary_positions.push_back(pos_temp);
            }
            positions.push_back(pos_temp);
        }
    }

    // Adjust array size
    if (ct_pos.size() < (size_t)nprobes*nspots)
        ct_pos.resize(nprobes*nspots);
    if (ct_pos.size() > (size_t)nprobes*nspots) {
        ct_pos.resize(nprobes*nspots);
        ct_pos.shrink_to_fit();
    }

    // Sample nprobes voxels from array
    std::default_random_engine generator(0);
    std::uniform_int_distribution<uint> distribution(0, positions.size()-1);
    std::vector<uint> rand_indexes(nprobes);
    for (ushort i = 0; i < nprobes; ++i) {
        rand_indexes.at(i) = distribution(generator);
    }
    std::sort (rand_indexes.begin(), rand_indexes.end());

    for (ushort i = 0; i < nspots; ++i) {
        for (ushort j=0; j < nprobes; ++j) {
            uint index = i*nspots + j;
            ct_pos.at(index).x = positions.at(rand_indexes.at(j)).x;
            ct_pos.at(index).y = positions.at(rand_indexes.at(j)).y;
            ct_pos.at(index).z = positions.at(rand_indexes.at(j)).z;
        }
    }

    // CBCT:
    // Adjust array size
    if (cbct_pos.size() < (size_t)nprobes*nspots)
        cbct_pos.resize(nprobes*nspots);
    if (cbct_pos.size() > (size_t)nprobes*nspots) {
        cbct_pos.resize(nprobes*nspots);
        cbct_pos.shrink_to_fit();
    }

    // Fill influence array for CBCT
    Array4<double> temp(ct_pos.begin(), ct_pos.begin()+nspots);
    Array4<double> temp2 = warper.apply_to_points(temp, pat_ct);
    for (ushort i = 0; i < nspots; ++i)
        std::copy (temp2.begin(), temp2.end(), cbct_pos.begin()+i*nspots);
}

template
void structure_sampler< Vector4_t<double> >(const std::string&, const ushort,
                                            const ushort, Warper_t, const CT_Dims_t&,
                                            std::vector< Vector4_t<double> >&, std::vector< Vector4_t<double> >&);