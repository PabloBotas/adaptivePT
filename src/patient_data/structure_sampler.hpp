#ifndef __STRUCTURE_SAMPLER_HPP__
#define __STRUCTURE_SAMPLER_HPP__

#include "warper.hpp"

#include <string>
#include <vector>

template<class T>
void structure_sampler (const std::string& file, const ushort nprobes,
                        const ushort nspots, Warper_t warper, const CT_Dims_t& pat_ct,
                        std::vector<T>& ct_pos, std::vector<T>& cbct_pos);

#endif