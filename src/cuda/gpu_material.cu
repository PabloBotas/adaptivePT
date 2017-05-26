#include "gpu_material.cuh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

std::vector<float> readDensityCorrect(const std::string file)
//  read density correction factor
{
    std::cout << "densityCorrect: Reading " << file << std::endl;
    std::ifstream stream(file);
    if (!stream.is_open()) {
        std::cerr << "Can't open file: " << file << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::string dummy;
    // dummy lines 
    std::getline(stream, line);
    // Get number of factors
    std::getline(stream, line);
    size_t const nfactors = stoi(line);
    std::getline(stream, line);
    std::getline(stream, line);
    std::istringstream ss(line);

    std::vector<float> dcfactor;
    dcfactor.reserve(nfactors);
    for (size_t i = 0; i < nfactors; i++)
    {
        std::string s;
        getline(ss, s, ',');
        dcfactor.push_back(stoi(s));
    }

    return dcfactor;
}


float HU2dens(const short val)
{
//    convert HU to dens, in g/cm^3
    //    MGH calibration curve
    if (val >= -1000 && val < -98)
        return 0.00121 + 0.001029700665188*(1000.0 + val);
    else if (val >= -98 && val < 15)
        return 1.018 + 0.000893*val;
    else if (val >= 15 && val < 23)
        return 1.03;
    else if (val >= 23 && val < 101)
        return 1.003 + 0.001169*val;
    else if (val >= 101 && val < 2001)
        return 1.017 + 0.000592*val;
    else if (val >= 2001 && val < 2995)
        return 2.201 + 0.0005*(-2000.0 + val);
    else
        return 4.54;
}

int HU2matId(const int val, const std::vector<int> hu_indexes)
{
//  conert HU to material id according to hu_indexes
    for(size_t i = 0; i < hu_indexes.size()-1; i++)
    {
        if(val >= hu_indexes[i] && val < hu_indexes[i+1])
            return i;
    }
    return hu_indexes.size()-1;
}

