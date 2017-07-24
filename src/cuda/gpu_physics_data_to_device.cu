#include "gpu_physics_data_to_device.cuh"

#include "gpu_device_globals.cuh"
#include "utils.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

void gpu_physics_to_device::sendMassStoppingPowerRatio()
{
    std::vector<int> HU_starting_values;
    gpu_physics_to_device::sendMassStoppingPowerRatio(HU_starting_values);
}

void gpu_physics_to_device::sendMassStoppingPowerRatio(std::vector<int>& HU_starting_values)
{
    //  read mass stopping power ratio
    std::string file = std::string(INSTALLATION_PATH) + "/src/phys_data/mass_stopping_power_ratio.dat";
    std::cout << "sendMassStoppingPowerRatio: Reading " << file << std::endl;
    std::ifstream stream(file);
    utils::check_fs(stream, file, "to read mass stopping powers.");

    std::string line;
    std::string dummy;
    // Two dummy lines 
    std::getline(stream, line);
    std::getline(stream, line);
    // Get number of materials
    std::getline(stream, line);
    size_t const n_materials = stoi(line);

    HU_starting_values.resize(n_materials);
    std::vector<float> stp_ratios;
    size_t n_energies;
    float minimum_energy;
    float delta_energy;
    
    // Read data
    for (size_t imat = 0; imat < n_materials; imat++)
    {
        // Get number of energies per material
        std::getline(stream, line);
        std::getline(stream, line);
        std::istringstream ss(line);
        ss >> dummy >> HU_starting_values.at(imat) >>
              dummy >> n_energies >> minimum_energy >> delta_energy;
        std::getline(stream, line);

        if (imat == 0)
        {
            stp_ratios.resize(n_energies*n_materials);
        }

        for (size_t i = 0; i < n_energies; i++)
        {
            std::getline(stream, line);
            std::istringstream ss(line);
            ss >> dummy >> stp_ratios.at(i + imat*n_energies);
        }
    }

    //  transfer to GPU
    minimum_energy *= MeV2eV;
    delta_energy *= MeV2eV;
    cudaMemcpyToSymbol(stp_ratio_min_e, &minimum_energy, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(stp_ratio_delta_e, &delta_energy, sizeof(float), 0, cudaMemcpyHostToDevice);

    cudaMallocArray(&stp_ratio_array, &stp_ratio_tex.channelDesc, n_energies, n_materials);
    cudaMemcpyToArray(stp_ratio_array, 0,0, &stp_ratios[0], sizeof(float)*stp_ratios.size(), cudaMemcpyHostToDevice);
    stp_ratio_tex.filterMode = cudaFilterModeLinear;
    cudaBindTextureToArray(stp_ratio_tex, stp_ratio_array);
}

void gpu_physics_to_device::sendWaterRestrictedSPower()
{
    std::string file = std::string(INSTALLATION_PATH) + "/src/phys_data/nist_stopping_power_water.dat";
    std::cout << "sendWaterRestrictedSPower: Reading " << file << std::endl;
    std::ifstream stream(file);
    utils::check_fs(stream, file, "to read restricted stopping powers.");

    std::string line;
    // Two dummy lines 
    std::getline(stream, line);
    std::getline(stream, line);
    // Get next
    std::getline(stream, line);
    float dummy, min_energy, delta_energy;
    size_t ndata;
    std::istringstream ss(line);
    ss >> dummy >> min_energy >> dummy >> delta_energy >> ndata;

    min_energy   *= MeV2eV;
    delta_energy *= MeV2eV;

    //  read
    std::vector<float> stp_w(ndata);
    std::vector<float> stp_w_b(ndata);
    //  read info
    std::getline(stream, line);
    std::getline(stream, line);

    for(size_t i = 0; i < ndata; i++)
    {
        std::getline(stream, line, '\n');
        ss.clear();
        ss.str(line);
        ss >> dummy >> dummy >> dummy >> stp_w[i] >> stp_w_b[i];
        stp_w[i] *= MeV2eV;
    }

    cudaMemcpyToSymbol(stp_w_min_e, &min_energy, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(stp_w_delta_e, &delta_energy, sizeof(float), 0, cudaMemcpyHostToDevice);

    //  pass to GPU
    cudaMallocArray(&stp_w_array, &stp_w_tex.channelDesc, ndata, 1);
    cudaMemcpyToArray(stp_w_array, 0, 0, &stp_w[0], sizeof(float)*ndata, cudaMemcpyHostToDevice);
    stp_w_tex.filterMode = cudaFilterModeLinear;
    cudaBindTextureToArray(stp_w_tex, stp_w_array);

    cudaMallocArray(&stp_w_b_coeff_array, &stp_w_b_coeff_tex.channelDesc, ndata, 1);
    cudaMemcpyToArray(stp_w_b_coeff_array, 0, 0, stp_w_b.data(), sizeof(float)*ndata, cudaMemcpyHostToDevice);
    stp_w_b_coeff_tex.filterMode = cudaFilterModeLinear;
    cudaBindTextureToArray(stp_w_b_coeff_tex, stp_w_b_coeff_array);
}

void freePhysicsMemory()
{
    cudaFreeArray(stp_w_array);
    cudaUnbindTexture(stp_w_tex);
    cudaFreeArray(stp_w_b_coeff_array);
    cudaUnbindTexture(stp_w_b_coeff_tex);
    cudaFreeArray(stp_ratio_array);
    cudaUnbindTexture(stp_ratio_tex);
}
