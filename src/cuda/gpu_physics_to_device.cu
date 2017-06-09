#include "gpu_physics_to_device.cuh"

#include "gpu_device_globals.cuh"
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
    std::string file = "../src/phys_data/mass_stopping_power_ratio.dat";
    std::cout << "sendMassStoppingPowerRatio: Reading " << file << std::endl;
    std::ifstream stream(file);
    if (!stream.is_open()) {
        std::cerr << "Can't open file: " << file << std::endl;
        exit(EXIT_FAILURE);
    }

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

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaMallocArray(&stp_ratio_array, &channelDesc, n_energies, n_materials);
    cudaMemcpyToArray(stp_ratio_array, 0,0, &stp_ratios[0], sizeof(float)*stp_ratios.size(), cudaMemcpyHostToDevice);
    stp_ratio_tex.filterMode = cudaFilterModeLinear;
    cudaBindTextureToArray(stp_ratio_tex, stp_ratio_array, channelDesc);
}

void freePhysicsMemory()
{
    cudaFreeArray(stp_ratio_array);
    cudaUnbindTexture(stp_ratio_tex);
}
