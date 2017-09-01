// TODO read doses
// TODO read tramps
// TODO transfer data to GPU (thurst)
// TODO Get Dijs

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <cublas_v2.h>
#include <curand.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "tramp.hpp"

#define IDX2C(i,j,ld) (((j)*(ld))+(i))

void read_dose (std::vector<float>& dose, std::string& file, uint nvox=0);
void process_commands (int argc, char* argv[],
                       std::vector<std::string>& tramp_files,
                       std::vector<std::string>& dose_files,
                       std::string& out_dir,
                       std::vector<ushort>& dims);
void cublas_operations(cublasHandle_t& handle,
                       float *dij, const float *dose, const float *weights,
                       const int nvox, const int nspots);

int main(int argc, char* argv[])
{
    std::vector<std::string> tramp_files;
    std::vector<std::string> dose_files;
    std::string out_dir;
    std::vector<ushort> dims;
    process_commands(argc, argv, tramp_files, dose_files, out_dir, dims);
    
    size_t nfiles = tramp_files.size();
    uint nx = dims.at(0);
    uint ny = dims.at(1);
    uint nz = dims.at(2);
    uint nvox = nx*ny*nz;

    // Create a handle for CUBLAS
    cublasHandle_t handle;
    cublasCreate(&handle);

    std::vector<float> dose(nvox);
    for (size_t i = 0; i < nfiles; i++)
    {
        Tramp_t tramp(tramp_files.at(i));
        std::vector<float> weights = tramp.get_weights();
        const uint nspots = weights.size();
        read_dose(dose, dose_files.at(i), nvox);
        std::vector<float> dij(nvox*nspots, 0);

        cublas_operations(handle,
            &dij[0], dose.data(), weights.data(),
            nvox, nspots);
    }

    cublasDestroy(handle);

    return EXIT_SUCCESS;
}

// Multiply the arrays A and B on GPU and save the result in C
// W(m,n) = w(m,k)*(w(m,k)*w(m,k))'
// D(m,n) = A(m,k) * W(k,n)
void cublas_operations(cublasHandle_t& handle,
                       float *dij, const float *dose, const float *weights,
                       const int nvox, const int nspots)
{
    const float alpha = 1.0;
    const float beta = 0.0;

    cublasStatus_t status;

    // Set main device matrices
    float* d_weights;
    float* d_dose;
    float* d_dij;
    cudaMalloc ((void**)&d_weights, nspots*sizeof(float));
    cudaMalloc ((void**)&d_dose, nvox*sizeof(float));
    cudaMalloc ((void**)&d_dij, nvox*nspots*sizeof(float));
    cublasSetMatrix (nspots, 1, sizeof(float), weights, nspots, d_weights, nspots);
    cublasSetMatrix (nvox, 1, sizeof(float), dose, nvox, d_dose, nvox);
    cublasSetMatrix (nvox, nspots, sizeof(float), dij, nvox, d_dij, nvox);

    // Set helper device matrices
    std::vector<float> h_wwt(nspots*nspots, 0);
    std::vector<float> h_wt_wwtinv(nspots, 0);
    float* wwt;
    float* wwtinv;
    float* wt_wwtinv;
    cudaMalloc ((void**)&wwt, nspots*nspots*sizeof(float));
    cudaMalloc ((void**)&wwtinv, nspots*nspots*sizeof(float));
    cudaMalloc ((void**)&wt_wwtinv, nspots*sizeof(float));
    cublasSetMatrix (nspots, nspots, sizeof(float), wwt, nspots,
        h_wwt.data(), nspots);
    cublasSetMatrix (nspots, nspots, sizeof(float), wwtinv, nspots,
        h_wwt.data(), nspots);
    cublasSetMatrix (1, nspots, sizeof(float), wt_wwtinv, 1,
        h_wt_wwtinv.data(), 1);

    // Perform calculations --------------------
    // Multiply w and wt
    status = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T,
                         nspots, nspots, 1,
                         &alpha,
                         d_weights, nspots,
                         d_weights, 1,
                         &beta,
                         wwt, nspots);
    // Invert wwt
    int *d_info;
    int info;
    float* wwtinv_array[] = { wwtinv };
    status = cublasSmatinvBatched(handle,
                                  nspots, &wwt, nspots,
                                  wwtinv_array, nspots,
                                  d_info, 1);
    cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    // Multiply wt and inverse(wwt)
    status = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T,
                         nspots, nspots, 1,
                         &alpha,
                         d_weights, nspots,
                         d_weights, 1,
                         &beta,
                         wwt, nspots);

    cudaFree(d_weights); cudaFree(d_dose); cudaFree(d_dij);
    cudaFree(wwt); cudaFree(wwtinv); cudaFree(wt_wwtinv);
}

void read_dose (std::vector<float>& dose, std::string& file, uint nvox)
{
    std::ifstream ifs(file, std::ios::in | std::ios::binary);
    if (!ifs.is_open()) {
        std::cerr << "Cannot open file: " << file << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Reading input file ..." << std::endl;

    if (nvox == 0)
    {
        ifs.seekg (0, ifs.end);
        nvox = ifs.tellg() / sizeof(float);
        ifs.seekg (0, ifs.beg);
        dose.resize(nvox);
    }

    ifs.read((char*)&dose[0], nvox*sizeof(float));
    ifs.close();
}

void process_commands (int argc, char* argv[],
                       std::vector<std::string>& tramp_files,
                       std::vector<std::string>& dose_files,
                       std::string& out_dir,
                       std::vector<ushort>& dims)
{
    try 
    {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "Produce this help message.")
        ("tramps", po::value< std::vector<std::string> >(&tramp_files)->multitoken()->required(), "Tramp files.")
        ("doses",  po::value< std::vector<std::string> >(&dose_files)->multitoken()->required(), "Dose files.")
        ("outdir", po::value<std::string>(&out_dir)->required(), "Output directory to write Dijs to.")
        ("dims",   po::value<std::vector<ushort>>(&dims)->multitoken()->required(), "Dose cube dimensions.");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cerr << desc << std::endl;
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);

        if ( tramp_files.size() != dose_files.size() ) {
            std::cerr << "ERROR! Number of tramp and dose files are not consistent." << std::endl;
            std::cerr << desc << std::endl;
            exit(EXIT_SUCCESS);
        }
    }
    catch(std::exception& e) {
        std::cerr << "ERROR! " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}