#include "gpu_ray_kernel.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_ray_class.cuh"
#include "gpu_geometry_tools.cuh"

__global__ void raytrace_plan_kernel(const int num,
                                     const short* spots_per_field,
                                     float4 *pos_scorer,
                                     float4 *ini_scorer,
                                     short2 *meta_scorer,
                                     float* traces)
{
    const int id = blockIdx.x*blockDim.x + threadIdx.x;
    if(id < num)
    {
        // float const max_energy_loss = 0.2; // % of pre-step energy
        Ray ray(xdata[id], vxdata[id], ixdata[id]);
        size_t ind = get_endpoints_index(ray.beam_id, ray.spot_id, spots_per_field);
        int4 vox;
        vox.x = floor(ray.pos.x/ctVoxSize.x);
        vox.y = floor(ray.pos.y/ctVoxSize.y);
        vox.z = floor(ray.pos.z/ctVoxSize.z);
        vox.w = getVoxelIndex(vox);

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;

        if(ini_scorer)
        {
            ini_scorer[ind].x = ray.pos.x;
            ini_scorer[ind].y = ray.pos.y;
            ini_scorer[ind].z = ray.pos.z;
            ini_scorer[ind].w = ray.energy;
        }

        if(meta_scorer)
        {
            meta_scorer[ind].x = ray.beam_id;
            meta_scorer[ind].y = ray.spot_id;
        }

        while (ray.isAlive() && vox.w != -1)
        {
            if(traces)
                atomicAdd(&traces[vox.w], 1.f);
            float max_step = inters(ray, vox, voxUpdater, voxStepper);
            float step_water, step;
            getWaterStep(step, step_water, max_step, ray.energy, ray.wepl, vox);
            ray.move(step, step_water);
            changeVoxel(vox, voxUpdater, voxStepper);
        }

        pos_scorer[ind].x = ray.pos.x;
        pos_scorer[ind].y = ray.pos.y;
        pos_scorer[ind].z = ray.pos.z;
        pos_scorer[ind].w = ray.initial_wepl;
    }
}

__global__ void backtrace_endpoints_kernel(const int num,
                                           const short* spots_per_field,
                                           float4 *pos_scorer,
                                           float* traces)
{
    const int id = blockIdx.x*blockDim.x + threadIdx.x;
    if(id < num)
    {
        Ray ray(xdata[id], vxdata[id], ixdata[id]);
        ray.reverse();
        ray.setMaxWEPL(ray.initial_wepl);
        ray.wepl = 0;
        int4 vox;
        vox.x = floor(ray.pos.x/ctVoxSize.x);
        vox.y = floor(ray.pos.y/ctVoxSize.y);
        vox.z = floor(ray.pos.z/ctVoxSize.z);
        vox.w = getVoxelIndex(vox);

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;

        while (ray.isAlive() && vox.w != -1)
        {
            if(traces)
                atomicAdd(&traces[vox.w], 1.f);
            float max_step = inters(ray, vox, voxUpdater, voxStepper);
            float step_water, step;
            getWaterStep(step, step_water, max_step, ray.energy, ray.wepl, vox);
            ray.move_back(step, step_water);
            changeVoxel(vox, voxUpdater, voxStepper);
        }

        size_t ind = get_endpoints_index(ray.beam_id, ray.spot_id, spots_per_field);
        pos_scorer[ind].x = ray.pos.x;
        pos_scorer[ind].y = ray.pos.y;
        pos_scorer[ind].z = ray.pos.z;
        pos_scorer[ind].w = ray.wepl;
    }
}

__device__ void getWaterStep(float& step,
                             float& step_water,
                             const float max_step,
                             const float energy,
                             const float avail_wepl,
                             const int4& vox)
{
    // Get density
    float const density  = tex3D(dens_tex, vox.z, vox.y, vox.x);
    // Get stp ratio
    float const mass_stp_ratio = massStpRatio(energy, vox);

    // Set steps
    step = max_step;
    step_water = mass_stp_ratio*density*max_step;
    // Verify it's not too much
    if (step_water > avail_wepl)
    {
        step_water = avail_wepl;
        step = step_water/mass_stp_ratio/density;
    }
}

__device__ float massStpRatio(const float energy, const int4& vox)
//  mass stopping power ratio wrt water for a material
{
    float const energy_index = (energy - stp_ratio_min_e)/stp_ratio_delta_e + 0.5;
    int const material_id = tex3D(matid_tex, vox.z, vox.y, vox.x);
    int const material_index = material_id + 0.5;
    return tex2D(stp_ratio_tex, energy_index , material_index);
}


__device__ unsigned int get_endpoints_index(const short beam_id,
                                            const short spot_id,
                                            const short* spots_per_field)
{
    unsigned int index = spot_id;
    for (int ibeam = 0; ibeam < beam_id; ibeam++)
    {
        index += spots_per_field[ibeam];
    }
    return index;
}
