#include "gpu_ray_kernel.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_ray_class.cuh"
#include "gpu_geometry_tools.cuh"

__global__ void raytrace_plan_kernel(const int num,
                                     const short* spots_per_field,
                                     const float4* orig_endpoints,
                                     float4 *pos_scorer,
                                     float* traces)
{
    const int id = blockIdx.x*blockDim.x + threadIdx.x;
    if(id < num)
    {
        // float const max_energy_loss = 0.2; // % of pre-step energy
        Ray ray(xdata[id], vxdata[id], ixdata[id]);
        size_t ind = get_endpoints_index(ray.get_beam_id(), ray.get_spot_id(), spots_per_field);
        int4 vox = get_voxel (ray.get_position());

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;

        pos_scorer[ind].w = ray.get_energy();

        // ray.print();
        
        while (ray.is_alive() && vox.w != -1)
        {
            if(traces)
                atomicAdd(&traces[vox.w], 1.0f);

            float step_water, step;
            float max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                         vox, voxUpdater, voxStepper);
            getWaterStep(step, step_water, max_step,
                         ray.get_energy(), ray.get_wepl(), vox);
            ray.move(step, step_water);
            if (step == max_step)
                changeVoxel(vox, voxUpdater, voxStepper);
        }

        // Save scorer
        pos_scorer[ind].x = ray.pos.x;
        pos_scorer[ind].y = ray.pos.y;
        pos_scorer[ind].z = ray.pos.z;

        if (orig_endpoints)
        {
            // printf("Converging to endpoint!!!\n");
            const float sample_energy = 160*MeV2eV;
            const float sample_wepl   = 17.82000; // Janni table for sample energy
            float3 plan_endpoint = make_float3(orig_endpoints[id]);
            ray.set_energy(sample_energy); // sample energy
            ray.set_wepl(sample_wepl);     // Janni table for sample energy
            ray.dir *= ahead_or_behind(ray.dir, plan_endpoint, ray.pos);

            while (ray.is_alive() && vox.w != -1)
            {
                float step_water, step;
                float max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                             vox, voxUpdater, voxStepper,
                                             plan_endpoint);
                getWaterStep(step, step_water, max_step,
                             ray.get_energy(), ray.get_wepl(), vox);
                ray.move(step, step_water);
                if (voxUpdater != NONE)
                    changeVoxel(vox, voxUpdater, voxStepper);
                else
                    break;
            }

            pos_scorer[ind].w = sample_energy - ray.get_energy();
        }

        // ray.print();
    }
}

__device__ void getWaterStep(float& step,
                             float& step_water,
                             const float max_step,
                             const float energy_in,
                             const float avail_wepl,
                             const int4& vox)
{
    // Get density
    float const density  = tex3D(dens_tex, vox.z, vox.y, vox.x);
    // Get stp ratio
    float const mass_stp_ratio = massStpRatio(energy_in, vox);

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


__device__ size_t get_endpoints_index(const short beam_id,
                                      const short spot_id,
                                      const short* spots_per_field)
{   
    unsigned int index = spot_id;
    for (short ibeam = 0; ibeam < beam_id; ibeam++)
    {
        index += spots_per_field[ibeam];
    }
    return index;
}
