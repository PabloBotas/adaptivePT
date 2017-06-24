#include "gpu_ray_kernel.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"

__global__ void raytrace_plan_kernel(const short num,
                                     const short* spots_per_field,
                                     const float4* const orig_endpoints,
                                     float4 *pos_scorer,
                                     float* traces)
{
    const int thread = blockIdx.x*blockDim.x + threadIdx.x;
    if(thread < num)
    {
        Ray ray(xdata[thread], vxdata[thread], ixdata[thread]);
        size_t ind = get_endpoints_index(ray.get_beam_id(),
                                         ray.get_spot_id(),
                                         spots_per_field);
        int4 vox = get_voxel (ray.get_position());

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;

        pos_scorer[ind].w = ray.get_energy();

        // ray.print();
        while (ray.is_alive() && vox.w != -1)
        {
            if(traces)
                atomicAdd(&traces[vox.w], 1.0f);

            float step_water = 0, step = 0, de = 0;
            float max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                         vox, voxUpdater, voxStepper);
            get_water_step(step, step_water, de, max_step,
                       ray.get_energy(), vox);

            ray.move(step, step_water, de);

            if (step == max_step)
                changeVoxel(vox, voxUpdater, voxStepper);
        }
        // Save scorer
        pos_scorer[ind].x = ray.pos.x;
        pos_scorer[ind].y = ray.pos.y;
        pos_scorer[ind].z = ray.pos.z;

        if(orig_endpoints)
        {
            // printf("Converging to endpoint!!!\n");
            const float sample_energy = 160*MeV2eV;
            const float sample_wepl   = 17.82000; // Janni table for sample energy
                    
            float3 plan_endpoint = make_float3(orig_endpoints[thread]);
            ray.set_energy(sample_energy); // sample energy
            ray.set_wepl(sample_wepl);     // Janni table for sample energy
            ray.dir *= ahead_or_behind(ray.dir, plan_endpoint, ray.pos);

            while (ray.is_alive() && vox.w != -1)
            {
                float step_water, step;
                float max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                             vox, voxUpdater, voxStepper,
                                             plan_endpoint);
                float de = 0;
                get_water_step(step, step_water, de, max_step,
                               ray.get_energy(), vox);
                ray.move(step, step_water, de);
                if (voxUpdater != NONE)
                    changeVoxel(vox, voxUpdater, voxStepper);
                else
                    break;
            }

            pos_scorer[ind].w = sample_energy - ray.get_energy();
        }
    }

    
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
