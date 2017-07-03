#include "gpu_ray_kernel.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"

#include <curand_kernel.h>

template<class T>
__global__ void raytrace_plan_kernel(const short num,
                                     const short* spots_per_field,
                                     const float4* const orig_endpoints,
                                     float4 *pos_scorer,
                                     T* traces)
{
    const int thread = blockIdx.x*blockDim.x + threadIdx.x;
    if(thread < num)
    {
        curandState localState;
        curand_init(0, thread, 0, &localState);

        Ray ray(xdata[thread], vxdata[thread], ixdata[thread]);
        size_t const ind = get_endpoints_index(ray.get_beam_id(),
                                         ray.get_spot_id(),
                                         spots_per_field);
        int4 vox = get_voxel (ray.get_position());

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;
        const float initial_energy = ray.get_energy();
        pos_scorer[ind].w = initial_energy;

        // ray.print();
        while (ray.is_alive() && vox.w >= -1)
        {
            if(traces)
                score_traces<T>(vox.w, localState, ray.get_beam_id(), ray.get_spot_id(), traces);

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
            const float sample_energy = initial_energy;
            const float sample_wepl   = 17.82000; // Janni table for 160MeV
                    
            float3 plan_endpoint = make_float3(orig_endpoints[thread]);
            ray.set_energy(sample_energy); // sample energy
            ray.set_wepl(sample_wepl);     // Janni table for sample energy
            int sign = ahead_or_behind(ray.dir, plan_endpoint, ray.pos);
            ray.dir *= sign;

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

            pos_scorer[ind].w = sign*(sample_energy - ray.get_energy());
        }
    }
}
template
__global__ void raytrace_plan_kernel(const short num,
                                     const short* spots_per_field,
                                     const float4* const orig_endpoints,
                                     float4 *pos_scorer,
                                     float* traces);
template
__global__ void raytrace_plan_kernel(const short num,
                                     const short* spots_per_field,
                                     const float4* const orig_endpoints,
                                     float4 *pos_scorer,
                                     unsigned long long int* traces);

template<class T>
__device__ void score_traces(const int& voxnum, curandState& localState,
                             const short& beamid, const short& spotid, T *traces)
{
    int rand_index = curand_uniform(&localState)*ctTotalVoxN;
    if (sizeof(traces[rand_index]) == sizeof(unsigned long long int))
    {
        unsigned long long int val = ((unsigned long long int)beamid) |
                                     ((unsigned long long int)spotid) << 4 |
                                     ((unsigned long long int)voxnum) << (4+12);
        atomicExch(&traces[voxnum], val);
    }
    else
    {
        atomicAdd(&traces[voxnum], 1.0f);
    }
}

template __device__ void score_traces(const int& voxnum, curandState& localState, const short& beamid, const short& spotid, float *traces);
template __device__ void score_traces(const int& voxnum, curandState& localState, const short& beamid, const short& spotid, unsigned long long int *traces);

__device__ size_t get_endpoints_index(const short beam_id,
                                      const short spot_id,
                                      const short* spots_per_field)
{   
    size_t accu_spots = 0;
    for (short i = 0; i < beam_id; i++)
        accu_spots += spots_per_field[i];
    size_t index = accu_spots + spot_id;
    return index;
}
