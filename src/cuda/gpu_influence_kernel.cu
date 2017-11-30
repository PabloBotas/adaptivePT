#include "gpu_influence_kernel.cuh"

#include "gpu_ray_kernel.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"


__global__ void get_influence_kernel(const uint nspots,
                                     const uint nprobes,
                                     float4* influence,
                                     float* spot_weights,
                                     float* inf_volume,
                                     float *new_energies)
{
    const uint thread = blockIdx.x*blockDim.x + threadIdx.x;
    const uint spot_i = thread / nprobes;
    const uint probe_j = thread % nprobes;
    // Influence of spot i on probe j
    if(spot_i < nspots && probe_j < nprobes) {
        float wepl_r = -1.0, wepl_d = -1.0;

        // Initialize probe ray
        float4 tempvxdata = vxdata[spot_i];
        if (new_energies)
            tempvxdata.w = new_energies[spot_i];
        Ray ray(xdata[spot_i], tempvxdata, ixdata[spot_i]);
        const float initial_energy = ray.get_initial_energy();
        // Stop condition: projection of endpoint on ray trajectory
        float3 point0, point1;
        point1 = make_float3(influence[probe_j]);
        point0 = ray.get_position() +
                 dot((point1 - ray.get_position()), ray.get_direction()) *
                 ray.get_direction();

        float dose = 0;
        if (point0.x >= 0 && point0.x < ctVoxSize.x*ctVox.x &&
            point0.y >= 0 && point0.y < ctVoxSize.y*ctVox.y &&
            point0.z >= 0 && point0.z < ctVoxSize.z*ctVox.z) {
            
            wepl_d = wepl_to_point (ray, point0);
            if (wepl_d >= 0) {
                ray.set_direction(point1 - ray.get_position());
                wepl_r = wepl_to_point (ray, point1, true);
            }
            dose = get_dose_at (initial_energy, wepl_d, wepl_r);
        }

        influence[nprobes*spot_i+probe_j].x = influence[probe_j].x;
        influence[nprobes*spot_i+probe_j].y = influence[probe_j].y;
        influence[nprobes*spot_i+probe_j].z = influence[probe_j].z;
        influence[nprobes*spot_i+probe_j].w = dose*spot_weights[spot_i];

        // Accumulate influence of spot i on probe j position
        atomicAdd(&inf_volume[get_voxel_abs (point1)], dose);
    }
}


__device__ float wepl_to_point (Ray& ray, float3 stop_point, bool lateral_wepl)
{
    float wepl = -1;
    ray.set_wepl(0);
    if (lateral_wepl)
        ray.set_energy(160000000);

    // Get voxel
    int4 vox = get_voxel (ray.get_position());
    VoxelUpdater voxUpdater;
    VoxelStepper voxStepper;

    // Raytrace to stop_point
    while (vox.w != -1) {
        float step_water = 0, step = 0, de = 0;
        float max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                      vox, voxUpdater, voxStepper, stop_point);
        if (lateral_wepl)
            get_step(step, step_water, de, max_step, ray.get_energy(), vox);
        else
#if defined(__STEP_CENTRAL_AXIS__)
            get_step(step, step_water, de, max_step, ray.get_energy(), vox);
#elif defined(__STEP_Q50__)
            get_q50_blur_step(step, step_water, de, max_step, ray, vox);
#else
            get_average_blur_step(step, step_water, de, max_step, ray, vox);
#endif
        if (lateral_wepl)
            ray.move(step, step, 0);
        else
            ray.move(step, step_water, de);

        if (voxUpdater == NONE) {
            wepl = ray.get_wepl();
            break;
        }
        else if (!ray.is_alive())
            break;
        else
            changeVoxel(vox, voxUpdater, voxStepper);
    }

    return wepl;
}

