#include "gpu_influence_kernel.cuh"

#include "gpu_ray_kernel.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"


__global__ void get_influence_kernel(const short num,
                                     const double4* const endpos,
                                     double4* influence,
                                     float* spot_weights,
                                     float* inf_volume)
{
    const uint thread = blockIdx.x*blockDim.x + threadIdx.x;
    const uint spot_i = thread / num;
    const uint probe_j = thread % num;
    // Influence of spot i on probe j
    if(spot_i < num && probe_j < num)
    {
        double wepl_r = -1.0, wepl_d = -1.0;

        // Initialize probe ray
        Ray ray(xdata[spot_i], vxdata[spot_i], ixdata[spot_i]);
        // Stop condition: projection of endpoint on ray trajectory
        double3 point0, point1;
        point1 = make_double3(endpos[probe_j]);
        point0 = ray.get_position() +
                 dot((point1 - ray.get_position()), ray.get_direction()) *
                 ray.get_direction();

        // Get WEPL to point0
        wepl_d = wepl_to_point (ray, point0);

        if (wepl_d >= 0)
        {
            ray.set_direction(point1 - ray.get_position());
            wepl_r = wepl_to_point (ray, point1, true);
        }

        double influ = get_influence (vxdata[spot_i].w, wepl_r, wepl_d);

        influence[num*spot_i+probe_j].x = vxdata[spot_i].w;
        influence[num*spot_i+probe_j].y = wepl_r;
        influence[num*spot_i+probe_j].z = wepl_d;
        influence[num*spot_i+probe_j].w = influ*spot_weights[spot_i];

        // Accumulate influence of spot i on probe j position
        int4 vox = get_voxel (point1);
        atomicAdd(&inf_volume[vox.w], influ*spot_weights[spot_i]);

        // if(spot_i < 10 && probe_j < 10)
            // printf("%u %u %f %f\n", spot_i, probe_j, vxdata[spot_i].w, influence[num*spot_i+probe_j].w);
    }
}


__device__ double wepl_to_point (Ray& ray, double3 stop_point, bool overwrite_energy)
{
    double wepl = -1;
    ray.set_wepl(0);
    if (overwrite_energy)
        ray.set_energy(190000000);

    // Get voxel
    int4 vox = get_voxel (ray.get_position());
    VoxelUpdater voxUpdater;
    VoxelStepper voxStepper;

    // Raytrace to stop_point
    while (vox.w != -1)
    {
        double step_water = 0, step = 0, de = 0;
        double max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                      vox, voxUpdater, voxStepper, stop_point);
        get_step(step, step_water, de, max_step, ray.get_energy(), vox);
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


__device__ double get_influence (double const energy, double const wepl_r, double const wepl_d)
{
    if (wepl_d < 0 && wepl_r < 0)
        return 0;
    // Get nominal influence at center axis of integrated peak from LUT
    double const energy_index = (energy - bp_energy_min)/bp_energy_delta + 0.5;
    double const depth_index = wepl_d/bp_depth_delta + 0.5;
    double b = tex2D(bp_b_tex, energy_index, depth_index);
    double n = tex2D(bp_n_tex, energy_index, depth_index);
    double s = tex2D(bp_s_tex, energy_index, depth_index);
    double w = tex2D(bp_w_tex, energy_index, depth_index);

    // Get radial decay
    double influence = n * (w/(sqrtf(2*PI)*s) * exp(-0.5 * wepl_r*wepl_r/(s*s)) +
        (1-w) * 2*pow(b, 1.5)/PI *1/pow(wepl_r*wepl_r + b, 2));

    return fabs(influence);
}

