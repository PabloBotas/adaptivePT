#include "gpu_influence_kernel.cuh"

#include "gpu_ray_kernel.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"


__global__ void get_influence_kernel(const short num,
                                     const short* spots_per_field,
                                     const double4* const endpos,
                                     double4* influence)
{
    const uint thread = blockIdx.x*blockDim.x + threadIdx.x;
    const uint spot_i = thread / num;
    const uint probe_j = thread % num;
    if(spot_i < num && probe_j < num)
    {
        double wepl_r = -1.0, wepl_d = -1.0;

        if (true)
        {
            // Initialize probe ray
            Ray ray(xdata[spot_i], vxdata[spot_i], ixdata[spot_i]);
            // Stop condition: projection of endpoint on ray trajectory
            double3 point0, point1;
            point1 = make_double3(endpos[probe_j]);
            point0 = ray.get_position() +
                     dot((point1 - ray.get_position()), ray.get_direction()) *
                     ray.get_direction();

            // if (spot_i < 1 && probe_j < 40)
            //     printf("%u - %u %u - %f %f %f - %f %f %f - %f %f %f\n",
            //         spot_i, probe_j, num*spot_i+probe_j,
            //         xdata[spot_i].x, xdata[spot_i].y, xdata[spot_i].z,
            //         point0.x, point0.y, point0.z,
            //         point1.x, point1.y, point1.z);

            // Get WEPL to point0
            wepl_d = wepl_to_point (ray, point0);

            if (wepl_d > 0)
            {
                ray.set_direction(point1 - ray.get_position());
                wepl_r = wepl_to_point (ray, point1);
            }

            if (spot_i < 1 && probe_j < 40)
                printf("%u - %u %u - %f %f %f - %f %f %f - %f %f\n",
                    spot_i, probe_j, num*spot_i+probe_j,
                    point0.x, point0.y, point0.z,
                    point1.x, point1.y, point1.z,
                    wepl_r, wepl_d);
        }

        double influ = get_influence (wepl_r, wepl_d);

        influence[num*spot_i+probe_j].x = wepl_r;
        influence[num*spot_i+probe_j].y = wepl_r;
        influence[num*spot_i+probe_j].z = wepl_d;
        influence[num*spot_i+probe_j].w = influ;
    }
}


__device__ double wepl_to_point (Ray& ray, double3 stop_point)
{
    double wepl = -1;
    ray.set_wepl(0);
    ray.set_energy(190000000);

    // Get voxel
    int4 vox = get_voxel (ray.get_position());
    VoxelUpdater voxUpdater;
    VoxelStepper voxStepper;

    // Raytrace to stop_point
    while (vox.w != -1)
    {
        double step_water = 0, step = 0, dummy = 0;
        double max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                      vox, voxUpdater, voxStepper, stop_point);
        get_step(step, step_water, dummy, max_step, ray.get_energy(), vox);
        ray.move(step, step_water, 0);

        if (voxUpdater == NONE)
        {
            wepl = ray.get_wepl();
            break;
        }
        changeVoxel(vox, voxUpdater, voxStepper);
    }

    return wepl;
}


__device__ double get_influence (double wepl_r, double wepl_d)
{
    double influence = -1.0;
    // Get nominal influence at center axis of integrated peak from LUT

    // Get radial decay
    if (influence > 0.0)
    {
        // Get initial sigma
        // double sigma_0 = ;
        // Get halo sigma
        double sigma_H = 6.5-3.4*wepl_d+0.078*wepl_d*wepl_d; // wepl in cm, sigma in mm
        // Get MCS sigma
        // double sigma_MCS = ;
        // Total sigmas
        double sigma1_2 = sigma_MCS*sigma_MCS + sigma_0*sigma_0;
        double sigma2_2 = sigma_H*sigma_H + sigma_0*sigma_0;
        // Alpha
        double t = d/R;
        double a0 =  1.002E+0 - 5.900E-4*R;
        double a1 =  2.128E-3 - 2.044E-2*R + 3.178E-4*R*R;
        double a2 = -2.459E-3 + 2.125E-2*R - 3.788E-4*R*R;
        double alpha = min(1, a0 + a1*t + a2*t*t);
    }

    return influence;
}
