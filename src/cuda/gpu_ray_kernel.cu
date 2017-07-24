#include "gpu_ray_kernel.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"


__global__ void raytrace_plan_kernel(const short num,
                                     const short* spots_per_field,
                                     const double4* const orig_endpoints,
                                     double4 *pos_scorer,
                                     float* traces)
{
    const int thread = blockIdx.x*blockDim.x + threadIdx.x;
    if(thread < num)
    {
        Ray ray(xdata[thread], vxdata[thread], ixdata[thread]);
        size_t const ind = get_endpoints_index(ray.get_beam_id(),
                                               ray.get_spot_id(),
                                               spots_per_field);
        ray.set_wepl(0);
        int4 vox = get_voxel (ray.get_position());

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;
        const double initial_energy = ray.get_energy();
        pos_scorer[ind].w = initial_energy;

        // ray.print();
        while (ray.is_alive() && vox.w >= -1)
        {
            double step_water = 0, step = 0, de = 0;
            double max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                          vox, voxUpdater, voxStepper);
            get_step(step, step_water, de, max_step, ray.get_energy(), vox);
            ray.move(step, step_water, de);

            if(traces)
                score_traces(traces, vox.w, !ray.is_alive());

            if (step == max_step)
                changeVoxel(vox, voxUpdater, voxStepper);
        }

        // Save scorer
        pos_scorer[ind].x = ray.pos.x;
        pos_scorer[ind].y = ray.pos.y;
        pos_scorer[ind].z = ray.pos.z;

        if(orig_endpoints)
        {
            const double accu_wepl = ray.get_wepl();
            const double sample_energy = initial_energy;
            const double sample_wepl   = 0;
                    
            double3 plan_endpoint = make_double3(orig_endpoints[thread]);
            ray.set_energy(sample_energy);
            ray.set_wepl(sample_wepl);
            int sign = ahead_or_behind(ray.dir, plan_endpoint, ray.pos);
            ray.dir *= sign;

            while (ray.is_alive() && vox.w != -1)
            {
                double step_water = 0, step = 0, de = 0;
                double max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                              vox, voxUpdater, voxStepper, plan_endpoint);
                get_step(step, step_water, de, max_step, ray.get_energy(), vox);
                ray.move(step, step_water, de);

                if (voxUpdater != NONE)
                {
                    if (traces)
                    {
                        int tempvox = sign < 0 ? -vox.w : vox.w;
                        score_traces(traces, tempvox, false);
                    }
                    changeVoxel(vox, voxUpdater, voxStepper);
                }
                else
                {
                    if (traces)
                        score_traces(traces, vox.w, true);
                    break;
                }
            }

            const double total_wepl = accu_wepl + sign*ray.get_wepl();
            const double alpha1 = 1.63455120e-05;
            const double alpha2 = 7.72957942e-04;
            const double alpha3 = 3.07077098e-14;
            const double p1 = 1.80084932;
            const double p2 = 0.669078092;
            const double p3 = 1.80100517;
            const double E = pow(total_wepl/alpha1, 1/p1) + pow(total_wepl/alpha2, 1/p2) + pow(total_wepl/alpha3, 1/p3);
            const double delta_E = E - initial_energy;
            pos_scorer[ind].w = delta_E;
        }
    }
}


__device__ void score_traces(float *traces, int voxnum, bool last)
{
    bool del_content = voxnum < 0;
    if (del_content)
    {
        atomicExch(&traces[abs(voxnum)], 0.0f);
    }
    else
    {
        float val = last ? 50.0f : 1.0f;
        atomicExch(&traces[voxnum], val);
    }
}


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
