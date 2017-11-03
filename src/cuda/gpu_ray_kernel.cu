#include "gpu_ray_kernel.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"


__global__ void raytrace_plan_kernel(const ushort num,
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
            get_average_blur_step(step, step_water, de, max_step, ray, initial_energy, vox);
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
            double3 plan_endpoint = make_double3(orig_endpoints[thread]);
            ray.set_energy(initial_energy);
            ray.set_wepl(0);
            int sign = ahead_or_behind(ray.dir, plan_endpoint, ray.pos);
            ray.dir *= sign;

            while (ray.is_alive() && vox.w != -1)
            {
                double step_water = 0, step = 0, de = 0;
                double max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                              vox, voxUpdater, voxStepper, plan_endpoint);
#if defined(__STEP_Q50__)
                get_q50_blur_step(step, step_water, de, max_step, ray, initial_energy, vox);
#else
                get_average_blur_step(step, step_water, de, max_step, ray, initial_energy, vox);
#endif
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
            const float alpha1 = 1.63455120e-05;
            const float alpha2 = 7.72957942e-04;
            const float alpha3 = 3.07077098e-14;
            const float p1 = 1.80084932;
            const float p2 = 0.669078092;
            const float p3 = 1.80100517;
            const double E = pow(total_wepl/alpha1, 1/p1) + pow(total_wepl/alpha2, 1/p2) + pow(total_wepl/alpha3, 1/p3);
            double dev = 100*(E - initial_energy)/initial_energy;
            // Coerce change to 0.25% steps
            const float energy_grid_step = 0.25;
            dev = floor(dev/energy_grid_step + 0.5)*energy_grid_step;
            double delta = dev*initial_energy/100;
            pos_scorer[ind].w = delta;
        }
    }
}

// run --patient Opt4D/P01_mcgrm2425358/cbct_1 --cbct Opt4D/P01_mcgrm2425358/cbct_1/cbct_1.mha --ct_target ../patients/P01_mcgrm2425358/contours/base_plan/CTV\ OP.mha --vf ../patients/P01_mcgrm2425358/transforms/xform_deform_cCBCT1-pCT.mha --outdir Opt4D/P01_mcgrm2425358/cbct_1/adapt_free --out_vf --out_shifts --free --opt4D_out Opt4D/P01_mcgrm2425358/cbct_1/adapt_free/opt_reoptimization_files --report Opt4D/P01_mcgrm2425358/cbct_1/adapt_free/adapt_free_P01_cbct_1.pdf
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
