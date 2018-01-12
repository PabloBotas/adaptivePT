#include "gpu_ray_kernel.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"


__global__ void raytrace_plan_kernel(const ushort num,
                                     const short* spots_per_field,
                                     const float4* const orig_endpoints,
                                     float4 *pos_scorer,
                                     float* traces)
{
    const int thread = blockIdx.x*blockDim.x + threadIdx.x;
    if (thread < num) {
        Ray ray(xdata[thread], vxdata[thread], ixdata[thread]);
        size_t const ind = get_endpoints_index(ray.get_beam_id(),
                                               ray.get_spot_id(),
                                               spots_per_field);
        ray.set_wepl(0);
        int4 vox = get_voxel(ray.get_position());
        bool if_correct_energy = orig_endpoints ? true : false;
        bool entered_mask = false;
        bool vox_in_mask = false;
        float3 sample_pos = ray.get_position();;
        float3 initial_pos = ray.get_position();;
        int sample_vox = vox.w;

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;
        const float initial_energy = ray.get_energy();
        pos_scorer[ind].w = initial_energy;

        while (ray.is_alive() && vox.w >= -1) {
            float step_water = 0, step = 0, de = 0;
            float max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                         vox, voxUpdater, voxStepper);
#if defined(__STEP_CENTRAL_AXIS__)
            get_step(step, step_water, de, max_step, ray.get_energy(), vox);
#elif defined(__STEP_Q50__)
            get_q50_blur_step(step, step_water, de, max_step, ray, vox);
#else
            get_average_blur_step(step, step_water, de, max_step, ray, vox);
#endif
            ray.move(step, step_water, de);
            if (masking_vf) {
                if (vox_in_mask || (!vox_in_mask && !entered_mask)) {
                    sample_pos = ray.get_position();
                    sample_vox = vox.w;
                }
                if (vox_in_mask && !entered_mask) {
                    entered_mask = true;
                }
            }
            // This line returns to the old behavior without the mask
            // sample_pos = ray.get_position();

            if (traces && (thread == 0 || thread == 300)) {
                score_traces(thread, traces, vox.w, !ray.is_alive());
            }
            if (step == max_step) {
                changeVoxel(vox, voxUpdater, voxStepper);
                if (masking_vf && vox.w != -1) {
                    vox_in_mask = tex3D(vf_mask_tex, vox.z, vox.y, vox.x) > 0.5;
                }
            }
        }

        // Save scorer
        if (traces && (thread == 0 || thread == 300)) {
            score_traces(thread, traces, sample_vox, true);
        }
        if (!masking_vf) {
            sample_pos = ray.get_position();
        }
        pos_scorer[ind].x = sample_pos.x;
        pos_scorer[ind].y = sample_pos.y;
        pos_scorer[ind].z = sample_pos.z;

        if (if_correct_energy) {
            const float accu_wepl = ray.get_wepl();
            float3 plan_endpoint = make_float3(orig_endpoints[thread]);
            ray.set_energy(initial_energy);
            ray.set_wepl(0);
            int sign = ahead_or_behind(ray.dir, plan_endpoint, ray.pos);
            ray.dir *= sign;

            while (ray.is_alive() && vox.w != -1) {
                float step_water = 0, step = 0, de = 0;
                float max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                             vox, voxUpdater, voxStepper, plan_endpoint);
#if defined(__STEP_CENTRAL_AXIS__)
                get_step(step, step_water, de, max_step, ray.get_energy(), vox);
#elif defined(__STEP_Q50__)
                get_q50_blur_step(step, step_water, de, max_step, ray, vox);
#else
                get_average_blur_step(step, step_water, de, max_step, ray, vox);
#endif
                ray.move(step, step_water, de);

                if (voxUpdater != NONE) {
                    if (traces) {
                        int tempvox = sign < 0 ? -vox.w : vox.w;
                        score_traces(thread, traces, tempvox, false);
                    }
                    changeVoxel(vox, voxUpdater, voxStepper);
                } else {
                    if (traces)
                        score_traces(thread, traces, vox.w, true);
                    break;
                }
            }

            const float total_wepl = accu_wepl + sign*ray.get_wepl();
            const float i_alpha1 = 6.11788728306583E+04;
            const float i_alpha2 = 1.29373145117383E+03;
            const float i_alpha3 = 3.25651117101543E+13;
            const float i_p1 = 5.55293543382075E-01;
            const float i_p2 = 1.49459384779856E+00;
            const float i_p3 = 5.55245491049867E-01;
            const float E = powf(total_wepl*i_alpha1, i_p1) + powf(total_wepl*i_alpha2, i_p2) +
                            powf(total_wepl*i_alpha3, i_p3);

            // float dev = 100*(E - initial_energy)/initial_energy;
            // // Coerce change to 0.25% steps
            // const float energy_grid_step = 0.0f;
            // float delta = E - initial_energy;
            // if (energy_grid_step > 0.0f) {
            //     dev = floor(dev/energy_grid_step + 0.5)*energy_grid_step;
            //     delta = dev*initial_energy/100;
            // }
            pos_scorer[ind].w = E;
        }
    }
}

__device__ void score_traces(int thread, float *traces, int voxnum, bool last)
{
    bool del_content = voxnum < 0;
    if (del_content) {
        atomicExch(&traces[abs(voxnum)], 0.0f);
    } else {
        float val = last ? 100000.f : float(thread);
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
