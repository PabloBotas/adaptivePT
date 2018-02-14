#include "gpu_ray_kernel.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_physics.cuh"
#include "gpu_ray_class.cuh"


__global__ void raytrace_plan_kernel(const ushort num,
                                     const short* spots_per_field,
                                     const float4* const orig_endpoints,
                                     float4* pos_scorer,
                                     float* traces,
                                     float3* traces_data)
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
        float3 sample_pos = ray.get_position();
        float3 initial_pos = ray.get_position();
        int sample_vox = vox.w;
        float sample_wepl = 0;

        VoxelUpdater voxUpdater;
        VoxelStepper voxStepper;
        const float initial_energy = ray.get_energy();

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
                    sample_wepl = ray.get_wepl();
                }
                if (vox_in_mask && !entered_mask) {
                    entered_mask = true;
                }
            }

            // This lines return to the old behavior without the mask
            // sample_pos = ray.get_position();
            // sample_vox = vox.w;
            // sample_wepl = ray.get_wepl();

            if (traces && (!masking_vf || length(sample_pos - ray.get_position()) < 0.001)) {
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
        if (!masking_vf) {
            sample_pos = ray.get_position();
            sample_vox = vox.w;
            sample_wepl = ray.get_wepl();
        }
        if (traces) {
            score_traces(thread, traces, sample_vox, true);
        }
        pos_scorer[ind].x = sample_pos.x;
        pos_scorer[ind].y = sample_pos.y;
        pos_scorer[ind].z = sample_pos.z;
        pos_scorer[ind].w = sample_wepl;

        if (if_correct_energy) {
            traces_data[thread*2] = make_float3(xdata[thread].x, xdata[thread].y, xdata[thread].z);
            traces_data[thread*2+1] = ray.get_position();

            const float accu_wepl = ray.get_wepl();
            float3 stop_pos = make_float3(orig_endpoints[thread]);
            float orig_wepl = orig_endpoints[thread].w;
            ray.set_energy(initial_energy);
            ray.set_wepl(0);
            int sign = ahead_or_behind(ray.dir, stop_pos, ray.pos);
            ray.dir *= sign;

            while (ray.is_alive() && vox.w != -1) {
                float step_water = 0, step = 0, de = 0;
                float max_step = to_boundary(ray.get_position(), ray.get_direction(),
                                             vox, voxUpdater, voxStepper, stop_pos);
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

            float wepl_stop = accu_wepl + sign*ray.get_wepl();
            float frac = orig_wepl/accu_wepl;
            float new_wepl = wepl_stop/frac;

            float e_r[8] = {6.10690676e-09, -1.02370849e-06, 7.18732504e-05, -2.76315988e-03,
                            6.43961648e-02, -9.85218453e-01, 1.49352594e+01, 2.27490185e+01};
            // The E0/initial factor will calibrate the curve to the vicinity of the energy
            // My observations say that I have errors of about ~0.3 % in energy,
            // but this simple approach cuts them down to 0.0001703...
            float E0 = e_r[7] + e_r[6]*accu_wepl + e_r[5]*std::pow(accu_wepl, 2) +
                       e_r[4]*std::pow(accu_wepl, 3) + e_r[3]*std::pow(accu_wepl, 4) +
                       e_r[2]*std::pow(accu_wepl, 5) + e_r[1]*std::pow(accu_wepl, 6) +
                       e_r[0]*std::pow(accu_wepl, 7);
            float E = e_r[7] + e_r[6]*new_wepl + e_r[5]*std::pow(new_wepl, 2) +
                      e_r[4]*std::pow(new_wepl, 3) + e_r[3]*std::pow(new_wepl, 4) +
                      e_r[2]*std::pow(new_wepl, 5) + e_r[1]*std::pow(new_wepl, 6) +
                      e_r[0]*std::pow(new_wepl, 7);
            E *= 1e6;
            E *= initial_energy/E0/1e6;

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
        float val = last ? 100000.f : float(thread+1);
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
