#include "gpu_source_positioning.cuh"

#include "initialize_rays.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_utils.cuh"

void virtual_src_to_treatment_plane(const unsigned int& num,
                                    const std::vector<BeamAngles_t>& angles,
                                    const float3& ct_offsets)
{
    std::vector<float2> temp(angles.size());
    for (size_t i = 0; i < angles.size(); i++) {
        temp[i].x = angles.at(i).gantry;
        temp[i].y = angles.at(i).couch;
    }

    float2* angles_gpu;
    array_to_device<float2>(angles_gpu, temp.data(), angles.size());

    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    virtual_src_to_treatment_plane_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, angles_gpu, ct_offsets);
    check_kernel_execution(__FILE__, __LINE__);

    cudaFree(angles_gpu);
}

__global__ void virtual_src_to_treatment_plane_kernel(const int num,
                                                      const float2* angles,
                                                      const float3 ct_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num) {
        float4 pos  = xdata[tid];
        float4 vel  = vxdata[tid];
        short2 meta = ixdata[tid]; // x = beam_id, y = spot_id
        short beamid = meta.x;

        // Adjust to internal coordinates
        pos = pos_to_int_coordinates(pos);
        vel = pos_to_int_coordinates(vel);

        //  rotate location and direction using gantry and couch angles
        float gantry = angles[beamid].x;
        float couch  = angles[beamid].y;
        pos = rotate(pos, gantry, couch);
        vel = rotate(vel, gantry, couch);

        // Add offsets
        pos.x -= ct_offsets.x;
        pos.y -= ct_offsets.y;
        pos.z -= ct_offsets.z;

        // Initialize them inside the CT
        pos = ray_trace_to_CT_volume(pos, vel);
        xdata[tid]  = pos;
        vxdata[tid] = vel;
    }
}

void correct_offsets(const unsigned int& num,
                     const float3& offsets,
                     const float3& original_offsets)
{
    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    correct_offsets_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, offsets, original_offsets);
    check_kernel_execution(__FILE__, __LINE__);
}

__global__ void correct_offsets_kernel(const int num,
                                       const float3 offsets,
                                       const float3 original_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num) {
        float4 pos  = xdata[tid];
        float4 vel  = vxdata[tid];

        // Add offsets
        pos.x += original_offsets.x - offsets.x;
        pos.y += original_offsets.y - offsets.y;
        pos.z += original_offsets.z - offsets.z;

        // Initialize them inside the CT
        pos = ray_trace_to_CT_volume(pos, vel);
        xdata[tid]  = pos;
        vxdata[tid] = vel;
    }
}

Array4<float> offset_endpoints(const Array4<float>& orig_endpoints,
                                const float3& offsets,
                                const float3& original_offsets)
{
    Array4<float> off_endpoints = orig_endpoints;
    float4* dev_orig_endpoints = NULL;
    array_to_device<float4, Vector4_t<float>>(dev_orig_endpoints, off_endpoints.data(), orig_endpoints.size());

    size_t num = off_endpoints.size();
    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    correct_offsets_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num,
                             dev_orig_endpoints,
                             offsets,
                             original_offsets);
    check_kernel_execution(__FILE__, __LINE__);

    retrieve_scorer<float, float4>(&off_endpoints[0].x, dev_orig_endpoints, num);

    cudaFree(dev_orig_endpoints);

    return off_endpoints;
}

__global__ void correct_offsets_kernel(const int num,
                                       float4* dev_orig_endpoints,
                                       const float3 offsets,
                                       const float3 original_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num) {
        dev_orig_endpoints[tid].x += original_offsets.x - offsets.x;
        dev_orig_endpoints[tid].y += original_offsets.y - offsets.y;
        dev_orig_endpoints[tid].z += original_offsets.z - offsets.z;

        // printf("%d - 0 - %f %f %f\n", 
               // tid, dev_orig_endpoints[tid].x, dev_orig_endpoints[tid].y, dev_orig_endpoints[tid].z);
    }
}

__host__ void treatment_plane_to_virtual_src(Array4<float>& pos,
                                             Array4<float> dir,
                                             const Patient_Parameters_t& pat,
                                             const Vector3_t<float>& isocenter_shift)
{
    std::vector<float2> temp(pat.angles.size());
    for (size_t i = 0; i < pat.nbeams; i++) {
        temp[i].x = pat.angles.at(i).gantry;
        temp[i].y = pat.angles.at(i).couch;
    }
    float3 offset = make_float3(pat.original_ct.offset.x, pat.original_ct.offset.y, pat.original_ct.offset.z);
    offset -= make_float3(isocenter_shift.x, isocenter_shift.y, isocenter_shift.z);

    for (size_t i = 0; i < dir.size(); i++) {
        // turn dir, containing a point, to a direction
        dir.at(i) -= pos.at(i);
        dir.at(i).normalize();

        // beam info
        short2 meta = get_beam_spot_id(i, pat.spots_per_field);
        pos.at(i).w = meta.x; // beamid
        dir.at(i).w = meta.y; // spotid
    }

    float2 *angles_gpu;
    float4 *pos_gpu, *dir_gpu;
    float3 *plane_pos_gpu, *plane_dir_gpu;
    array_to_device<float2>(angles_gpu, temp.data(), pat.nbeams);
    array_to_device<float4>(pos_gpu, pos.data(), pos.size());
    array_to_device<float4>(dir_gpu, dir.data(), dir.size());
    array_to_device<float3>(plane_pos_gpu, pat.treatment_planes.p.data(), pat.nbeams);
    array_to_device<float3>(plane_dir_gpu, pat.treatment_planes.dir.data(), pat.nbeams);

    int num = pos.size();
    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    treatment_plane_to_virtual_src_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(
                                            num,
                                            pat.nbeams,
                                            pos_gpu, dir_gpu,
                                            angles_gpu,
                                            plane_dir_gpu,
                                            plane_pos_gpu,
                                            offset);
    check_kernel_execution(__FILE__, __LINE__);

    retrieve_scorer<float, float4>(&pos[0].x, pos_gpu, pos.size());

    cudaFree(angles_gpu);
    cudaFree(pos_gpu);
    cudaFree(dir_gpu);
    cudaFree(plane_pos_gpu);
    cudaFree(plane_dir_gpu);
}

__global__ void treatment_plane_to_virtual_src_kernel(const int num,
                                                      const int nbeams,
                                                      float4* pos_,
                                                      const float4* dir_,
                                                      const float2* angles,
                                                      const float3* plane_dir,
                                                      const float3* plane_pos,
                                                      const float3 ct_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num) {
        float3 pos  = make_float3(pos_[tid]);
        float3 vel  = make_float3(dir_[tid]);
        short beamid = pos_[tid].w;
        float3 p_pos = plane_pos[beamid];
        float3 p_dir = plane_dir[beamid];

        // Ray trace to virtual source plane if necessary
        float d = dot((p_pos - pos), p_dir) / dot(vel, p_dir);
        if ( fabs(d) > 0.000001f ) {
            pos.x += d*vel.x;
            pos.y += d*vel.y;
            pos.z += d*vel.z;
        }

        // Add offsets
        pos.x += ct_offsets.x;
        pos.y += ct_offsets.y;
        pos.z += ct_offsets.z;

        //  rotate location using gantry and couch
        float gantry = angles[beamid].x;
        float couch  = angles[beamid].y;
        pos = rotate(pos, -gantry, -couch);
        // vel = rotate(vel, -gantry, -couch);

        // Adjust to external coordinates
        pos = pos_to_ext_coordinates(pos);
        // vel = int_to_ext_coordinates(vel);

        // Return to array
        pos_[tid] = make_float4(pos, pos_[tid].w);
        // dir_[tid] = vel;
    }
}

__device__ float4 ray_trace_to_CT_volume(const float4& p,
                                          const float4& v)
{
    return ray_trace_to_CT_volume(p, v, ctVox, ctVoxSize);
}

__device__ __host__ float4 ray_trace_to_CT_volume(const float4& p,
                                                   const float4& v,
                                                   const int3 nvox,
                                                   const float3 dvox)
{
    float4 out = p;

    float3 CT_size = nvox*dvox;
    if ((p.x > dvox.x && p.x < CT_size.x) &&
        (p.y > dvox.y && p.y < CT_size.y) &&
        (p.z > dvox.z && p.z < CT_size.z))
        return out;

    // 0.001f is to start a fraction of a voxel inside the CT
    // Distances to faces of the CT
    float3 d_1 = (0.0001f*dvox - p)/v;
    float3 d_n = (CT_size - 0.0001f*dvox - p)/v;

    if ((d_1.x < 0.0f && d_n.x < 0.0f) ||
        (d_1.y < 0.0f && d_n.y < 0.0f) ||
        (d_1.z < 0.0f && d_n.z < 0.0f)) {

    } else if ((d_1.x*d_n.x <= 0.0f) &&
              (d_1.y*d_n.y <= 0.0f) &&
              (d_1.z*d_n.z <= 0.0f)) {

    } else {
        float temp = min(d_1.x, d_n.x);
        float alphaMin = -1.0f;
        alphaMin = max(alphaMin, temp);

        temp = min(d_1.y, d_n.y);
        alphaMin = max(alphaMin, temp);

        temp = min(d_1.z, d_n.z);
        alphaMin = max(alphaMin, temp);

        out.x = p.x + v.x*alphaMin;
        out.y = p.y + v.y*alphaMin;
        out.z = p.z + v.z*alphaMin;
    }

    return out;
}

__device__ __host__ float3 ray_trace_to_CT_volume(const float3& p,
                                                   const float3& v,
                                                   const int3 nvox,
                                                   const float3 dvox)
{
    float4 p_4 = make_float4(p, 0);
    float4 v_4 = make_float4(v, 0);
    float4 out = ray_trace_to_CT_volume(p_4, v_4, nvox, dvox);
    return make_float3(out);
}

__device__ __host__ float3 pos_to_int_coordinates(float3 a)
{
    return make_float3(-a.y, -a.x, a.z);
}

__device__ __host__ float4 pos_to_int_coordinates(float4 a)
{
    return make_float4(-a.y, -a.x, a.z, a.w);
}

__device__ __host__ float3 pos_to_ext_coordinates(float3 a)
{
    return make_float3(-a.y, -a.x, a.z);
}

__device__ __host__ float4 pos_to_ext_coordinates(float4 a)
{
    return make_float4(-a.y, -a.x, a.z, a.w);
}

__device__ float3 rotate(const float3& p, const float& gantry, const float& couch)
{
    return make_float3(rotate(make_float4(p), gantry, couch));
}

__device__ float4 rotate(const float4& p, const float& gantry, const float& couch)
{
    float c_couch = __cosf(couch);
    float s_couch = __sinf(couch);
    float c_gantry = __cosf(gantry);
    float s_gantry = __sinf(gantry);

    float4 res;
    res.x = p.x*c_couch - s_couch*(p.y*s_gantry + p.z*c_gantry);
    res.y = p.y*c_gantry - p.z*s_gantry;
    res.z = p.x*s_couch + c_couch*(p.y*s_gantry + p.z*c_gantry);
    res.w = p.w;

    return res;
}
