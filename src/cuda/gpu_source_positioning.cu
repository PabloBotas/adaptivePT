#include "gpu_source_positioning.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_utils.cuh"

void virtual_src_to_treatment_plane(const unsigned int& num,
                                    const std::vector<BeamAngles_t>& angles,
                                    const float3& ct_offsets)
{
    std::vector<float2> temp(angles.size());
    for (size_t i = 0; i < angles.size(); i++)
    {
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

void correct_offsets(const unsigned int& num,
                     const float3& offsets,
                     const float3& original_offsets)
{
    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    correct_offsets_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, offsets, original_offsets);
    check_kernel_execution(__FILE__, __LINE__);
}

__host__ Array4<float>
get_treatment_planes (const std::vector<BeamAngles_t>& angles)
{
    Array4<float> planes(angles.size());
    for (size_t i = 0; i < angles.size(); i++)
    {
        planes.at(i) = rotate(Vector4_t<float>(0, 0, 1, 0),
                              angles.at(i).gantry, angles.at(i).couch);
    }

    return planes;
}

__global__ void virtual_src_to_treatment_plane_kernel(const int num,
                                                      const float2* angles,
                                                      const float3 ct_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
        float4 pos  = xdata[tid];
        float4 vel  = vxdata[tid];
        short2 meta = ixdata[tid]; // x = beam_id, y = spot_id
        short beamid = meta.x;

        // Adjust to internal coordinates
        pos = ext_to_int_coordinates(pos);
        vel = ext_to_int_coordinates(vel);

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

__global__ void correct_offsets_kernel(const int num,
                                       const float3 offsets,
                                       const float3 original_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
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

__global__ void treatment_plane_to_virtual_src_kernel(const int num,
                                                      const float2* angles,
                                                      const float3 ct_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
        float4 pos  = xdata[tid];
        float4 vel  = vxdata[tid];
        short2 meta = ixdata[tid]; // x = beam_id, y = spot_id
        short beamid = meta.x;

        // TODO: Ray trace to plane outside of CT volume if necessary????

        pos.x += ct_offsets.x;
        pos.y += ct_offsets.y;
        pos.z += ct_offsets.z;

        //  rotate location using gantry and couch
        float gantry = angles[beamid].x;
        float couch  = angles[beamid].y;
        pos = rotate(pos, -gantry, -couch);
        vel = rotate(vel, -gantry, -couch);

        // Adjust to external coordinates
        pos = int_to_ext_coordinates(pos);
        vel = int_to_ext_coordinates(vel);

        xdata[tid]  = pos;
        vxdata[tid] = vel;
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

    if((d_1.x < 0.0f && d_n.x < 0.0f) ||
       (d_1.y < 0.0f && d_n.y < 0.0f) ||
       (d_1.z < 0.0f && d_n.z < 0.0f))
    {

    }
    else if((d_1.x*d_n.x <= 0.0f) &&
            (d_1.y*d_n.y <= 0.0f) &&
            (d_1.z*d_n.z <= 0.0f))
    {

    }
    else
    {
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

__device__ __host__ float3 ext_to_int_coordinates(float3 a)
{
    return make_float3(-a.y, -a.x, a.z);
}

__device__ __host__ float4 ext_to_int_coordinates(float4 a)
{
    return make_float4(-a.y, -a.x, a.z, a.w);
}

__device__ __host__ float3 int_to_ext_coordinates(float3 a)
{
    return make_float3(-a.y, -a.x, a.z);
}

__device__ __host__ float4 int_to_ext_coordinates(float4 a)
{
    return make_float4(-a.y, -a.x, a.z, a.w);
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

__host__ Vector4_t<float> rotate(const Vector4_t<float>& p, const float& gantry, const float& couch)
{
    float c_couch  = cos(couch);
    float s_couch  = sin(couch);
    float c_gantry = cos(gantry);
    float s_gantry = sin(gantry);

    Vector4_t<float> res;
    res.x = p.x*c_couch - s_couch*(p.y*s_gantry + p.z*c_gantry);
    res.y = p.y*c_gantry - p.z*s_gantry;
    res.z = p.x*s_couch + c_couch*(p.y*s_gantry + p.z*c_gantry);
    res.w = p.w;

    return res;
}

