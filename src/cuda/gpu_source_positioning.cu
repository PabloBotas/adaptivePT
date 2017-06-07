#include "gpu_source_positioning.cuh"

#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_utils.cuh"

void virtual_src_to_treatment_plane(const unsigned int num,
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

        //  rotate location and direction using gantry and couch
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

        // TODO: Ray trace to outside of CT volume if necessary????

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
    float4 out = p;

    float3 CT_size = ctVox*ctVoxSize;
    if ((p.x > ctVoxSize.x && p.x < CT_size.x) &&
        (p.y > ctVoxSize.y && p.y < CT_size.y) &&
        (p.z > ctVoxSize.z && p.z < CT_size.z))
        return out;

    // 0.1f is to start a fraction of a voxel inside the CT
    // Distances to faces of the CT
    float d_1x = (0.1f*ctVoxSize.x - p.x)/v.x;
    float d_1y = (0.1f*ctVoxSize.y - p.y)/v.y;
    float d_1z = (0.1f*ctVoxSize.z - p.z)/v.z;
    float d_nx = (CT_size.x - 0.1f*ctVoxSize.x - p.x)/v.x;
    float d_ny = (CT_size.y - 0.1f*ctVoxSize.y - p.y)/v.y;
    float d_nz = (CT_size.z - 0.1f*ctVoxSize.z - p.z)/v.z;

    if((d_1x < 0.0f && d_nx < 0.0f) ||
       (d_1y < 0.0f && d_ny < 0.0f) ||
       (d_1z < 0.0f && d_nz < 0.0f))
    {

    }
    else if((d_1x*d_nx <= 0.0f) &&
            (d_1y*d_ny <= 0.0f) &&
            (d_1z*d_nz <= 0.0f))
    {

    }
    else
    {
        float temp = min(d_1x, d_nx);
        float alphaMin = -1.0f;
        alphaMin = max(alphaMin, temp);

        temp = min(d_1y, d_ny);
        alphaMin = max(alphaMin, temp);

        temp = min(d_1z, d_nz);
        alphaMin = max(alphaMin, temp);

        out.x = p.x + v.x*alphaMin;
        out.y = p.y + v.y*alphaMin;
        out.z = p.z + v.z*alphaMin;
    }

    return out;
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

