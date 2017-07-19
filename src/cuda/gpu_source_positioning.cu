#include "gpu_source_positioning.cuh"

#include "initialize_rays.cuh"
#include "gpu_device_globals.cuh"
#include "gpu_geometry_tools.cuh"
#include "gpu_utils.cuh"

void virtual_src_to_treatment_plane(const unsigned int& num,
                                    const std::vector<BeamAngles_t>& angles,
                                    const double3& ct_offsets)
{
    std::vector<double2> temp(angles.size());
    for (size_t i = 0; i < angles.size(); i++)
    {
        temp[i].x = angles.at(i).gantry;
        temp[i].y = angles.at(i).couch;
    }

    double2* angles_gpu;
    array_to_device<double2>(angles_gpu, temp.data(), angles.size());

    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    virtual_src_to_treatment_plane_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, angles_gpu, ct_offsets);
    check_kernel_execution(__FILE__, __LINE__);

    cudaFree(angles_gpu);
}

__global__ void virtual_src_to_treatment_plane_kernel(const int num,
                                                      const double2* angles,
                                                      const double3 ct_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
        double4 pos  = xdata[tid];
        double4 vel  = vxdata[tid];
        short2 meta = ixdata[tid]; // x = beam_id, y = spot_id
        short beamid = meta.x;

        // Adjust to internal coordinates
        pos = ext_to_int_coordinates(pos);
        vel = ext_to_int_coordinates(vel);

        //  rotate location and direction using gantry and couch angles
        double gantry = angles[beamid].x;
        double couch  = angles[beamid].y;
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
                     const double3& offsets,
                     const double3& original_offsets)
{
    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    correct_offsets_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num, offsets, original_offsets);
    check_kernel_execution(__FILE__, __LINE__);
}

__global__ void correct_offsets_kernel(const int num,
                                       const double3 offsets,
                                       const double3 original_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
        double4 pos  = xdata[tid];
        double4 vel  = vxdata[tid];

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

Array4<double> offset_endpoints(const Array4<double>& orig_endpoints,
                                const double3& offsets,
                                const double3& original_offsets)
{
    Array4<double> off_endpoints = orig_endpoints;
    double4* dev_orig_endpoints = NULL;
    array_to_device<double4, Vector4_t<double> >(dev_orig_endpoints, off_endpoints.data(), orig_endpoints.size());

    size_t num = off_endpoints.size();
    int nblocks = 1 + (num-1)/NTHREAD_PER_BLOCK_SOURCE;
    correct_offsets_kernel<<<nblocks, NTHREAD_PER_BLOCK_SOURCE>>>(num,
                             dev_orig_endpoints,
                             offsets,
                             original_offsets);
    check_kernel_execution(__FILE__, __LINE__);

    retrieve_scorer<double, double4>(&off_endpoints[0].x, dev_orig_endpoints, num);

    cudaFree(dev_orig_endpoints);

    return off_endpoints;
}

__global__ void correct_offsets_kernel(const int num,
                                       double4* dev_orig_endpoints,
                                       const double3 offsets,
                                       const double3 original_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
        dev_orig_endpoints[tid].x += original_offsets.x - offsets.x;
        dev_orig_endpoints[tid].y += original_offsets.y - offsets.y;
        dev_orig_endpoints[tid].z += original_offsets.z - offsets.z;

        // printf("%d - 0 - %f %f %f\n", 
               // tid, dev_orig_endpoints[tid].x, dev_orig_endpoints[tid].y, dev_orig_endpoints[tid].z);
    }
}

__host__ void treatment_plane_to_virtual_src(Array4<double>& pos,
                                             const Array4<double>& pos2,
                                             const Patient_Parameters_t& pat)
{
    std::vector<double2> temp(pat.angles.size());
    for (size_t i = 0; i < pat.nbeams; i++)
    {
        temp[i].x = pat.angles.at(i).gantry;
        temp[i].y = pat.angles.at(i).couch;
    }
    double3 offset = make_double3(pat.original_ct.offset.x, pat.original_ct.offset.y, pat.original_ct.offset.z);
    Array4<double> dir = pos2;
    for (size_t i = 0; i < dir.size(); i++)
    {
        // turn dir containing a point to a direction
        dir.at(i).x -= pos.at(i).x;
        dir.at(i).y -= pos.at(i).y;
        dir.at(i).z -= pos.at(i).z;
        double d = sqrt(dir.at(i).x*dir.at(i).x +
                       dir.at(i).y*dir.at(i).y +
                       dir.at(i).z*dir.at(i).z);
        dir.at(i).x /= d;
        dir.at(i).y /= d;
        dir.at(i).z /= d;

        // beam info
        short2 meta = get_beam_spot_id(i, pat.spots_per_field);
        pos.at(i).w = meta.x; // beamid
        dir.at(i).w = meta.y; // spotid
    }

    double2* angles_gpu;
    array_to_device<double2>(angles_gpu, temp.data(), pat.nbeams);

    double4* pos_gpu, *dir_gpu;
    array_to_device<double4>(pos_gpu, pos.data(), pos.size());
    array_to_device<double4>(dir_gpu, dir.data(), dir.size());

    double4* plane_pos_gpu, *plane_dir_gpu;
    array_to_device<double4>(plane_pos_gpu, pat.treatment_planes.p.data(), pat.nbeams);
    array_to_device<double4>(plane_dir_gpu, pat.treatment_planes.dir.data(), pat.nbeams);

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

    retrieve_scorer<double, double4>(&pos[0].x, pos_gpu, pos.size());

    cudaFree(angles_gpu);
    cudaFree(pos_gpu);
    cudaFree(dir_gpu);
    cudaFree(plane_pos_gpu);
    cudaFree(plane_dir_gpu);
}

__global__ void treatment_plane_to_virtual_src_kernel(const int num,
                                                      const int nbeams,
                                                      double4* pos_,
                                                      const double4* dir_,
                                                      const double2* angles,
                                                      const double4* plane_dir,
                                                      const double4* plane_pos,
                                                      const double3 ct_offsets)
//  set source direction
{
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num)
    {
        double4 pos  = pos_[tid];
        double4 vel  = dir_[tid];
        short beamid = pos_[tid].w;
        double4 p_pos = plane_pos[beamid];
        double4 p_dir = plane_dir[beamid];

        // p_pos -= make_double4(ct_offsets);
        // Ray trace to virtual source plane if necessary
        double d = dot((p_pos - pos), p_dir) / dot(vel, p_dir);
        if ( fabs(d) > 0.0001f )
        {
            pos.x += d*vel.x;
            pos.y += d*vel.y;
            pos.z += d*vel.z;
        }

        // Add offsets
        pos.x += ct_offsets.x;
        pos.y += ct_offsets.y;
        pos.z += ct_offsets.z;

        //  rotate location using gantry and couch
        double gantry = angles[beamid].x;
        double couch  = angles[beamid].y;
        pos = rotate(pos, -gantry, -couch);
        // vel = rotate(vel, -gantry, -couch);

        // Adjust to external coordinates
        pos = int_to_ext_coordinates(pos);
        // vel = int_to_ext_coordinates(vel);

        // Return to array
        pos_[tid] = pos;
        // dir_[tid] = vel;
    }
}

__device__ double4 ray_trace_to_CT_volume(const double4& p,
                                          const double4& v)
{
    return ray_trace_to_CT_volume(p, v, ctVox, ctVoxSize);
}

__device__ __host__ double4 ray_trace_to_CT_volume(const double4& p,
                                                   const double4& v,
                                                   const int3 nvox,
                                                   const double3 dvox)
{
    double4 out = p;

    double3 CT_size = nvox*dvox;
    if ((p.x > dvox.x && p.x < CT_size.x) &&
        (p.y > dvox.y && p.y < CT_size.y) &&
        (p.z > dvox.z && p.z < CT_size.z))
        return out;

    // 0.001f is to start a fraction of a voxel inside the CT
    // Distances to faces of the CT
    double3 d_1 = (0.0001f*dvox - p)/v;
    double3 d_n = (CT_size - 0.0001f*dvox - p)/v;

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
        double temp = min(d_1.x, d_n.x);
        double alphaMin = -1.0f;
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

__device__ __host__ double3 ray_trace_to_CT_volume(const double3& p,
                                                   const double3& v,
                                                   const int3 nvox,
                                                   const double3 dvox)
{
    double4 p_4 = make_double4(p, 0);
    double4 v_4 = make_double4(v, 0);
    double4 out = ray_trace_to_CT_volume(p_4, v_4, nvox, dvox);
    return make_double3(out);
}

__device__ __host__ double3 ext_to_int_coordinates(double3 a)
{
    return make_double3(-a.y, -a.x, a.z);
}

__device__ __host__ double4 ext_to_int_coordinates(double4 a)
{
    return make_double4(-a.y, -a.x, a.z, a.w);
}

__device__ __host__ double3 int_to_ext_coordinates(double3 a)
{
    return make_double3(-a.y, -a.x, a.z);
}

__device__ __host__ double4 int_to_ext_coordinates(double4 a)
{
    return make_double4(-a.y, -a.x, a.z, a.w);
}

__device__ double4 rotate(const double4& p, const double& gantry, const double& couch)
{
    double c_couch = __cosf(couch);
    double s_couch = __sinf(couch);
    double c_gantry = __cosf(gantry);
    double s_gantry = __sinf(gantry);

    double4 res;
    res.x = p.x*c_couch - s_couch*(p.y*s_gantry + p.z*c_gantry);
    res.y = p.y*c_gantry - p.z*s_gantry;
    res.z = p.x*s_couch + c_couch*(p.y*s_gantry + p.z*c_gantry);
    res.w = p.w;

    return res;
}

