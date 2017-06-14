#include "warper.hpp"

#include "special_types.hpp"
#include "utils.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

void
warp_data (std::vector< Vector4_t<float> >& endpoints,
           std::vector< Vector4_t<float> >& init_pos,
           const std::string vf_file,
           const CT_Dims_t& ct,
           std::vector< Vector4_t<float> > treatment_plane)
{
    flip_positions_X(endpoints, ct);
    flip_positions_X(init_pos, ct);
    flip_direction_X(treatment_plane);

    warp_vector_projection(init_pos, treatment_plane, vf_file);
    warp_vector(endpoints, vf_file);

    flip_positions_X(endpoints, ct);
    flip_positions_X(init_pos, ct);
}

void 
warp_vector (std::vector< Vector4_t<float> >& p,
             const std::string vf_file)
{
    std::vector< Vector3_t<float> > vf;
    probe_vf(vf, p, vf_file);
    apply_vf(p, vf);
}

void 
warp_vector_projection (std::vector< Vector4_t<float> >& p,
                        const std::vector< Vector4_t<float> >& treatment_plane,
                        const std::string vf_file)
{
    std::vector< Vector3_t<float> > vf;
    probe_vf(vf, p, vf_file);
    project_vector_on_plane(vf, treatment_plane);
    apply_vf(p, vf);
}

void 
apply_vf (std::vector< Vector4_t<float> >& p,
          const std::vector< Vector3_t<float> >& vf)
{
    for (size_t i = 0; i < vf.size(); i++)
    {
        p.at(i).x += vf.at(i).x;
        p.at(i).y += vf.at(i).y;
        p.at(i).z += vf.at(i).z;
    }
}

void 
get_unitary_vector (std::vector< Vector3_t<float> >& r,
                    const std::vector< Vector4_t<float> >& p,
                    const std::vector< Vector4_t<float> >& p2)
{
    for (size_t i = 0; i < r.size(); i++)
    {
        r.at(i).x = p2.at(i).x - p.at(i).x;
        r.at(i).y = p2.at(i).y - p.at(i).y;
        r.at(i).z = p2.at(i).z - p.at(i).z;
        float norm = sqrt(r.at(i).x*r.at(i).x +
                          r.at(i).y*r.at(i).y +
                          r.at(i).z*r.at(i).z);
        r.at(i).x /= norm;
        r.at(i).y /= norm;
        r.at(i).z /= norm;
    }
}

void 
project_vector_on_plane (std::vector< Vector3_t<float> >& p,
                         const std::vector< Vector4_t<float> >& n)
{
    for (size_t i = 0; i < p.size(); i++)
    {
        float inner = p.at(i).x*n.at(i).x +
                      p.at(i).y*n.at(i).y +
                      p.at(i).z*n.at(i).z;
        float mod_squared = n.at(i).x*n.at(i).x +
                            n.at(i).y*n.at(i).y +
                            n.at(i).z*n.at(i).z;
        float norm = inner/mod_squared;

        p.at(i).x -= norm*n.at(i).x;
        p.at(i).y -= norm*n.at(i).y;
        p.at(i).z -= norm*n.at(i).z;
    }
}

void
probe_vf (std::vector< Vector3_t<float> >& vf,
          const std::vector< Vector4_t<float> >& p,
          std::string vf_file)
{
    //  plastimatch probe --location "0 0 0; 0.5 0.5 0.5; 1 1 1" infile.nrrd
    // Build command
    std::string locations;
    size_t elements = p.size();
    for (size_t i = 0; i < elements; i++)
    {
        locations += to_location_str(p.at(i), i+1 == elements);
    }

    // Get vector field image
    if ( ! (utils::ends_with_string(vf_file, ".mha") || 
            utils::ends_with_string(vf_file, ".mhd")) )
    {
        std::string ext = utils::get_file_extension(vf_file);
        std::string trans_cmd;
        trans_cmd = "plastimatch xf-convert --input " + vf_file;
        vf_file = utils::replace_substring(vf_file, ext, "mha");
        trans_cmd += " --output " + vf_file +
                     " --output-type vf";

        std::string temp = utils::run_command(trans_cmd);
        std::cout << temp << std::endl;
    }

    std::string cmd = "plastimatch probe --location \"";
    cmd += locations + "\" ";
    cmd += vf_file;

    // Run command
    std::string stdout = utils::run_command(cmd);
    // std::cout << stdout << std::endl;

    vf = get_vf_from_stdout(stdout);
}

std::vector< Vector3_t<float> >
probe_vf (const std::vector< Vector4_t<float> >& p,
          std::string vf_file)
{
    std::vector< Vector3_t<float> > vf;
    probe_vf(vf, p, vf_file);
    return vf;
}

std::vector< Vector3_t<float> >
get_vf_from_stdout (std::string str)
{
    std::vector< Vector3_t<float> > v;

    std::stringstream ss_per_line(str);
    std::string line;

    if (!str.empty())
    {
        while(std::getline(ss_per_line, line, '\n'))
        {
            std::istringstream ss_per_space(line);

            Vector3_t<float> this_vf;
            while (ss_per_space)
            {
                std::string dummy;
                ss_per_space >> dummy;
                ss_per_space >> dummy >> dummy >> dummy;
                ss_per_space >> dummy >> dummy >> dummy;
                ss_per_space >> this_vf.x >> this_vf.y >> this_vf.z;
            }
            v.push_back(Vector3_t<float>(this_vf.x, this_vf.y, this_vf.z));
        }
    }

    return v;
}

void
flip_positions_X (std::vector< Vector4_t<float> >& vec,
                  const CT_Dims_t dims)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec.at(i).x = dims.n.x*dims.d.x - vec.at(i).x;
    }
}

void
flip_direction_X (std::vector< Vector4_t<float> >& vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec.at(i).x *= -1;
    }
}


std::string
to_location_str (const Vector3_t<float>& p,
                 const bool last)
{
    std::string s;
    s += " " + std::to_string(p.x) +
         " " + std::to_string(p.y) +
         " " + std::to_string(p.z);
    if (!last)
        s += ";";
    return s;
}
