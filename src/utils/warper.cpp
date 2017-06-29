#include "warper.hpp"

#include "special_types.hpp"
#include "utils.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

Warper_t::Warper_t (const std::string vf_file,
                    const std::string output_vf)
{
    file = vf_file;
    output = output_vf;
    if (output_vf.empty())
        exp_file = false;
    else
        exp_file = true;

    set_vf_origins();
}

void Warper_t::apply_to (Array4<float>& endpoints,
           Array4<float>& init_pos,
           const CT_Dims_t& ct,
           Array4<float> treatment_plane,
           const std::vector<short>& spots_per_field)
{
    flip_positions_X (endpoints, ct);
    flip_positions_X (init_pos, ct);
    flip_direction_X (treatment_plane);

    probe (endpoints);
    if(exp_file)
        write_to_file (endpoints, spots_per_field);
    warp_points (endpoints);
    project_vf_on_plane (treatment_plane, spots_per_field);
    warp_points (init_pos);

    flip_positions_X(endpoints, ct);
    flip_positions_X(init_pos, ct);
}

void Warper_t::set_vf_origins()
{
    std::string cmd = "plastimatch header " + file;
    std::string str = utils::run_command(cmd);

    std::stringstream ss_per_line(str);
    std::string line;

    if (!str.empty())
    {
        while(std::getline(ss_per_line, line, '\n'))
        {
            if (line.find("Origin") == std::string::npos)
                continue;
            std::istringstream ss_per_space(line);
            while (ss_per_space)
            {
                std::string dummy;
                ss_per_space >> dummy >> dummy;
                ss_per_space >> origins.x >> origins.y >> origins.z;
            }
        }
    }
    // Correct units
    origins.x /= 10;
    origins.y /= 10;
    origins.z /= 10;
}

void Warper_t::project_vf_on_plane (const Array4<float>& n,
                                        const std::vector<short>& spots_per_field)
{
    for (size_t i = 0, ibeam = 0; i < vf.size(); i++)
    {
        if (i == (size_t)spots_per_field.at(ibeam))
            ibeam += 1;
        float inner = vf.at(i).x*n.at(ibeam).x +
                      vf.at(i).y*n.at(ibeam).y +
                      vf.at(i).z*n.at(ibeam).z;
        float mod_squared = n.at(ibeam).x*n.at(ibeam).x +
                            n.at(ibeam).y*n.at(ibeam).y +
                            n.at(ibeam).z*n.at(ibeam).z;
        float norm = inner/mod_squared;

        vf.at(i).x -= norm*n.at(ibeam).x;
        vf.at(i).y -= norm*n.at(ibeam).y;
        vf.at(i).z -= norm*n.at(ibeam).z;
    }
}


void Warper_t::warp_points (Array4<float>& p)
{
    for (size_t i = 0; i < p.size(); i++)
    {
        p.at(i).x += vf.at(i).x;
        p.at(i).y += vf.at(i).y;
        p.at(i).z += vf.at(i).z;
    }
}

void Warper_t::write_to_file(const Array4<float>& p,
                             const std::vector<short>& spots_per_field)
{
    std::ofstream ofs;
    ofs.open (output, std::ios::out | std::ios::binary);
    if( !ofs.is_open() )
    {
        std::cerr << "Can't open file " << output << " to write vector field." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Writting probed VF to " << output << std::endl;

    int beamid = 0;
    ofs << "vx vy vz x y z beamid spotid\n";
    for (short i = 0, spotid = 0; i < (short)vf.size(); i++, spotid++)
    {
        if (i == spots_per_field.at(beamid))
        {
            spotid -= spots_per_field.at(beamid);
            beamid += 1;
        }
        ofs << vf.at(i).x << " " << vf.at(i).y << " " << vf.at(i).z << " ";
        ofs << p.at(i).x + origins.x << " " << p.at(i).y + origins.y;
        ofs << " " << p.at(i).z + origins.z << " ";
        ofs << beamid << " " << spotid << "\n";
    }
    ofs.close();
}

void Warper_t::probe (const Array4<float>& p)
{
    // plastimatch probe --location "0 0 0; 0.5 0.5 0.5; 1 1 1" infile.nrrd
    // Build command
    std::string locations;
    size_t elements = p.size();
    for (size_t i = 0; i < elements; i++)
    {
        Vector4_t<float> temp;
        temp.x = p.at(i).x + origins.x;
        temp.y = p.at(i).y + origins.y;
        temp.z = p.at(i).z + origins.z;
        locations += to_location_str(temp, i+1 == elements);
    }

    // Get vector field image
    if ( ! (utils::ends_with_string(file, ".mha") || 
            utils::ends_with_string(file, ".mhd")) )
    {
        std::string ext = utils::get_file_extension(file);
        std::string trans_cmd;
        trans_cmd = "plastimatch xf-convert --input " + file;
        file = utils::replace_substring(file, ext, "mha");
        trans_cmd += " --output " + file +
                     " --output-type vf";

        std::string temp = utils::run_command(trans_cmd);
        std::cout << temp << std::endl;
    }

    std::string cmd = "plastimatch probe --location \"";
    cmd += locations + "\" ";
    cmd += file;

    // Run command
    std::string str = utils::run_command(cmd);
    // std::cout << str << std::endl;

    // Get vf from stdout
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
            vf.push_back(Vector3_t<float>(this_vf.x/10,
                                          this_vf.y/10, this_vf.z/10));
        }
    }
}

// Utils ----------------------------
void Warper_t::flip_positions_X (Array4<float>& vec,
                                 const CT_Dims_t dims)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec.at(i).x = dims.n.x*dims.d.x - vec.at(i).x;
    }
}

void Warper_t::flip_direction_X (Array4<float>& vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec.at(i).x *= -1;
    }
}


std::string Warper_t::to_location_str (const Vector3_t<float>& p,
                                       const bool last)
{
    std::string s;
    s += std::to_string(p.x) + " " + 
         std::to_string(p.y) + " " + 
         std::to_string(p.z);
    if (!last)
        s += "; ";
    return s;
}
