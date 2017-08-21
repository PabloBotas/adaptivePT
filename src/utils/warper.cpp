#include "warper.hpp"

#include "program_options.hpp"
#include "special_types.hpp"
#include "utils.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>
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

void Warper_t::apply_to (Array4<double>& endpoints,
                         Array4<double>& init_pos,
                         const CT_Dims_t& ct,
                         Planes_t treatment_plane,
                         const std::vector<BeamAngles_t>& angles,
                         const std::vector<short>& spots_per_field,
                         const Warp_opts_t options)
{
    probe (endpoints, ct);
    apply_position_options(options, vf, spots_per_field);
    if(exp_file)
        write_to_file (endpoints, spots_per_field);

    warp_points (endpoints);
    warp_init_points (init_pos, treatment_plane, spots_per_field, angles);
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


void Warper_t::warp_init_points (Array4<double>& init_pos,
                                 const Planes_t& pln,
                                 const std::vector<short>& spots_per_field,
                                 const std::vector<BeamAngles_t>& angles)
{
    /*        a2     P  The initial position (b) has to be shifted taking into account the source position (O).
     * O ---a--|-->--|  The angle theta describes the triangle OO'd. Once this is obtained, the point ? can be
          \`   | u   |  determined. The initial data is a, b, c, d and the vector describing the plane.
           \ ` |b    |
            \  `     |
             \ | `   |   r_xy is the unitary vector going from x to y.
              \|   ` |
               o     |c
             X? \    |
                 \   | This part is not necessarily coplanar with the top part
                  \  |
                   \ |
                    \|d
     */

    std::cout << "Working warp_init_points" << std::endl;
    project_vf_on_plane (pln, spots_per_field);

    for (size_t i = 0, ibeam = 0; i < vf.size(); i++)
    {
        if (i == (size_t)spots_per_field.at(ibeam))
            ibeam += 1;

        // Vector3_t<double> a = pln.p.at(ibeam); // has not been raytraced to CT volume
        // Vector3_t<double> u = pln.dir.at(ibeam);
        Vector3_t<double> b = init_pos.at(i); // has been raytraced to CT volume
        Vector3_t<double> c = probes.at(i);
        Vector3_t<double> d = probes.at(i) + vf.at(i);

        Vector3_t<double> Oa = pln.source_a.at(ibeam);
        Vector3_t<double> Ob = pln.source_b.at(ibeam);
        // Vector3_t<double> P  = utils::closest_point(u, a, c);
        // Vector3_t<double> a2 = utils::closest_point(u, a, b);

        // Undo rotation and coords adjustment to correct for SAD
        Vector3_t<double> temp;
        temp = d.get_rotated(-angles.at(ibeam).gantry, -angles.at(ibeam).couch);
        Vector3_t<double> d2(-temp.y, -temp.x, temp.z);
        temp = c.get_rotated(-angles.at(ibeam).gantry, -angles.at(ibeam).couch);
        Vector3_t<double> c2(-temp.y, -temp.x, temp.z);
        temp = b.get_rotated(-angles.at(ibeam).gantry, -angles.at(ibeam).couch);
        Vector3_t<double> b2(-temp.y, -temp.x, temp.z);
 
        double X_x = (b2.z-Oa.z)/(c2.z-Oa.z)*(d2.x-c2.x) + b2.x;
        double X_y = (b2.z-Ob.z)/(c2.z-Ob.z)*(d2.y-c2.y) + b2.y;
        Vector3_t<double> X(-X_y, -X_x, b2.z);
        X.rotate(angles.at(ibeam).gantry, angles.at(ibeam).couch);
        // if (i == 0)
        // {
        //     std::cout << "a:  " << a.x << " " << a.y << " " << a.z << std::endl;
        //     std::cout << "b:  " << b.x << " " << b.y << " " << b.z << std::endl;
        //     std::cout << "c:  " << c.x << " " << c.y << " " << c.z << std::endl;
        //     std::cout << "d:  " << d.x << " " << d.y << " " << d.z << std::endl;
        //     std::cout << "u:  " << u.x << " " << u.y << " " << u.z << std::endl;
        //     std::cout << "vf: " << vf.at(i).x << " " << vf.at(i).y << " " << vf.at(i).z << std::endl;
        //     std::cout << "Oa: " << Oa.x << " " << Oa.y << " " << Oa.z << std::endl;
        //     std::cout << "Ob: " << Ob.x << " " << Ob.y << " " << Ob.z << std::endl;
        //     // std::cout << "P:  " << P.x << " " << P.y << " " << P.z << std::endl;
        //     // std::cout << "a2: " << a2.x << " " << a2.y << " " << a2.z << std::endl;
        //     std::cout << "b2: " << b2.x << " " << b2.y << " " << b2.z << std::endl;
        //     std::cout << "c2: " << c2.x << " " << c2.y << " " << c2.z << std::endl;
        //     std::cout << "d2: " << d2.x << " " << d2.y << " " << d2.z << std::endl;
        //     std::cout << "X:  " << X.x << " " << X.y << " " << X.z << std::endl;
        //     std::cout << "X_wrong: " << b.x+vf.at(i).x << " " << b.y+vf.at(i).y << " " << b.z+vf.at(i).z << std::endl << std::endl;
        // }

        init_pos.at(i) = X;
        // init_pos.at(i) += vf.at(i);
    }
}


void Warper_t::project_vf_on_plane (const Planes_t& pln,
                                    const std::vector<short>& spots_per_field)
{
    for (size_t i = 0, ibeam = 0; i < vf.size(); i++)
    {
        if (i == (size_t)spots_per_field.at(ibeam))
            ibeam += 1;
        double inner = vf.at(i).x*pln.dir.at(ibeam).x +
                      vf.at(i).y*pln.dir.at(ibeam).y +
                      vf.at(i).z*pln.dir.at(ibeam).z;
        double mod_squared = pln.dir.at(ibeam).x*pln.dir.at(ibeam).x +
                            pln.dir.at(ibeam).y*pln.dir.at(ibeam).y +
                            pln.dir.at(ibeam).z*pln.dir.at(ibeam).z;
        double norm = inner/mod_squared;

        vf.at(i).x -= norm*pln.dir.at(ibeam).x;
        vf.at(i).y -= norm*pln.dir.at(ibeam).y;
        vf.at(i).z -= norm*pln.dir.at(ibeam).z;
    }
}


void Warper_t::warp_points (Array4<double>& p)
{
    for (size_t i = 0; i < p.size(); i++)
    {
        p.at(i).x += vf.at(i).x;
        p.at(i).y += vf.at(i).y;
        p.at(i).z += vf.at(i).z;
    }
}

void Warper_t::write_to_file(const Array4<double>& p,
                             const std::vector<short>& spots_per_field)
{
    std::string dir = output.substr(0, output.find_last_of('/'));
    mkdir(dir.c_str(), 0774);
    std::ofstream ofs;
    ofs.open (output, std::ios::out | std::ios::binary);
    utils::check_fs(ofs, output, "to write vector field.");
    std::cout << "Writting probed VF to " << output << " ... ";

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
    std::cout << "done!" << std::endl;
}

void Warper_t::set_probes (const Array4<double>& p)
{
    probes.resize(p.size());
    for (size_t i = 0; i < p.size(); i++)
    {
        probes.at(i).x = p.at(i).x;
        probes.at(i).y = p.at(i).y;
        probes.at(i).z = p.at(i).z;
    }
}

void Warper_t::probe (const Array4<double>& p, const CT_Dims_t& ct)
{
    set_probes(p);

    std::cout << "Probing vector field ..." << std::endl;

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

    // plastimatch probe --location "0 0 0; 0.5 0.5 0.5; 1 1 1" infile.nrrd
    // Build command
    std::string str;
    size_t elements = probes.size();
    const size_t locs_per_batch = 3000;
    size_t nbatches = size_t(elements/locs_per_batch)+1;
    for (size_t batch = 0; batch < nbatches; batch++)
    {
        std::string locations;
        for (size_t i = batch*locs_per_batch; i < elements && i < (batch+1)*locs_per_batch; i++)
        {
            Vector4_t<double> temp;
            temp.x = ct.n.x*ct.d.x - probes.at(i).x + origins.x;
            temp.y = probes.at(i).y + origins.y;
            temp.z = probes.at(i).z + origins.z;
            locations += to_location_str(temp, i+1 == elements && i+1 == (batch+1)*locs_per_batch);
        }
        std::string cmd = "plastimatch probe --location \"";
        cmd += locations + "\" ";
        cmd += file + " 2>&1";

        // Run command
        str += utils::run_command(cmd);
    }

    // Get vf from stdout
    std::stringstream ss_per_line(str);
    std::string line;
    vf.reserve(p.size());

    if (!str.empty())
    {
        while(std::getline(ss_per_line, line, '\n'))
        {
            std::istringstream ss_per_space(line);

            double z, y, x;
            while (ss_per_space)
            {
                std::string dummy;
                ss_per_space >> dummy;
                ss_per_space >> dummy >> dummy >> dummy;
                ss_per_space >> dummy >> dummy >> dummy;
                ss_per_space >> z >> y >> x;
            }
            vf.push_back(Vector3_t<double>(x/10, -y/10, -z/10));
        }
    }
}

// Utils ----------------------------
void Warper_t::flip_positions_X (Array4<double>& vec,
                                 const CT_Dims_t dims)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec.at(i).x = dims.n.x*dims.d.x - vec.at(i).x;
    }
}

void Warper_t::flip_direction_X (Array4<double>& vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        vec.at(i).x *= -1;
    }
}


std::string Warper_t::to_location_str (const Vector3_t<double>& p,
                                       const bool last)
{
    std::string s;
    s += std::to_string(p.z) + " " + 
         std::to_string(p.y) + " " + 
         std::to_string(p.x);
    if (!last)
        s += "; ";
    return s;
}
