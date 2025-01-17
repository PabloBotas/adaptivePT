#include "warper.hpp"

#include "program_options.hpp"
#include "special_types.hpp"
#include "utils.hpp"
#include "volume_vf_reader.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <vector>

Warper_t::Warper_t ()
{
}

Warper_t::Warper_t (const std::string vf_file,
                    const std::string data_vf_file)
{
    set(vf_file, data_vf_file);
}

void Warper_t::set (const std::string vf_file,
                    const std::string data_vf_file)
{
    file = vf_file;
    output = data_vf_file;
    if (data_vf_file.empty())
        exp_file = false;
    else
        exp_file = true;

    set_vf_origins();
}

void Warper_t::apply_to_plan (Array4<float>& endpoints,
                              Array4<float>& init_pos,
                              const CT_Dims_t& ct,
                              Planes_t treatment_plane,
                              const std::vector<BeamAngles_t>& angles,
                              const std::vector<short>& spots_per_field,
                              const Adapt_constraints_t options)
{
    if (options != Adapt_constraints_t::FIXED &&
        options != Adapt_constraints_t::FIXED_POS) {
        // Normal execution of the adaptation!!
        probe_plastimatch (endpoints, ct);
        set_average();
        apply_position_options(options, spots_per_field);
        if (exp_file)
            write_to_file (endpoints, spots_per_field, angles);

        warp_points (endpoints);
        warp_init_points (init_pos, treatment_plane, spots_per_field, angles);
        // print_vf();
    } else {
        // No shifts are allowed!!
        vf = Array3<float>(endpoints.size(), Vector3_t<float>(0,0,0));
        set_average();
        if (exp_file)
            write_to_file (endpoints, spots_per_field, angles);
    }
}

Array4<float> Warper_t::apply_to_points (const Array4<float>& pos,
                                         const CT_Dims_t& ct,
                                         Vector3_t<float>* ave)
{
    probe_plastimatch (pos, ct);
    Array4<float> newpos(pos.size());
    newpos.assign(pos.begin(), pos.end());
    warp_points (newpos, ave);
    return newpos;
}

void Warper_t::apply_position_options (Adapt_constraints_t options,
                                       const std::vector<short>& spots_per_field)
{
    if (options == Adapt_constraints_t::ISOCENTER_SHIFT ||
        options == Adapt_constraints_t::ISOCENTER_SHIFT_RANGE_SHIFTER ||
        options == Adapt_constraints_t::ISOCENTER_SHIFT_V_RANGE_SHIFTER)
        apply_rigid_positions();
    else if (options == Adapt_constraints_t::FIELD_ISOCENTER_SHIFT ||
             options == Adapt_constraints_t::FIELD_ISOCENTER_SHIFT_RANGE_SHIFTER ||
             options == Adapt_constraints_t::FIELD_ISOCENTER_SHIFT_V_RANGE_SHIFTER)
        apply_rigid_positions_per_beam(spots_per_field);
    
}

void Warper_t::apply_rigid_positions ()
{
    std::cout << "Applying rigid positions!!" << std::endl;
    // Set average
    for (size_t i = 0; i < vf.size(); i++)
        vf.at(i) = vf_ave;
}

void Warper_t::apply_rigid_positions_per_beam (const std::vector<short>& spots_per_field)
{
    std::cout << "Applying rigid positions per field!!" << std::endl;
    // Calculate average per field
    float accu_spots = 0;
    Array3<float> avgs(spots_per_field.size());
    for (size_t ibeam = 0; ibeam < avgs.size(); ibeam++) {
        // std::cout << "BEAM " << ibeam << std::endl;
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++) {
            size_t idx;
            if (ibeam > 0)
                idx = ispot + accu_spots;
            else
                idx = ispot;
            avgs.at(ibeam) += vf.at(idx);
            // std::cout << "Index " << idx << " " << vf.at(idx).x << " " << avgs.at(ibeam).x << std::endl;
        }
        accu_spots += spots_per_field.at(ibeam);
        avgs.at(ibeam) /= spots_per_field.at(ibeam);
        // std::cout << "Average " << avgs.at(ibeam).x << std::endl;
    }

    // Set average per field
    accu_spots = 0;
    for (size_t ibeam = 0; ibeam < avgs.size(); ibeam++) {
        for (short ispot = 0; ispot < spots_per_field.at(ibeam); ispot++) {
            size_t idx;
            if (ibeam > 0)
                idx = ispot+accu_spots;
            else
                idx = ispot;
            vf.at(idx) = avgs.at(ibeam);
        }
        accu_spots += spots_per_field.at(ibeam);
    }
}


void Warper_t::set_vf_origins()
{
    std::string cmd = "plastimatch header " + file;
    std::string str = utils::run_command(cmd);

    std::stringstream ss_per_line(str);
    std::string line;

    if (!str.empty()) {
        while(std::getline(ss_per_line, line, '\n')) {
            if (line.find("Origin") != std::string::npos) {
                std::istringstream ss_per_space(line);
                std::string dummy;
                ss_per_space >> dummy >> dummy;
                ss_per_space >> origins.x >> origins.y >> origins.z;
            } else if (line.find("Size") != std::string::npos) {
                std::istringstream ss_per_space(line);
                std::string dummy;
                ss_per_space >> dummy >> dummy;
                ss_per_space >> vf_dim_n.x >> vf_dim_n.y >> vf_dim_n.z;
            } else if (line.find("Spacing") != std::string::npos) {
                std::istringstream ss_per_space(line);
                std::string dummy;
                ss_per_space >> dummy >> dummy;
                ss_per_space >> vf_dim_d.x >> vf_dim_d.y >> vf_dim_d.z;
                vf_dim_d /= 10;
            }
            
        }
    }
    // Correct units
    origins.x /= 10;
    origins.y /= 10;
    origins.z /= 10;
}


void Warper_t::warp_init_points (Array4<float>& init_pos,
                                 const Planes_t& pln,
                                 const std::vector<short>& spots_per_field,
                                 const std::vector<BeamAngles_t>& angles)
{
    /*        a2     P  The initial position (b) has to be shifted taking into
     * O ---a--|-->--|  account the source position (O).
     *    \`   | u   |  The angle theta describes the triangle OPd. Once this is
     *     \ ` |b    |  obtained, the point X can be determined. The initial data
     *      \  `     |  is a, b, c, d and the vector describing the plane.
     *       \ | `   |  r_xy is the unitary vector going from x to y.
     *        \|   `.|
     *         o    c|
     *       X? \    |
     *           \   |  This part is not necessarily coplanar with the top part
     *            \  |
     *             \ |
     *              \V d
     */

    std::cout << "Working warp_init_points " << std::endl;
    std::vector<short> accu(spots_per_field.size());
    std::partial_sum(spots_per_field.begin(), spots_per_field.end(), accu.begin());
    project_vf_on_plane (pln, spots_per_field, accu);

    for (size_t i = 0, ibeam = 0; i < vf_planes.size(); i++) {
        if (i == (size_t)accu.at(ibeam))
            ibeam += 1;

        // Vector3_t<float> a = pln.p.at(ibeam); // has not been raytraced to CT volume
        // Vector3_t<float> u = pln.dir.at(ibeam);
        Vector3_t<float> b = init_pos.at(i); // has been raytraced to CT volume
        Vector3_t<float> c = probes.at(i);
        Vector3_t<float> d = probes.at(i) + vf_planes.at(i);

        Vector3_t<float> Ox = pln.source_a.at(ibeam);
        Vector3_t<float> Oy = pln.source_b.at(ibeam);
        // Vector3_t<float> P  = utils::closest_point(u, a, c);
        // Vector3_t<float> a2 = utils::closest_point(u, a, b);

        // Offset with VF average (isocenter shift)
        b += vf_ave_planes.at(ibeam);
        c += vf_ave_planes.at(ibeam);
        // d += vf_ave_planes.at(ibeam);
        Ox += vf_ave_planes.at(ibeam);
        Oy += vf_ave_planes.at(ibeam);

        // Undo rotation and coords adjustment to correct for SAD
        Vector3_t<float> temp;
        temp = d.get_rotated(-angles.at(ibeam).gantry, -angles.at(ibeam).couch);
        Vector3_t<float> d2(-temp.y, -temp.x, temp.z);
        temp = c.get_rotated(-angles.at(ibeam).gantry, -angles.at(ibeam).couch);
        Vector3_t<float> c2(-temp.y, -temp.x, temp.z);
        temp = b.get_rotated(-angles.at(ibeam).gantry, -angles.at(ibeam).couch);
        Vector3_t<float> b2(-temp.y, -temp.x, temp.z);

        float X_x = (b2.z-Ox.z)/(c2.z-Ox.z)*(d2.x-c2.x) + b2.x;
        float X_y = (b2.z-Oy.z)/(c2.z-Oy.z)*(d2.y-c2.y) + b2.y;
        Vector3_t<float> X(-X_y, -X_x, b2.z);
        X.rotate(angles.at(ibeam).gantry, angles.at(ibeam).couch);
        
        // if (i == 0) {
        //     // std::cout << "a:  " << a.x << " " << a.y << " " << a.z << std::endl;
        //     std::cout << "b:  " << b.x << " " << b.y << " " << b.z << std::endl;
        //     std::cout << "c:  " << c.x << " " << c.y << " " << c.z << std::endl;
        //     std::cout << "d:  " << d.x << " " << d.y << " " << d.z << std::endl;
        //     // std::cout << "u:  " << u.x << " " << u.y << " " << u.z << std::endl;
        //     std::cout << "planes: " << vf_planes.at(ibeam).x << " " << vf_planes.at(ibeam).y << " " << vf_planes.at(ibeam).z << std::endl;
        //     std::cout << "vf: " << vf_ave_planes.at(ibeam).x << " " << vf_ave_planes.at(ibeam).y << " " << vf_ave_planes.at(ibeam).z << std::endl;
        //     std::cout << "Ox: " << Ox.x << " " << Ox.y << " " << Ox.z << std::endl;
        //     std::cout << "Oy: " << Oy.x << " " << Oy.y << " " << Oy.z << std::endl;
        //     // std::cout << "P:  " << P.x << " " << P.y << " " << P.z << std::endl;
        //     // std::cout << "a2: " << a2.x << " " << a2.y << " " << a2.z << std::endl;
        //     std::cout << "b2: " << b2.x << " " << b2.y << " " << b2.z << std::endl;
        //     std::cout << "c2: " << c2.x << " " << c2.y << " " << c2.z << std::endl;
        //     std::cout << "d2: " << d2.x << " " << d2.y << " " << d2.z << std::endl;
        //     std::cout << "X:  " << X.x << " " << X.y << " " << X.z << std::endl;
        // }

        init_pos.at(i) = Vector4_t<float>(X, init_pos.at(i).w);
    }
}


void Warper_t::project_vf_on_plane (const Planes_t& pln,
                                    const std::vector<short>& spots_per_field,
                                    const std::vector<short>& accu)
{
    vf_planes.resize(vf.size());
    vf_ave_planes.resize(spots_per_field.size());

    for (size_t i = 0, ibeam = 0; i < vf.size(); i++) {
        if (i == (size_t)accu.at(ibeam))
            ibeam += 1;
        float inner = vf.at(i).x*pln.dir.at(ibeam).x +
                       vf.at(i).y*pln.dir.at(ibeam).y +
                       vf.at(i).z*pln.dir.at(ibeam).z;
        float mod_squared = pln.dir.at(ibeam).x*pln.dir.at(ibeam).x +
                             pln.dir.at(ibeam).y*pln.dir.at(ibeam).y +
                             pln.dir.at(ibeam).z*pln.dir.at(ibeam).z;
        float norm = inner/mod_squared;

        vf_planes.at(i).x = vf.at(i).x - norm*pln.dir.at(ibeam).x;
        vf_planes.at(i).y = vf.at(i).y - norm*pln.dir.at(ibeam).y;
        vf_planes.at(i).z = vf.at(i).z - norm*pln.dir.at(ibeam).z;

        vf_ave_planes.at(ibeam) += vf_planes.at(i)/spots_per_field.at(ibeam);
    }
}


void Warper_t::warp_points (Array4<float>& p, Vector3_t<float>* ave)
{
    for (size_t i = 0; i < p.size(); i++) {
        p.at(i).x += vf.at(i).x;
        p.at(i).y += vf.at(i).y;
        p.at(i).z += vf.at(i).z;
        if (ave) {
            ave->x += vf.at(i).x;
            ave->y += vf.at(i).y;
            ave->z += vf.at(i).z;
        }
    }
    if (ave) {
        ave->x /= p.size();
        ave->y /= p.size();
        ave->z /= p.size();
    }
}

void Warper_t::write_to_file(const Array4<float>& p,
                             const std::vector<short>& spots_per_field,
                             const std::vector<BeamAngles_t>& angles)
{
    std::string dir = output.substr(0, output.find_last_of('/'));
    mkdir(dir.c_str(), 0774);
    std::ofstream ofs;
    ofs.open (output, std::ios::out | std::ios::binary);
    utils::check_fs(ofs, output, "to write vector field.");
    std::cout << "Writting probed VF to " << output << " ... ";

    ofs << "vx vy vz x y z beamid spotid gantryangle couchangle\n";
    uint accu_spot = 0;
    for (uint ibeam = 0; ibeam < spots_per_field.size(); ++ibeam) {
        for (short i=0; i < spots_per_field.at(ibeam); ++i) {
            uint idx = i + accu_spot;
            ofs << vf.at(idx).x << " " << vf.at(idx).y << " " << vf.at(idx).z << " ";
            ofs << p.at(idx).x + origins.x << " " << p.at(idx).y + origins.y;
            ofs << " " << p.at(idx).z + origins.z << " ";
            ofs << ibeam << " " << i << " " << angles.at(ibeam).gantry << " ";
            ofs << angles.at(ibeam).couch << "\n";
        }
        accu_spot += spots_per_field.at(ibeam);
    }
    ofs.close();
    std::cout << "done!" << std::endl;
}

void Warper_t::set_probes (const Array4<float>& p)
{
    probes.resize(p.size());
    for (size_t i = 0; i < p.size(); i++) {
        probes.at(i).x = p.at(i).x;
        probes.at(i).y = p.at(i).y;
        probes.at(i).z = p.at(i).z;
    }
}

void Warper_t::probe_plastimatch (const Array4<float>& p, const CT_Dims_t& ct)
{
    set_probes(p);

    std::cout << "Probing vector field ..." << std::endl;
    // // TEMP!!
    // vf.reserve(p.size());
    // Vector3_t<float> v(0, 0, -2);
    // for (size_t i = 0; i < p.size(); ++i)
    //     vf.push_back(v);
    // // TEMP!!

    // Get vector field image
    if ( ! (utils::ends_with_string(file, ".mha") || 
            utils::ends_with_string(file, ".mhd")) ) {
        std::string ext = utils::get_file_extension(file);
        std::string trans_cmd;
        trans_cmd = "plastimatch xf-convert --input " + file;
        file = utils::replace_string(file, ext, "mha");
        trans_cmd += " --output " + file +
                     " --output-type vf";

        std::string temp = utils::run_command(trans_cmd);
        std::cout << temp << std::endl;
    }

    // Verify VF binning
    if (vf_dim_n.z != ct.n.x || vf_dim_n.y != ct.n.y || vf_dim_n.x != ct.n.z) {
        std::cerr << "ERROR! Wrong VF binning in file: " << file << std::endl;
        std::cerr << "VF.x " << vf_dim_n.x << " CT.x " << ct.n.x << std::endl;
        std::cerr << "VF.y " << vf_dim_n.y << " CT.y " << ct.n.y << std::endl;
        std::cerr << "VF.z " << vf_dim_n.z << " CT.z " << ct.n.z << std::endl << std::endl;
        std::cerr << "The current (bad) implementation probes the VF with voxel indexes. Either "
                     "change the implementation to do it by absolute position or change the VF "
                     "binning" << std::endl;
        exit(EXIT_FAILURE);
    }

    // plastimatch probe --location "0 0 0; 0.5 0.5 0.5; 1 1 1" infile.nrrd
    // Build command
    // TODO: this approach works if VF and CT share dimensions and binning
    std::string str;
    size_t elements = probes.size();
    const size_t locs_per_batch = 3000;
    size_t nbatches = size_t(elements/locs_per_batch)+1;
    for (size_t batch = 0; batch < nbatches; batch++) {
        std::string locations;
        for (size_t i = batch*locs_per_batch; i < elements && i < (batch+1)*locs_per_batch; i++) {
            Vector3_t<uint> idx;
            idx.x = probes.at(i).x/ct.d.x;
            idx.y = probes.at(i).y/ct.d.y;
            idx.z = probes.at(i).z/ct.d.z;
            idx.x = ct.n.x - idx.x - 1;
            std::swap(idx.x, idx.z);
            locations += to_location_str(idx, i+1 == elements && i+1 == (batch+1)*locs_per_batch);
        }
        std::string cmd = "plastimatch probe --index \"";
        cmd += locations + "\" ";
        cmd += file + " 2>&1";

        // Run command
        str += utils::run_command(cmd);
    }

    // Get vf from stdout
    std::stringstream ss_per_line(str);
    std::string line;
    if (vf.size() > 0)
        vf.clear();
    vf.reserve(p.size());

    if (!str.empty()) {
        while(std::getline(ss_per_line, line, '\n')) {
            std::istringstream ss_per_space(line);

            float z, y, x;
            Vector3_t<float> v;
            while (ss_per_space) {
                std::string dummy;
                ss_per_space >> dummy;
                ss_per_space >> dummy >> dummy >> dummy;
                ss_per_space >> dummy >> dummy >> dummy;
                ss_per_space >> x >> y >> z;
            }
            // Change cordinate system and units
            v.x = z/10;
            v.y = -y/10;
            v.z = -x/10;
            // Store
            // printf("%f %f %f || ", v.x, v.y, v.z);
            vf.push_back(v);
        }
        // printf("\n");
    }
}

void Warper_t::probe_internal (const Array4<float>& p, const CT_Dims_t& ct)
{
    std::cout << "Probing vector field ..." << std::endl;
    // Get vector field image
    if ( ! (utils::ends_with_string(file, ".mha") || 
            utils::ends_with_string(file, ".mhd")) ) {
        std::string ext = utils::get_file_extension(file);
        std::string trans_cmd;
        trans_cmd = "plastimatch xf-convert --input " + file;
        file = utils::replace_string(file, ext, "mha");
        trans_cmd += " --output " + file +
                     " --output-type vf";

        std::string temp = utils::run_command(trans_cmd);
        std::cout << temp << std::endl;
    }

    set_probes(p);
    Vf_reader_t vol(file);
    vol.to_int_coordinates();

    if (vf.size() > 0)
        vf.clear();
    vf.reserve(p.size());

    for (size_t i = 0; i < p.size(); ++i) {
        Vector3_t<float> temp;
        temp.x = p.at(i).x + (ct.origin.x - ct.n.x*ct.d.x + ct.d.x/2)
                           - (vol.origin.x - vol.n.x*vol.d.x + vol.d.x/2);
        temp.y = p.at(i).y - (ct.origin.y + ct.n.y*ct.d.y - ct.d.y/2)
                           + (vol.origin.y + vol.n.y*vol.d.y - vol.d.y/2);
        temp.z = p.at(i).z - (ct.origin.z + ct.n.z*ct.d.z - ct.d.z/2)
                           + (vol.origin.z + vol.n.z*vol.d.z - vol.d.z/2);
        Vector3_t<int> vox;
        vox.x = int(floor(temp.x/vol.d.x)+0.5);
        vox.y = int(floor(temp.y/vol.d.y)+0.5);
        vox.z = int(floor(temp.z/vol.d.z)+0.5);
        size_t absvox = vox.z + vox.y*vol.n.z + vox.x*vol.n.z*vol.n.y;
        if (i == 5) {
            Vector3_t<float> temp2 = vol.data.at(absvox);
            std::cerr << "VF: " << temp2.x << " " << temp2.y << " " << temp2.z << std::endl;
        }
        vf.push_back(vol.data.at(absvox));
    }
}

void Warper_t::set_average ()
{
    vf_ave = Vector3_t<float>(0,0,0);
    for (size_t i = 0; i < vf.size(); ++i)
        vf_ave += vf.at(i);
    vf_ave /= vf.size();
}


// Utils ----------------------------
void Warper_t::flip_positions_X (Array4<float>& vec,
                                 const CT_Dims_t dims)
{
    for (size_t i = 0; i < vec.size(); i++)
        vec.at(i).x = dims.n.x*dims.d.x - vec.at(i).x;
}

void Warper_t::flip_direction_X (Array4<float>& vec)
{
    for (size_t i = 0; i < vec.size(); i++)
        vec.at(i).x *= -1;
}


std::string Warper_t::to_location_str (const Vector3_t<uint>& p,
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

void Warper_t::print_vf (unsigned int n)
{
    unsigned int m = vf.size() < n ? vf.size() : n;
    for (unsigned int i = 0; i < m; ++i)
        std::cout << vf.at(i).x << " " << vf.at(i).y << " " << vf.at(i).z << std::endl;
}

