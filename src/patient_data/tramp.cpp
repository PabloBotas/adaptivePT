#include "tramp.hpp"
#include "spot.hpp"
#include "energy_range.hpp"
#include "vector4.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <sys/stat.h>

#define MM2CM 0.1
#define CM2MM 10

Tramp_t::Tramp_t()
{
}

Tramp_t::Tramp_t(std::string f, std::string machine_) : machine(machine_), internal_energy_set(false), file(f)
{
    defaults();
    read_();
    setEnergies();
    energy_to_internal();
    setWEPLs();
}

Tramp_t::Tramp_t(std::string f) : internal_energy_set(false), file(f)
{
    defaults();
    read_();
    setEnergies();
    setWEPLs();
}

Tramp_t::~Tramp_t()
{
}

void Tramp_t::defaults()
{
    machine = "";
    gigaprotons = 0;
    nspots = 0;
    patient_id = "";
    patient_first_name = "";
    patient_middle_initial = "";
    patient_last_name = "";
    astroid_id = "";
    course_name = "";
    beam_name = "";
    gantry = "";
    couch_rotation = "";
    z = 0;
    zeff = 0;
}

template<class T>
T Tramp_t::getHeaderValue(std::string line)
{
    std::istringstream ss(line);
    std::string dummy;
    T out;
    ss >> dummy >> dummy >> out;
    return out;
}

template<class T>
T Tramp_t::getHeaderValue(std::ifstream &stream)
{
    std::string line;
    std::getline(stream, line);
    T out = getHeaderValue<T>(line);
    return out;
}

void Tramp_t::read_()
{
    // std::cout << "Reading file " << file << std::endl;
    std::ifstream stream(file);
    if (!stream.is_open()) {
        std::cerr << "Can't open file: " << file << std::endl;
        return;
    }

    patient_id                 = getHeaderValue<std::string>(stream);
    patient_first_name         = getHeaderValue<std::string>(stream);
    patient_middle_initial     = getHeaderValue<std::string>(stream);
    patient_last_name          = getHeaderValue<std::string>(stream);
    astroid_id                 = getHeaderValue<std::string>(stream);
    course_name                = getHeaderValue<std::string>(stream);
    beam_name                  = getHeaderValue<std::string>(stream);
    gantry                     = getHeaderValue<std::string>(stream);
    couch_rotation             = getHeaderValue<std::string>(stream);
    double gigaprotons_header  = getHeaderValue<double>(stream);
    unsigned int nspots_header = getHeaderValue<unsigned int>(stream);

    std::string line;
    std::getline(stream, line);
    while ( std::getline(stream, line) )
    {
        if( line.empty() || line.find_first_not_of(' ') == std::string::npos )
            continue;

        Spot_t thisSpot(line);
        thisSpot.x *= MM2CM;
        thisSpot.y *= MM2CM;
        spots.push_back(thisSpot);
        gigaprotons += thisSpot.w;
        nspots++;
    }

    if( nspots != nspots_header )
    {
        std::cerr << "WARNING! Inconsistent tramp file. ";
        std::cerr << "Number of spots in the header (" << nspots_header << ") differs ";
        std::cerr << "from number of spots (" << nspots << ")" << std::endl;
    }
    if( int(1000.*gigaprotons) != int(1000.*gigaprotons_header) )
    {
        std::cerr << "WARNING! Inconsistent tramp file." << std::endl;
        std::cerr << "GigaProtons in the header (" << gigaprotons_header << ") differs ";
        std::cerr << "from the summed gigaprotons (" << gigaprotons << ")" << std::endl;
    }
}

void Tramp_t::read_file_header(std::string f)
{
    std::ifstream stream(f);
    if (!stream.is_open()) {
        std::cerr << "Can't open file: " << f << std::endl;
        return;
    }

    patient_id             = getHeaderValue<std::string>(stream);
    patient_first_name     = getHeaderValue<std::string>(stream);
    patient_middle_initial = getHeaderValue<std::string>(stream);
    patient_last_name      = getHeaderValue<std::string>(stream);
    astroid_id             = getHeaderValue<std::string>(stream);
    course_name            = getHeaderValue<std::string>(stream);
    beam_name              = getHeaderValue<std::string>(stream);
    gantry                 = getHeaderValue<std::string>(stream);
    couch_rotation         = getHeaderValue<std::string>(stream);
    gigaprotons            = getHeaderValue<double>(stream);
    nspots                 = getHeaderValue<unsigned int>(stream);
}

void Tramp_t::print(unsigned int n)
{
    if (n == 0)
        n = nspots;

    print(0,n);
}

void Tramp_t::print(unsigned int n0, unsigned int n1)
{
    assert ( n0 < n1 );

    std::cout << "Energy (MeV)  |  x (cm)  |  y (cm) |  Weight\n";
    for (unsigned int i = n0; i < n1; i++)
    {
        std::cout << spots.at(i) << '\n';
    }
}


void Tramp_t::to_file(std::string f, std::string dir)
{
    if(!dir.empty())
    {
        mkdir(dir.c_str(), 0774);
        f = f.substr(f.find_last_of("/"));
        f = dir + "/" + f;
    }

    std::ofstream os(f);
    os << "# patient_id "              << patient_id             << '\n';
    os << "# patient_first_name "      << patient_first_name     << '\n';
    os << "# patient_middle_initial "  << patient_middle_initial << '\n';
    os << "# patient_last_name "       << patient_last_name      << '\n';
    os << "# astroid_id "              << astroid_id             << '\n';
    os << "# course_name "             << course_name            << '\n';
    os << "# beam_name "               << beam_name              << '\n';
    os << "# gantry "                  << gantry                 << '\n';
    os << "# couch_rotation "          << couch_rotation         << '\n';
    os << std::fixed << std::setprecision(10) << "# gigaproton_total " << gigaprotons << '\n';
    os << "# rows_total "              << nspots                 << '\n';
    os << "# E(MeV) X(mm) Y(mm) N(Gp)" << '\n';

    for (size_t i = 0; i < nspots; i++)
    {
        Spot_t thisspot = spots.at(i);
        thisspot.x *= CM2MM;
        thisspot.y *= CM2MM;
        os << std::fixed << std::setprecision(10) << thisspot << "\n";
    }
}

void Tramp_t::setEnergies()
{
    energies.reserve(nspots);
    for (size_t i = 0; i < nspots; i++)
    {
        energies.push_back(spots[i].e);
    }
}


void Tramp_t::shift_energies(const std::vector<double>& e_, bool units)
{
    if (e_.size() != nspots)
    {
        std::cerr << "Number of energies to shift " << e_.size();
        std::cerr << " is bigger than the number of spots " << nspots;
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }

    double conv = units ? 1/1e6 : 1;
    for (size_t i = 0; i < nspots; i++)
    {
        energies.at(i) += e_.at(i)*conv;
        spots.at(i).e += e_.at(i)*conv;
    }
}

void Tramp_t::set_pos(const std::vector< Vector4_t<double> > p)
{
    if (p.size() != nspots)
    {
        std::cerr << "Number of positions to set " << p.size();
        std::cerr << " is bigger than the number of spots " << nspots;
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < nspots; i++)
    {
        spots.at(i).x = p.at(i).x;
        spots.at(i).y = p.at(i).y;
    }
}

void Tramp_t::setWEPLs()
{
    std::vector<double> temp;
    if(internal_energy_set)
        temp = energies_internal;
    else
        temp = energies;

    wepls.resize(nspots);
    Energy_Range_Calculator_t calculator;
    std::transform(temp.begin(), temp.end(), wepls.begin(), calculator);
}

std::vector<double> Tramp_t::getWEPLs()
{
    if(wepls.empty())
        setWEPLs();
    return wepls;
}

std::vector<float> Tramp_t::get_weights()
{
    std::vector<float> weights(spots.size());
    for (size_t i = 0; i < spots.size(); i++)
    {
        weights.at(i) = spots.at(i).w;
    }
    return weights;
}

std::vector<float> Tramp_t::get_energies()
{
    std::vector<float> temp(spots.size());
    for (size_t i = 0; i < spots.size(); i++)
    {
        temp.at(i) = spots.at(i).e;
    }
    return temp;
}

std::vector<float> Tramp_t::get_ypos()
{
    std::vector<float> temp(spots.size());
    for (size_t i = 0; i < spots.size(); i++)
    {
        temp.at(i) = spots.at(i).y;
    }
    return temp;
}

std::vector<float> Tramp_t::get_xpos()
{
    std::vector<float> temp(spots.size());
    for (size_t i = 0; i < spots.size(); i++)
    {
        temp.at(i) = spots.at(i).x;
    }
    return temp;
}


double InterpTable(double *vector_X, double *vector_Y, double x, int const npoints)
{
    double result;
    int order = 4; // order of the poly
    // Allocate enough space for any table we'd like to read.
    std::vector<double> lambda(npoints);
    // check order of interpolation
    if (order > npoints)
        order = npoints;
    // if x is ouside the vector_X[] interval
    if (x <= vector_X[0])
        return result = vector_Y[0];
    if (x >= vector_X[npoints-1])
        return result = vector_Y[npoints-1];
    // loop to find j so that x[j-1] < x < x[j]
    int j=0;
    while (j < npoints)
    {
        if (vector_X[j] >= x) break;
        j++;
    }
    // shift j to correspond to (npoint-1)th interpolation
    j = j - order/2;
    // if j is ouside of the range [0, ... npoints-1]
    if (j < 0) j=0;
    if (j+order > npoints ) j=npoints-order;
    result = 0.0;
    for (int is = j; is < j+order; is++)
    {
        lambda[is] = 1.0;
        for (int il=j; il < j+order; il++)
        {
            if(il != is) lambda[is] = lambda[is]*(x-vector_X[il])/(vector_X[is]-vector_X[il]);
        }
        result += vector_Y[is]*lambda[is];
    }
    return result;
}

void Tramp_t::energy_to_internal()
{
    energies_internal.reserve(energies.size());
    if(machine.compare("topasmghr4") == 0 ||
       machine.compare("topasmghr5") == 0)
    {
        double R80AstroidTopasb8[3][27] =
        {
            // Energy [MeV]
            { 91.015, 95.489, 101.50, 109.36, 117.30, 122.69,
              127.81, 134.27, 139.28, 146.68, 152.35, 156.37,
              161.56, 166.79, 171.42, 176.88, 181.97, 186.06,
              190.47, 195.37, 199.13, 204.11, 208.37, 212.67,
              216.58, 221.15, 223.58 },
            // R80 astroid [mm]
            { 65.420, 71.310, 79.590, 91.020, 103.06, 111.54,
              119.84, 130.77, 139.50, 152.79, 163.30, 170.89,
              180.91, 191.30, 200.60, 211.82, 222.48, 231.22,
              240.77, 251.57, 259.96, 271.24, 281.02, 291.04,
              300.26, 311.14, 317.04 },
            // R80 TOPASb8 [mm]
            { 65.040, 70.850, 78.980, 90.170, 102.08, 110.52,
              118.77, 129.54, 138.14, 151.25, 161.62, 169.13,
              179.04, 189.24, 198.44, 209.53, 220.06, 228.67,
              238.09, 248.74, 257.02, 268.16, 277.82, 287.70,
              296.79, 307.55, 313.34 }
        };

        for (size_t i = 0; i < energies.size(); i++)
        {
            // converts energy from MGHR4 machine to the one used by the virtual machines
            double range = InterpTable(&(R80AstroidTopasb8[0][0]),
                                       &(R80AstroidTopasb8[1][0]),
                                       energies.at(i), 27);
            double ECorr = InterpTable(&(R80AstroidTopasb8[2][0]),
                                       &(R80AstroidTopasb8[0][0]),
                                       range, 27);
            energies_internal.push_back(ECorr);
        }
    }
    else
    {
        for (size_t i = 0; i < energies.size(); i++)
        {
            energies_internal.push_back(energies[i]);
        }
    }
}

std::string toLower(std::string s) {
    transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

// void Tramp_t::scale(double ratio)
// {
//     if(previouslyScaled)
//         return;
//     for (size_t i = 0; i < nspots; i++)
//     {
//         spots.at(i).w *= ratio;
//     }
// }
