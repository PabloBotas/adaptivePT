#include "tramp.hpp"
#include "spot.hpp"
#include "energy_range.hpp"


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>

Tramp_t::Tramp_t()
{
}

Tramp_t::Tramp_t(std::string f) : file(f)
{
    read_();
    setEnergies();
}

Tramp_t::~Tramp_t()
{
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
    std::cout << "Reading file " << file << std::endl;
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
        spots.push_back(thisSpot);
        gigaprotons += thisSpot.w;
    }
    nspots = spots.size();

    if( round(1000.*gigaprotons) != round(1000.*gigaprotons_header) )
    {
        std::cerr << "WARNING! Inconsistent tramp file." << std::endl;
        std::cerr << "GigaProtons in the header (" << gigaprotons << ") differs ";
        std::cerr << "from the summed gigaprotons (" << gigaprotons << ")" << std::endl;
    }
    if( nspots != nspots_header )
    {
        std::cerr << "WARNING! Inconsistent tramp file. ";
        std::cerr << "Number of spots in the header (" << nspots_header << ") differs ";
        std::cerr << "from number of spots (" << nspots << ")" << std::endl;
    }
}

void Tramp_t::print(unsigned int n)
{
    if ( n == 0)
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


void Tramp_t::to_file(std::string f)
{
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
    os << "# gigaproton_total "        << gigaprotons            << '\n';
    os << "# rows_total "              << nspots                 << '\n';
    os << "# E(MeV) X(mm) Y(mm) N(Gp)" << '\n';

    for (size_t i = 0; i < nspots; i++)
    {
        os << spots.at(i) << "\n";
    }
}

void Tramp_t::setEnergies()
{
    energies.resize(nspots);
    for (size_t i = 0; i < nspots; i++)
    {
        energies[i] = spots[i].e;
    }
}

std::vector<double> Tramp_t::getWEPLs()
{
    std::vector<double> wepls(nspots);
    Energy_Range_Calculator_t calculator;
    std::transform(energies.begin(), energies.end(), wepls.begin(), calculator);
    return wepls;
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
