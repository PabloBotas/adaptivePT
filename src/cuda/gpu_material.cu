#include "gpu_material.cuh"

#include "gpu_errorcheck.cuh"

#include <iostream>
#include <string>
#include <vector>

std::vector<float> readDensityCorrect(std::string fname)
//  read density correction factor
{
    char buffer[200];
    std::cout << "densityCorrect: Reading " << fname << std::endl;
    FILE *fp;
    fp = fopen(fname.c_str(),"r");
    if (fp == NULL)
    {
        std::cout << "Couldn't open file: " << fname << std::endl;
        exit (EXIT_FAILURE);
    }

    if (!fgets(buffer,200,fp)) ioError("densityCorrect");

    int nfactors = 0;
    if (!fscanf(fp, "%d\n", &nfactors)) ioError("densityCorrect");
    if (!fgets(buffer,200,fp)) ioError("densityCorrect");
    
    //  allocate space and read data
    std::vector<float> dcfactor(nfactors);
    for(int i = 0; i < nfactors; i++)
    {
        if (!fscanf(fp, "%f ,", &dcfactor[i])) ioError("densityCorrect");
    }
    fclose(fp);

    return dcfactor;
}

float HU2dens(int huValue)
//    convert HU to dens, in g/cm^3
{
    float temp;
    //    MGH calibration curve
    if(huValue >= -1000 && huValue < -98)
        temp = 0.00121 + 0.001029700665188*(1000.0 + huValue);
    else if(huValue >= -98 && huValue < 15)
        temp = 1.018 + 0.000893*huValue;
    else if(huValue >= 15 && huValue < 23)
        temp = 1.03;
    else if(huValue >= 23 && huValue < 101)
        temp = 1.003 + 0.001169*huValue;
    else if(huValue >= 101 && huValue < 2001)
        temp = 1.017 + 0.000592*huValue;
    else if(huValue >= 2001 && huValue < 2995)
        temp = 2.201 + 0.0005*(-2000.0 + huValue);
    else
        temp = 4.54;
    return temp;
}

