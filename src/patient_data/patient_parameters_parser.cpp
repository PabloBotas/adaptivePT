#include "patient_parameters_parser.hpp"
#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <algorithm>
#include <dirent.h>
#include <cstring>
#include <sys/stat.h>

Patient_Parameters_Parser_t::Patient_Parameters_Parser_t()
    : separator ("=")
    , verbose (false)
    , convertUnits (true)
{
}

Patient_Parameters_Parser_t::Patient_Parameters_Parser_t(std::string file)
    : separator ("=")
    , verbose (false)
    , convertUnits (true)
{
    readFromFile(file);
}

void Patient_Parameters_Parser_t::readFromFile(std::string fileName) {
    input.clear();

    std::ifstream inputFile(fileName.c_str());
    std::vector<std::string> parameters;

    std::string line;
    while(getline(inputFile, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t comment = line.find("#");
        if (comment != std::string::npos) {
            line = line.substr(0, comment);
        }
        input.push_back(line);
    }
}

//------------------------------------------------------------
// Readers ---------------------------------------------------
//------------------------------------------------------------
template<class T>
T Patient_Parameters_Parser_t::readReal(std::string quantity, T defaultValue)
{
    std::vector<T> returnValues = readVector<T>(quantity, defaultValue, false);
    return returnValues[0];
}

template<class T>
T Patient_Parameters_Parser_t::readInteger(std::string quantity, T defaultValue)
{
    float defaultFloat = defaultValue;
    if (defaultValue == std::numeric_limits<T>::max()) defaultFloat = NAN;
    std::vector<float> returnValues = readVector<float>(quantity, defaultFloat, false);
    return static_cast<T>(returnValues[0]);
}

template<class T>
T Patient_Parameters_Parser_t::readLastRealInVector(std::string quantity, T defaultValue, bool firstIsSize)
{
    std::vector<T> vec = readVector(quantity, defaultValue, firstIsSize);
    return vec.back();
}

template<class T>
T Patient_Parameters_Parser_t::readLastIntInVector(std::string quantity, T defaultValue, bool firstIsSize)
{
    std::vector<T> vec = readVectorInts(quantity, defaultValue, firstIsSize);
    return static_cast<int>(vec[0]);
}

template<class T>
std::vector<T> Patient_Parameters_Parser_t::readVector(std::string quantity, T defaultValue, bool firstIsSize)
{
    // Parameter format: <Type>:<Quantity> = <Values> <Unit>
    std::string quantityLower = trim(utils::toLower(quantity));

    for (size_t line = 0; line < input.size(); line++) {

        size_t separatorPos = input[line].find(separator);

        if (separatorPos != std::string::npos) {

            std::string inputQuantity = trim(utils::toLower(input[line].substr(0, separatorPos)));

            size_t typeEnd = inputQuantity.find(":");
            if (typeEnd != std::string::npos) {
                inputQuantity = inputQuantity.substr(typeEnd + 1, std::string::npos);
            }

            if (inputQuantity == quantityLower) {

                std::string setting = input[line].substr(separatorPos + 1, std::string::npos);

                setting = trim(setting);

                std::string values;
                std::string unit;

                size_t unitStart = setting.rfind(" ");
                if (unitStart != std::string::npos) {
                    unit = trim(setting.substr(unitStart, std::string::npos));

                    char *unitStringPtr;
                    strtod(unit.c_str(), &unitStringPtr);

                    if (unit.c_str() == unitStringPtr) {
                        // Not a number, so unit
                        values = trim(setting.substr(0, unitStart));
                    } else {
                        values = trim(setting);
                        unit.clear();
                    }
                } else {
                    values = trim(setting);
                }

                unit = utils::toLower(unit);

                std::vector<T> returnValues;

                while (true) {
                    size_t valueSeparator = values.find(" ");
                    if (valueSeparator == std::string::npos) {
                        returnValues.push_back(atof(values.c_str()));
                        break;
                    } else {
                        std::string value = values.substr(0, valueSeparator);
                        returnValues.push_back(atof(value.c_str()));
                        values = values.substr(valueSeparator, std::string::npos);
                        values = trim(values);
                        if (values.empty()) break;
                    }
                }

                if (convertUnits) {
                    T multiplier = 1.0;

                    if (unit == "mm") multiplier *= 0.1;
                    if (unit == "deg") multiplier *= M_PI/180.0;

                    int skip = 0;
                    if (firstIsSize) skip = 1;

                    for (size_t i = skip; i < returnValues.size(); i++) {
                        returnValues[i] *= multiplier;
                    }
                }

                if (verbose) {
                    std::cout << quantity << " \t";
                    for (size_t i = 0; i< returnValues.size(); i++) {
                        std::cout << returnValues[i] << " ";
                    }
                    std::cout << std::endl;
                }

                if (firstIsSize)
                    returnValues.erase(returnValues.begin());

                return returnValues;
            }
        }
    }

    if (isnan(defaultValue)) {
        std::cerr << "Cannot find quantity " + quantity + " in parameter list." << std::endl;
        exit(EXIT_FAILURE);
    } else {
        std::vector<T> returnDefault;
        returnDefault.push_back(defaultValue);

        if (verbose) {
            std::cout << quantity << " \t";
            for (size_t i = 0; i< returnDefault.size(); i++) {
                std::cout << returnDefault[i] << " ";
            }
            std::cout << std::endl;
        }

        return returnDefault;
    }

}

std::string Patient_Parameters_Parser_t::readString(std::string quantity, std::string defaultValue)
{
    // Parameter format: <Type>:<Quantity> = <String>
    std::string quantityLower = trim(utils::toLower(quantity));

    for (size_t line = 0; line < input.size(); line++) {

        size_t separatorPos = input[line].find(separator);

        if (separatorPos != std::string::npos) {
            std::string inputQuantity = trim(utils::toLower(input[line].substr(0, separatorPos)));

            size_t typeEnd = inputQuantity.find(":");
            if (typeEnd != std::string::npos) {
                inputQuantity = inputQuantity.substr(typeEnd + 1, std::string::npos);
            }

            if (inputQuantity == quantityLower) {
                std::string setting = input[line].substr(separatorPos + 1, std::string::npos);
                setting = trim(setting);

                // Trim quotes
                if (!setting.empty() && setting[0] == '"') {
                    setting.erase(0, 1);
                }
                if (!setting.empty() && setting[setting.size() - 1] == '"') {
                    setting.erase(setting.size() - 1, 1);
                }

                if (verbose) {
                    std::cout << quantity << " \t" << setting << std::endl;
                }

                return setting;
            }
        }
    }

    if (defaultValue == "none") {
        std::cerr << "Cannot find quantity " + quantity + " in parameter list." << std::endl;
        exit(EXIT_FAILURE);
    } else {
        if (verbose) {
            std::cout << quantity << " \t" << defaultValue << std::endl;
        }
        return defaultValue;
    }
}

bool Patient_Parameters_Parser_t::readBool(std::string quantity, bool defaultValue)
{
    // Parameter format: <Type>:<Quantity> = <String>
    std::string defaultString = "false";
    if(defaultValue)
        defaultString = "true";
    
    std::string value = readString(quantity, defaultString);
    if(utils::toLower(value) == "true")
        return true;
    else
        return false;
    
}

template<class T>
std::vector<T> Patient_Parameters_Parser_t::readVectorInts(std::string quantity, T defaultValue, bool firstIsSize)
{
    float defaultFloat = defaultValue;
    if (defaultValue == std::numeric_limits<T>::max()) defaultFloat = NAN;

    std::vector<float> returnValues = readVector<float>(quantity, defaultFloat, firstIsSize);

    std::vector<T> returnIntegers;
    for (size_t i = 0; i < returnValues.size(); i++) {
        returnIntegers.push_back(static_cast<T>(returnValues[i]));
    }
    return returnIntegers;
}

//------------------------------------------------------------
// Getters, setters: Query from class ------------------------
//------------------------------------------------------------

std::string Patient_Parameters_Parser_t::getSeparator()
{
    return separator;
}

bool Patient_Parameters_Parser_t::getVerbose()
{
    return verbose;
}

bool Patient_Parameters_Parser_t::getConvertUnits()
{
    return convertUnits;
}

void Patient_Parameters_Parser_t::setSeparator(std::string newValue)
{
    separator = newValue;
}

void Patient_Parameters_Parser_t::setVerbose(bool newValue)
{
    verbose = newValue;
}

void Patient_Parameters_Parser_t::setConvertUnits(bool newValue)
{
    convertUnits = newValue;
}


// Explicit Instantiation
template float Patient_Parameters_Parser_t::readReal<float>(std::string, float);
template double Patient_Parameters_Parser_t::readReal<double>(std::string, double);
template unsigned int Patient_Parameters_Parser_t::readReal<unsigned int>(std::string, unsigned int);
template std::vector<float> Patient_Parameters_Parser_t::readVector(std::string, float, bool);
template std::vector<double> Patient_Parameters_Parser_t::readVector(std::string, double, bool);
template std::vector<unsigned int> Patient_Parameters_Parser_t::readVectorInts(std::string, unsigned int, bool);
template float Patient_Parameters_Parser_t::readLastRealInVector<float>(std::string, float, bool);
template unsigned int Patient_Parameters_Parser_t::readInteger<unsigned int>(std::string, unsigned int);
template unsigned int Patient_Parameters_Parser_t::readLastIntInVector<unsigned int>(std::string, unsigned int, bool);

std::vector<std::string> getFoldersWithFile(std::string folderpath, std::string name)
//      return all folder names of all folders
{
    std::vector<std::string> returnStrings; //empty returnvalue

    struct dirent *direntp = NULL;
    DIR *dirp = NULL;

    dirp = opendir(folderpath.c_str());
    if (dirp == NULL)
    {
        std::cerr << "ERROR, " << folderpath << " cannot be openned!" << std::endl;
        std::cerr << "Was searching for " << name  << std::endl;
        perror("Cannot open directory");
        exit(EXIT_FAILURE);
    }

    // For every directory entry...
    while ((direntp = readdir(dirp)) != NULL)
    {
        // Ignore special directories.
        if ((std::strcmp(direntp->d_name, ".") == 0) ||
                (std::strcmp(direntp->d_name, "..") == 0))
            continue;

        struct stat fstat;
        std::string full_name;

        full_name = std::string(folderpath);
        if(*(full_name.end()-1) != '/')
            full_name += std::string("/");
        full_name += std::string(direntp->d_name);

        // Get only if it is really directory.
        if (stat(full_name.c_str(), &fstat) < 0)
            continue;
        if (S_ISDIR(fstat.st_mode))
        {
            // Check name
            if (stat((full_name + "/" + name).c_str() , &fstat) < 0)
                continue;
            returnStrings.push_back(full_name);
        }
    }

    std::sort(returnStrings.begin(), returnStrings.end());

    return returnStrings;
}

std::vector<std::string> getFilesWithSuffix(std::string folderpath, std::string suffix, std::string contains)
//  get a list of files inside of the folderpath with given suffix
{
    std::vector<std::string> returnStrings;

    struct dirent *direntp = NULL;
    DIR *dirp = NULL;

    dirp = opendir(folderpath.c_str());
    if (dirp == NULL)
    {
        std::cerr << "ERROR, " << folderpath << " cannot be openned!" << std::endl;
        std::cerr << "Was searching for " << suffix  << std::endl;
        perror("Cannot open directory");
        exit(EXIT_FAILURE);
    }

    // For every directory entry...
    while ((direntp = readdir(dirp)) != NULL)
    {
        // Ignore special directories.
        if ((std::strcmp(direntp->d_name, ".") == 0) || (std::strcmp(direntp->d_name, "..") == 0))
            continue;

        struct stat fstat;
        std::string full_name;

        //  get full name
        full_name = std::string(folderpath);
        if (*(full_name.end()-1) != '/')
            full_name += std::string("/");
        full_name += std::string(direntp->d_name);

        //  get info for this entry(folder or file)
        if (stat(full_name.c_str(), &fstat) < 0)
            continue;

        //  if a file, check it
        if(!S_ISDIR(fstat.st_mode))
        {
            std::string temp(direntp->d_name);
            int n = temp.size();
            int start = n-suffix.size();
            std::string MacFilesPref = "._";
            if(start<0) continue;  //the current file name is shorter than suffix
            if(suffix == temp.substr(n-suffix.size(),n) && temp.substr(0,2) != MacFilesPref) {
                // Check name
                if (contains.empty() || temp.find(contains) != std::string::npos) {
                    returnStrings.push_back(full_name);
                }
            }
        }
        //  if a folder, recursive
        else
        {
            std::vector<std::string> temp = getFilesWithSuffix(full_name, suffix);
            returnStrings.insert(returnStrings.end(), temp.begin(), temp.end() ); //append
        }
    }

    std::sort(returnStrings.begin(), returnStrings.end());

    return returnStrings;
}
