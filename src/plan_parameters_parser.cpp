#include "plan_parameters_parser.hpp"

#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>

// Explicit Instantiation
template double Plan_Parameters_Parser_t::readReal<double>(std::string, double);
template std::vector<double> Plan_Parameters_Parser_t::readVector(std::string, double, bool);
template std::vector<unsigned int> Plan_Parameters_Parser_t::readVectorInts(std::string, unsigned int, bool);


Plan_Parameters_Parser_t::Plan_Parameters_Parser_t()
    : separator ("=")
    , verbose (false)
    , convertUnits (true)
{
}

Plan_Parameters_Parser_t::Plan_Parameters_Parser_t(std::string file)
    : separator ("=")
    , verbose (false)
    , convertUnits (true)
{
    readFromFile(file);
}

void Plan_Parameters_Parser_t::readFromFile(std::string fileName) {
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
T Plan_Parameters_Parser_t::readReal(std::string quantity, T defaultValue)
{
    std::vector<T> returnValues = readVector<T>(quantity, defaultValue, false);
    return returnValues[0];
}

int Plan_Parameters_Parser_t::readInteger(std::string quantity, int defaultValue)
{
    float defaultFloat = defaultValue;
    if (defaultValue == std::numeric_limits<int>::max()) defaultFloat = NAN;
    std::vector<float> returnValues = readVector<float>(quantity, defaultFloat, false);
    return static_cast<int>(returnValues[0]);
}

template<class T>
std::vector<T> Plan_Parameters_Parser_t::readVector(std::string quantity, T defaultValue, bool firstIsSize)
{
    // Parameter format: <Type>:<Quantity> = <Values> <Unit>
    std::string quantityLower = trim(toLower(quantity));

    for (size_t line = 0; line < input.size(); line++) {

        size_t separatorPos = input[line].find(separator);

        if (separatorPos != std::string::npos) {

            std::string inputQuantity = trim(toLower(input[line].substr(0, separatorPos)));

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

                unit = toLower(unit);

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

std::string Plan_Parameters_Parser_t::readString(std::string quantity, std::string defaultValue)
{
    // Parameter format: <Type>:<Quantity> = <String>
    std::string quantityLower = trim(toLower(quantity));

    for (size_t line = 0; line < input.size(); line++) {

        size_t separatorPos = input[line].find(separator);

        if (separatorPos != std::string::npos) {
            std::string inputQuantity = trim(toLower(input[line].substr(0, separatorPos)));

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

bool Plan_Parameters_Parser_t::readBool(std::string quantity, bool defaultValue)
{
    // Parameter format: <Type>:<Quantity> = <String>
    std::string defaultString = "false";
    if(defaultValue)
        defaultString = "true";
    
    std::string value = readString(quantity, defaultString);
    if(toLower(value) == "true")
        return true;
    else
        return false;
    
}

template<class T>
std::vector<T> Plan_Parameters_Parser_t::readVectorInts(std::string quantity, T defaultValue, bool firstIsSize)
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

std::string Plan_Parameters_Parser_t::getSeparator()
{
    return separator;
}

bool Plan_Parameters_Parser_t::getVerbose()
{
    return verbose;
}

bool Plan_Parameters_Parser_t::getConvertUnits()
{
    return convertUnits;
}

void Plan_Parameters_Parser_t::setSeparator(std::string newValue)
{
    separator = newValue;
}

void Plan_Parameters_Parser_t::setVerbose(bool newValue)
{
    verbose = newValue;
}

void Plan_Parameters_Parser_t::setConvertUnits(bool newValue)
{
    convertUnits = newValue;
}


