#ifndef __PATIENT_PARAMETERS_PARSER_HPP__
#define __PATIENT_PARAMETERS_PARSER_HPP__

#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <limits>


class Patient_Parameters_Parser_t
{
public:
    Patient_Parameters_Parser_t();
    Patient_Parameters_Parser_t(std::string file);

    void readFromFile(std::string fileName);

    template<class T>
    T              readReal(std::string quantity, T defaultValue = NAN);
    int            readInteger(std::string quantity, int defaultValue = std::numeric_limits<int>::max());
    template<class T>
    std::vector<T> readVector(std::string quantity, T defaultValue = NAN, bool firstIsSize = true);
    template<class T>
    std::vector<T> readVectorInts(std::string quantity, T defaultValue = std::numeric_limits<T>::max(), bool firstIsSize = true);
    std::string    readString(std::string quantity, std::string defaultValue = "none");
    bool           readBool(std::string quantity, bool defaultValue = false);

    std::string getSeparator();
    bool        getVerbose();
    bool        getConvertUnits();

    void setSeparator(std::string newValue);
    void setVerbose(bool newValue);
    void setConvertUnits(bool newValue);

private:
    std::vector<std::string> input;
    std::string separator;
    bool verbose;
    bool convertUnits;

    // trim from start
    static inline std::string ltrim(std::string s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    // trim from end
    static inline std::string rtrim(std::string s) {
        s.erase(std::find_if(s.rbegin(), s.rend(),
                std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    // trim from both ends
    static inline std::string trim(std::string s) {
        return ltrim(rtrim(s));
    }

};

std::vector<std::string> getFilesWithSuffix(std::string folderpath,
                                            std::string suffix,
                                            std::string contains="");
std::vector<std::string> getFoldersWithFile(std::string folderpath, std::string name);


#endif
