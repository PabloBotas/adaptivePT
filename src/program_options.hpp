#ifndef __WARP_OPTIONS__
#define __WARP_OPTIONS__

#include "special_types.hpp"
#include <algorithm>
#include <regex>
#include <string>
#include <vector>

enum class Adapt_constraints_t
{
    ISOCENTER_SHIFT_RANGE_SHIFTER,
    ISOCENTER_SHIFT_V_RANGE_SHIFTER,
    FIELD_ISOCENTER_SHIFT_RANGE_SHIFTER,
    FIELD_ISOCENTER_SHIFT_V_RANGE_SHIFTER,
    ISOCENTER_SHIFT,
    FIELD_ISOCENTER_SHIFT,
    RANGE_SHIFTER,
    V_RANGE_SHIFTER,
    FREE
};

enum class Adapt_methods_t
{
    GEOMETRIC,
    BEAM_MODEL,
    GPMC_DIJ,
    GPMC_DOSE
};

class RShifter_steps_t
{
    float max_thick;
    std::regex numrgx;
public:
    RShifter_steps_t () :
        max_thick(20),
        // Allow only positive numbers here!
        numrgx("((\\+)?[[:digit:]]+)(\\.(([[:digit:]]+)?))?((e|E)((\\+|-)?)[[:digit:]]+)?")
        {};
    ~RShifter_steps_t () {};
    enum options { MGH, CM, HALF_CM, FREE, INTERVAL, LIST };
    options mode;
    std::vector<float> wepl_steps;
    void set_mode (options opt, std::vector<std::string> str = std::vector<std::string>()) {
        mode = opt;
        fill_wepls(str);
    };
    float select_range_shifter (float R, int i=-1) const {
        if (mode == FREE) {
            return R;
        } else {
            if (i > 0) {
                return wepl_steps.at(i-1);
            }
            auto it = std::lower_bound (wepl_steps.begin(), wepl_steps.end(), R);
            return *it;
        }
    };
    int get_wepl_index (float R) const {
        if (mode == FREE) {
            return 0;
        } else {
            auto it = std::find(wepl_steps.begin(), wepl_steps.end(), R);
            return it-wepl_steps.begin();
        }
    };
private:
    void fill_wepls (std::vector<std::string> str = std::vector<std::string>()) {
        try {
            if (mode == MGH) {
                wepl_steps = {3.45, 5.405, 9.2};
            } else if (mode == CM) {
                fill_with_step(1.0);
            } else if (mode == HALF_CM) {
                fill_with_step(0.5);
            } else if (mode == LIST) {
                for (auto i: str) {
                    if (!std::regex_match(i, numrgx))
                        throw std::invalid_argument("Input value for range shifter list of WEPL "
                                                    "thicknesses is not valid: \""+i+"\"");
                    wepl_steps.push_back(std::stof(i));
                }
            } else if (mode == INTERVAL) {
                std::string temp = str.front();
                temp.erase(temp.begin());
                if (!std::regex_match(temp, numrgx))
                        throw std::invalid_argument("Input value for range shifter discrete WEPL "
                                                    "interval is not valid: \""+temp+"\"");
                float step = std::stof(temp);
                if (step == 0.f) {
                    std::cerr << "WARNING! Input interval for discrete range shifter thickness is 0. "
                                 " The range shifter is not discretized!." << std::endl;
                    set_mode (FREE);
                } 
                fill_with_step(step);
            }
        } catch(std::exception& e) {
            std::cerr << "ERROR! " << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
    };
    void fill_with_step(float step) {
        if (step > max_thick) {
            std::cerr << "WARNING! The maximum thickness of the range shifter had to be increased ";
            std::cerr << "from " << max_thick << " cm to " << 5*step << " cm (5 x step). ";
            std::cerr <<  "Are you sure the input data is in cm of water?" << std::endl;
            max_thick = 5*step;
        }
        float thick = step;
        wepl_steps.reserve(int(max_thick/step));
        while (thick <= max_thick) {
            wepl_steps.push_back(thick);
            thick += step;
        }
    };
};

void calculate_range_shifter (const RShifter_steps_t& rshifter_steps,
                              const Adapt_constraints_t& constr,
                              std::vector<float>& new_energies,
                              std::vector<RangeShifter_Dims_t>& rs,
                              std::string& machine,
                              const std::vector<float>& isocenter_to_beam_distance,
                              const std::vector<float>& orig_energies,
                              const std::vector<short>& spots_per_field);

void physical_range_shifter (const RShifter_steps_t& rshifter_steps,
                             std::vector<float>& new_energies,
                             std::vector<RangeShifter_Dims_t>& rs,
                             const std::vector<float>& isocenter_to_beam_distance,
                             const std::vector<float>& orig_energies,
                             const std::vector<short>& spots_per_field);

void virtual_range_shifter (const RShifter_steps_t& rshifter_steps,
                            std::vector<float>& new_energies,
                            const std::vector<float>& orig_energies,
                            const std::vector<short>& spots_per_field);

bool outside_machine_energies (std::vector<float>& below_thres,
                               std::vector<float>& energies,
                               const std::string& machine,
                               const std::vector<short>& spots_per_field);

void correct_energy_range (const RShifter_steps_t& rshifter_steps,
                           std::vector<float>& new_energies,
                           std::vector<RangeShifter_Dims_t>& rs,
                           const std::vector<float>& isocenter_to_beam_distance,
                           const std::vector<short>& spots_per_field,
                           const std::vector<float>& below_thres,
                           const std::string& machine);

#endif