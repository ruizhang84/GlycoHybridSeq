#ifndef ENGINE_SPECTRUM_SPECTRUM_SIM_H_
#define ENGINE_SPECTRUM_SPECTRUM_SIM_H_

#include <cmath>
#include <unordered_map>
#include <algorithm>
#include "../../model/spectrum/spectrum.h"

namespace engine {
namespace spectrum {

class SpectrumSim
{
public:
    SpectrumSim() = default;

    double ComputeCosine(const std::unordered_map<int, model::spectrum::Peak>& q1,
        const std::unordered_map<int, model::spectrum::Peak>& q2)
    {
        
        double numerator = DotProduct(q1, q2);
        double denominator = sqrt(DotProduct(q1, q1)) * sqrt(DotProduct(q2, q2));
        return numerator / denominator;
    }


protected:
    double DotProduct(const std::unordered_map<int, model::spectrum::Peak>& q1,
        const std::unordered_map<int, model::spectrum::Peak>& q2)
    {
        if (q1.empty() || q2.empty())
            return 0;
        
        double sum = 0;
        for(const auto& it : q1)
        {
            if (q2.find(it.first) != q2.end())
            {
                sum += it.second.Intensity() * q2.find(it.first)->second.Intensity();
            }
        }
        return sum;
    }

};


}  //  namespace calc
}  //  namespace util

#endif