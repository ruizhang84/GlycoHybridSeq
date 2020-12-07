#ifndef ENGINE_SPECTRUM_BINPACKING_H_
#define ENGINE_SPECTRUM_BINPACKING_H_
#include <vector>
#include <algorithm>  
#include <cmath>
#include <unordered_map>
#include "../../model/spectrum/spectrum.h"

namespace engine {
namespace spectrum {


class BinPacking
{
public:
    BinPacking(model::spectrum::ToleranceBy type, double tol, double lower=200.0, double upper=2000):
        type_(type), tolerance_(tol), lower_(lower), upper_(upper){}
    

    std::unordered_map<int, model::spectrum::Peak> Packing(const model::spectrum::Spectrum& spec)
    {
        if (type_ == model::spectrum::ToleranceBy::PPM)
            return PPMPacking(spec);
        return DaltonPacking(spec);
    }

    int Size()
    {
        if (type_ == model::spectrum::ToleranceBy::Dalton)
            return ceil((upper_ - lower_ + 1.0) / tolerance_);
        double ratio = 1.0/(1.0 - tolerance_ / 1000000);
        return  ceil(log(upper_ / lower_) / log(ratio));
    }

protected:
    std::unordered_map<int, model::spectrum::Peak> DaltonPacking(const model::spectrum::Spectrum& spec)
    {
        // allocate vector
        std::unordered_map<int, model::spectrum::Peak> result;
        const std::vector<model::spectrum::Peak>& peaks = spec.Peaks();
 
        // assign the value
        for(const auto& pk : peaks)
        {
            double mz = pk.MZ();
            if (mz < lower_ || mz > upper_)
                continue;
            int index = Index(mz);
            if (result.find(index) == result.end())
            {
                result[index] = pk;
            }
            else if (result[index].Intensity() < pk.Intensity())
            {
                result[index] = pk;
            }
        }

        return result;
    }


    std::unordered_map<int, model::spectrum::Peak> PPMPacking(const model::spectrum::Spectrum& spec)
    {
        // alocate vector
        std::unordered_map<int, model::spectrum::Peak> result;
        const std::vector<model::spectrum::Peak>& peaks = spec.Peaks();

        // assign the value
        for(const auto& pk : peaks)
        {
            double mz = pk.MZ();
            if (mz < lower_ || mz > upper_)
                continue;
            int index = Index(mz);
            if (result.find(index) == result.end())
            {
                result[index] = pk;
            }
            else if (result[index].Intensity() < pk.Intensity())
            {
                result[index] = pk;
            }
        }
        return result;
    }

    int Index(double expect)
    {
        if (type_ == model::spectrum::ToleranceBy::Dalton)
            return floor((expect - lower_) / tolerance_);
        double ratio = 1.0/(1.0 - tolerance_ / 1000000);
        return  floor(log(expect * 1.0/ lower_) / log(ratio));
    }

    model::spectrum::ToleranceBy type_;
    double tolerance_;
    double lower_;
    double upper_;
};

} // namespace spectrum
} // namespace engine

#endif