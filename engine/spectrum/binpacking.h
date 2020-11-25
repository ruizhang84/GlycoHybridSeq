#ifndef ENGINE_SPECTRUM_BINPACKING_H
#define ENGINE_SPECTRUM_BINPACKING_H
#include <vector>
#include <algorithm>  
#include <cmath>
#include "../../model/spectrum/spectrum.h"

namespace engine {
namespace spectrum {

class BinPacking
{
public:
    BinPacking(model::spectrum::ToleranceBy type, double tol, double lower, double upper):
        type_(type), tolerance_(tol), lower_(lower), upper_(upper){}
    

    std::vector<double> Packing(const model::spectrum::Spectrum& spec)
    {
        if (type_ == model::spectrum::ToleranceBy::PPM)
            return PPMPacking(spec);
        return DaltonPacking(spec);
    }

protected:
    std::vector<double> DaltonPacking(const model::spectrum::Spectrum& spec)
    {
        // allocate vector
        std::vector<double> result;
        const std::vector<model::spectrum::Peak>& peaks = spec.Peaks();
        
        // fill the bucket
        int size = ceil((upper_ - lower_ + 1) / tolerance_);
        result.assign(size, 0);

        // assign the value
        for(const auto& pk : peaks)
        {
            double mz = pk.MZ();
            if (mz < lower_ || mz > upper_)
                continue;
            int index = floor((mz - lower_) / tolerance_);
            result[index] = std::max(result[index], pk.Intensity());
        }

        return result;
    }


    std::vector<double> PPMPacking(const model::spectrum::Spectrum& spec)
    {
        // alocate vector
        std::vector<double> result;
        const std::vector<model::spectrum::Peak>& peaks = spec.Peaks();

        // fill the bucket
        double ratio = 2.0 * tolerance_ / 1000000;
        int size = ceil(log(upper_ / lower_) / log1p(ratio));
        result.assign(size, 0);

        // assign the value
        for(const auto& pk : peaks)
        {
            double mz = pk.MZ();
            if (mz < lower_ || mz > upper_)
                continue;
            int index = ceil(log(mz / lower_) / log1p(ratio));
            result[index] = std::max(result[index], pk.Intensity());
        }
        return result;
    }

    model::spectrum::ToleranceBy type_;
    double tolerance_;
    double lower_;
    double upper_;
};

} // namespace spectrum
} // namespace engine

#endif