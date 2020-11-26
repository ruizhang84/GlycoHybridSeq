#ifndef UTIL_MASS_SPECTRUM_H_
#define UTIL_MASS_SPECTRUM_H_

#include <cmath>

namespace util {
namespace mass {

class SpectrumMass
{
public:
    static double Compute(const double mz, const int charge)
    {
        return (mz - kIon) * charge;
    }
    static double ComputeMZ(const double mass, const int charge)
    {
        return (mass + kIon * charge) / charge;
    }
    static double ComputePPM(const double expected, const double observed)
    {
        return std::abs(expected - observed) / expected * 1000000.0;
    }

    static constexpr double kIon = 1.007825;
};

} // namespace mass
} // namespace util

#endif