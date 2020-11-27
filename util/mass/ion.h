#ifndef UTIL_MASS_ION_H
#define UTIL_MASS_ION_H

#include<string>
#include "peptide.h"

namespace util {
namespace mass {

enum class IonType { a, b, c, x, y, z};

class IonMass
{
public:
    static double Compute(const std::string& seq, const IonType ion)
    {
        double mass = PeptideMass::Compute(seq); //with an addtional h2o
        switch (ion)
        {
            case IonType::a:
                mass = mass - kOxygen * 2 - kHydrogen * 2 - kCarbon;
                break;
            case IonType::b:
                mass = mass - kOxygen - kHydrogen * 2;
                break;
            case IonType::c:
                mass = mass - kOxygen + kHydrogen + kNitrogen;
                break;
            case IonType::x:
                mass += kCarbon + kOxygen - kHydrogen * 2;
                break;
            case IonType::y:
                break;
            case IonType::z:
                mass = mass - kNitrogen - kHydrogen * 3;
                break;
        }
        return mass;
    }

    static constexpr double kCarbon = 12.0;
    static constexpr double kNitrogen = 14.003074;
    static constexpr double kOxygen = 15.99491463;
    static constexpr double kHydrogen = 1.007825;
};

} // namespace mass
} // namespace util

#endif