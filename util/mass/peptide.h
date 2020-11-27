#ifndef UTIL_MASS_PEPTIDE_H
#define UTIL_MASS_PEPTIDE_H

#include <string>
#include <cctype>

namespace util {
namespace mass {

class PeptideMass
{
public:
    static double Compute(const std::string& seq)
    {
        double mass = 18.0105;  //water
        for (char s : seq)
        {
            if (std::toupper(s) == 'C')
            {
                //Iodoacetamide
                mass += 57.02146;
            }
            mass += GetAminoAcidMW(s);
        }
        return mass;
    }


protected:
    static double GetAminoAcidMW(const char amino)
    {
        switch (std::toupper(amino))
        {
            case 'A':
                return 71.0371;
            case 'C':
                return 103.00919;
            case 'D':
                return 115.02694;
            case 'E':
                return 129.04259;
            case 'F':
                return 147.06841;
            case 'G':
                return 57.02146;
            case 'H':
                return 137.05891;
            case 'I':
                return 113.08406;
            case 'K':
                return 128.09496; //128.09497
            case 'L':
                return 113.08406;
            case 'M':
                return 131.04049;
            case 'N':
                return 114.04293;
            case 'P':
                return 97.05276;
            case 'Q':
                return 128.05858;
            case 'R':
                return 156.10111; //156.10112
            case 'S':
                return 87.03203;
            case 'T':
                return 101.04768;
            case 'V':
                return 99.06841; //99.06842
            case 'W':
                return 186.07931; //186.07932
            case 'Y':
                return 163.06333;
            default:
                return 118.9;   //Average molecular weight of an amino acid
        }
    }
};

} // namespace mass
} // namespace util

#endif