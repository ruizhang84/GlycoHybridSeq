#ifndef ENGINE_PROTEIN_PROTEIN_PTM_H_
#define ENGINE_PROTEIN_PROTEIN_PTM_H_

#include <vector>
#include <string>
#include <functional>
#include <iostream>

namespace engine {
namespace protein {

class ProteinPTM
{
public:
    static bool ContainsNGlycanSite(const std::string& sequence)
    {
        for (size_t i = 0; i < sequence.length() - 2; i++)
        {
            char s = std::toupper(sequence[i]);
            char nxs = std::toupper(sequence[i + 2]);
            if (s == 'N' && (nxs == 'S' || nxs == 'T')) return true;
        }

        return false;
    }

    static bool ContainsOGlycanSite(const std::string& sequence)
    {
        for (size_t i = 0; i < sequence.length(); i++)
        {
            char s = std::toupper(sequence[i]);
            if (s == 'S' || s == 'T') return true;
        }

        return false;
    }

    static std::vector<int> FindNGlycanSite(const std::string& sequence)
    {
        std::vector<int>  pos;
        for (size_t i = 0; i < sequence.length() - 2; i++)
        {
            char s = std::toupper(sequence[i]);
            char nxs = std::toupper(sequence[i + 2]);
            if (s == 'N' && (nxs == 'S' || nxs == 'T')) pos.push_back(i);
        }

        return pos;
    }

    static std::vector<int> FindOGlycanSite(const std::string& sequence)
    {
        std::vector<int> pos;
        for (size_t i = 0; i < sequence.length(); i++)
        {
            char s = std::toupper(sequence[i]);
            if (s == 'S' || s == 'T') pos.push_back(i);
        }

        return pos;
    }
    
};

} // namespace protein
} // namespace engine

#endif