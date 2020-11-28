#ifndef ENGINE_SEARCH_HELPER_H_
#define ENGINE_SEARCH_HELPER_H_

#include <string>
#include <memory>
#include <unordered_map>

#include "../../algorithm/search/search.h"
#include "../../model/glycan/glycan.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/ion.h"


namespace engine{
namespace search{

class SearchHelper
{
public:
    static std::string MakeKeySequence(const std::string& seq, const int pos)
    {
        return seq + "|" + std::to_string(pos);
    }

    static std::pair<std::string, int> ExtractSequence(std::string key)
    {
        std::string seq = key.substr(0, key.find("|"));
        int pos = std::stoi(key.substr(seq.length()+1));
        return std::make_pair(seq, pos);
    }

    // for computing the peptide ions
    static std::vector<double> ComputePTMPeptideMass(const std::string& seq, const int pos)
    {
        std::vector<double> mass_list;
        double mass = util::mass::PeptideMass::Compute(seq.substr(0, pos));
        for (int i = pos; i < (int) seq.length() - 1; i++) // seldom at n
        {
            const char amino = seq[i];
            mass += util::mass::PeptideMass::GetAminoAcidMW(amino);
            mass_list.push_back(util::mass::IonMass::Compute(mass, util::mass::IonType::b));
            mass_list.push_back(util::mass::IonMass::Compute(mass, util::mass::IonType::c));
        }
        mass = util::mass::PeptideMass::Compute(seq.substr(pos+1, seq.length()-pos));
        for (int i = pos; i >= 1; i--)
        {
            const char amino = seq[i];
            mass += util::mass::PeptideMass::GetAminoAcidMW(amino);
            mass_list.push_back(util::mass::IonMass::Compute(mass, util::mass::IonType::y));
            mass_list.push_back(util::mass::IonMass::Compute(mass, util::mass::IonType::z));
        }
        return mass_list;
    }

    static std::vector<double> ComputeNonePTMPeptideMass(const std::string& seq, const int pos)
    {
        std::vector<double> mass_list;
        double mass = 18.0105;
        for (int i = 0; i < pos; i++)
        {
            const char amino = seq[i];
            mass += util::mass::PeptideMass::GetAminoAcidMW(amino);
            mass_list.push_back(util::mass::IonMass::Compute(mass, util::mass::IonType::b));
            mass_list.push_back(util::mass::IonMass::Compute(mass, util::mass::IonType::c));
        }
        mass = 18.0105;
        for (int i = (int) seq.length()-1; i > pos; i--)
        {
            const char amino = seq[i];
            mass += util::mass::PeptideMass::GetAminoAcidMW(amino);
            mass_list.push_back(util::mass::IonMass::Compute(mass, util::mass::IonType::y));
            mass_list.push_back(util::mass::IonMass::Compute(mass, util::mass::IonType::z));
        }
        return mass_list;
    }
};



} // namespace engine
} // namespace search

#endif