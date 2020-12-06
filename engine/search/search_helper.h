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

    static std::string MakeKeyGlycoSequence(const std::string glycan_id, const std::string& seq)
    {
        return seq + "|" + glycan_id;
    }

    static std::pair<std::string, int> ExtractSequence(std::string key)
    {
        std::string seq = key.substr(0, key.find("|"));
        int pos = std::stoi(key.substr(seq.length()+1));
        return std::make_pair(seq, pos);
    }

    static std::pair<std::string, std::string> ExtractGlycoSequence(std::string key)
    {
        std::string seq = key.substr(0, key.find("|"));
        std::string glycan_id = key.substr(seq.length()+1);
        return std::make_pair(seq, glycan_id);
    }

    static double ComputePeakScore(const std::vector<model::spectrum::Peak>& peaks, 
        const std::unordered_set<int>& peak_indexes)
    {
        double score = 0;
        for(int index : peak_indexes)
        {
            score += log(peaks[index].Intensity());
        }
        return score;
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
        mass = util::mass::PeptideMass::Compute(seq.substr(pos+1, seq.length()-pos-1));
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

    static double ComputeGlycanMass(const std::vector<model::glycan::Glycan*>& glycans)
    {
        double sums = 0.0;
        int glycan_count = (int) glycans.size();
        for(const auto& glycan : glycans)
        {
            sums += glycan->Mass();
        }
        return sums * 1.0 / glycan_count;
    }

};



} // namespace engine
} // namespace search

#endif