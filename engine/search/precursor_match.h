#ifndef ENGINE_SEARCH_PRECURSOR_MATCH_H_
#define ENGINE_SEARCH_PRECURSOR_MATCH_H_

#include <string>
#include <vector>
#include <unordered_map>

#include "../../model/glycan/glycan.h"
#include "../../util/mass/spectrum.h"
#include "../../util/mass/peptide.h"
#include "../../algorithm/search/search.h"

namespace engine{
namespace search{

class PrecursorMatcher
{
public:
    PrecursorMatcher(std::unique_ptr<algorithm::search::ISearch<std::string>> searcher): 
        searcher_(std::move(searcher)){}

    void Init(const std::vector<std::string>& peptides, 
        const std::unordered_map<std::string, std::unique_ptr<model::glycan::Glycan>>& glycans)
    {
        for(const auto& it : peptides)
        {
            std::shared_ptr<algorithm::search::Point<std::string>> seq = 
                std::make_shared<algorithm::search::Point<std::string>>(util::mass::PeptideMass::Compute(it), it);
            peptides_.push_back(std::move(seq));
        }
        
        for(const auto& it : glycans)
        {
            glycans_.push_back(it.second.get());
        }
        searcher_->Init(peptides_);
    }

    std::unordered_map<std::string, std::vector<model::glycan::Glycan*>> Match(double precursor, int charge)
    {
        std::unordered_map<std::string, std::vector<model::glycan::Glycan*>> results;
        double mass = util::mass::SpectrumMass::Compute(precursor, charge);

        
        for(const auto& glycan : glycans_)
        {
            double target = mass - glycan->Mass();
            if (target <= 0)
                continue;
            std::vector<std::string> peptides = searcher_->Search(target, mass);

            for(const auto& seq : peptides)
            {
                if(results.find(seq) == results.end())
                {
                    results[seq] = std::vector<model::glycan::Glycan*>();
                }
                results[seq].push_back(glycan);
            }
        }
        return results;
    }

protected:
    std::unique_ptr<algorithm::search::ISearch<std::string>> searcher_;
    std::vector<std::shared_ptr<algorithm::search::Point<std::string>>> peptides_;
    std::vector<model::glycan::Glycan*> glycans_;

}; 

} // namespace engine
} // namespace search

#endif