#ifndef ENGINE_SEARCH_PRECURSOR_MATCH_H
#define ENGINE_SEARCH_PRECURSOR_MATCH_H

#include <string>
#include <vector>
#include <unordered_map>

#include "../../model/glycan/glycan.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/spectrum.h"
#include "../../util/mass/peptide.h"
#include "../../algorithm/search/bucket_search.h"

namespace engine{
namespace search{

class PrecursorMatcher
{
public:
    PrecursorMatcher(model::spectrum::ToleranceBy by, double tol): 
        tolerance_(tol), by_(by),
            searcher_(algorithm::search::BucketSearch<std::string>(by, tol)){}

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
        searcher_.Init(peptides_);
    }

    double Tolerance() const { return tolerance_; }
    model::spectrum::ToleranceBy ToleranceType() const { return by_; }
    void set_tolerance(double tol) { tolerance_ = tol;}
    void set_tolerance_by(model::spectrum::ToleranceBy by) { by_ = by; }

    std::unordered_map<std::string, std::vector<std::string>> Match(double precursor, int charge)
    {
        std::unordered_map<std::string, std::vector<std::string>> results;
        double mass = util::mass::SpectrumMass::Compute(precursor, charge);

        
        for(const auto& glycan : glycans_)
        {
            double target = mass - glycan->Mass();
            if (target <= 0)
                continue;
            std::vector<std::string> peptides = searcher_.Search(target, mass);

            for(const auto& seq : peptides)
            {
                if(results.find(seq) == results.end())
                {
                    results[seq] = std::vector<std::string>();
                }
                results[seq].push_back(glycan->Name());
            }
        }
        return results;
    }

protected:
    double tolerance_;
    model::spectrum::ToleranceBy by_;
    algorithm::search::BucketSearch<std::string> searcher_;
    std::vector<std::shared_ptr<algorithm::search::Point<std::string>>> peptides_;
    std::vector<model::glycan::Glycan*> glycans_;

}; 

} // namespace engine
} // namespace search

#endif