#ifndef ENGINE_SEARCH_GLYCAN_H_
#define ENGINE_SEARCH_GLYCAN_H_

#include <string>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "../../algorithm/search/search.h"
#include "../../model/glycan/glycan.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/ion.h"
#include "../../util/mass/spectrum.h"
#include "search_helper.h"


namespace engine{
namespace search{


class GlycanSearch
{
public:
    GlycanSearch(std::unique_ptr<algorithm::search::ISearch<model::spectrum::Peak>> searcher):
        searcher_(std::move(searcher)){}

    // peptide seq, glycan*
    void Search(
        const std::vector<model::spectrum::Peak>& peaks, int max_charge, 
        const std::unordered_map<std::string, std::vector<model::glycan::Glycan*>>& candidate)
    {
        // init search engine
        InitSearch(peaks, max_charge);

        

    }
    

protected:
    void InitSearch(const std::vector<model::spectrum::Peak>& peaks, int max_charge)
    {
        std::vector<std::shared_ptr<algorithm::search::Point<model::spectrum::Peak>>> peak_points; 
        for(const auto& peak : peaks)
        {
            for(int charge = 1; charge <= max_charge; charge++)
            {
                double mass = util::mass::SpectrumMass::Compute(peak.MZ(), charge);
                std::shared_ptr<algorithm::search::Point<model::spectrum::Peak>> point 
                    = std::make_shared<algorithm::search::Point<model::spectrum::Peak>>(mass, peak);
                peak_points.push_back(std::move(point));
            }
        }
        searcher_->Init(peak_points);
    }

    std::unique_ptr<algorithm::search::ISearch<model::spectrum::Peak>> searcher_;
};



} // namespace engine
} // namespace search

#endif