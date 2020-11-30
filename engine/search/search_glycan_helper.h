#ifndef ENGINE_SEARCH_GLYCAN_HELPER_H_
#define ENGINE_SEARCH_GLYCAN_HELPER_H_

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

class PeakNode
{
// peptide -> glycan_id -> list<peaks>
typedef std::unordered_map<std::string, 
        std::unordered_map<std::string, std::unordered_set<int>>> PeakMatch;
public:
    PeakNode() = default;
    
    int Missing() const { return miss_; }
    void set_miss(int miss) { miss_ = miss; }
    double Mass() const { return mass_; }
    void set_mass(double mass) { mass_ = mass; }

    PeakMatch Matches() { return matches_; }
    void set_matches(std::unordered_map<std::string, 
        std::unordered_map<std::string, std::unordered_set<int>>> matches)
        { matches_ = matches; }

    void Add(const std::string& peptide, const std::string& glycan_id, std::vector<int> peaks)
    {
        if(matches_.find(peptide) == matches_.end())
        {
            matches_[peptide] = std::unordered_map<std::string, std::unordered_set<int>>();
        }
        if (matches_[peptide].find(glycan_id) == matches_[peptide].end())
        {
            matches_[peptide].emplace(glycan_id, std::unordered_set<int>());
        }

        matches_[peptide][glycan_id].insert(peaks.begin(), peaks.end());
    }
    void Add(std::vector<int> peaks)
    {
        for(auto& it : matches_)
        {
            for(auto& g: it.second)
            {
                g.second.insert(peaks.begin(), peaks.end());
            }
        }
    }
    void Max(std::vector<model::spectrum::Peak> peaks)
    {
        PeakMatch best;
        for(const auto& it : matches_)
        {
            best[it.first] = std::unordered_map<std::string, std::unordered_set<int>>();
            double best_score = 0;
            for(const auto& g: it.second)
            {
                double score = SearchHelper::ComputePeakScore(peaks, g.second);
                if (best_score < score)
                {
                    best_score = score;
                    best[it.first].clear();
                    best[it.first][g.first] = g.second;
                }
                else if (best_score == score)
                {
                    best[it.first][g.first] = g.second;
                }
            }
        }
        matches_ = best;
    }


protected:
    int miss_ = 1;
    double mass_ = 0;
    PeakMatch matches_;
};

struct PeakNodeComparison
{
    bool operator()(PeakNode* node, PeakNode* other) const
    {
        return node->Mass() > other->Mass();
    }
};


} // namespace engine
} // namespace search

#endif