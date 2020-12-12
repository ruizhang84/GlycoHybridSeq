#ifndef ENGINE_SEARCH_GLYCAN_HELPER_H_
#define ENGINE_SEARCH_GLYCAN_HELPER_H_

#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "../../algorithm/search/search.h"
#include "../../model/glycan/glycan.h"
#include "../../model/glycan/nglycan_complex.h"
#include "../../model/glycan/nglycan_hybrid.h"
#include "../../model/glycan/highmannose.h"
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

    bool static IsComplex(const std::string& glycan_id)
        { return (int) glycan_id.size() == 48; }
    bool static IsHybrid(const std::string& glycan_id)
        { return (int) glycan_id.size() == 32; }
    bool static IsHighMannose(const std::string& glycan_id)
        { return (int) glycan_id.size() == 12; }

    PeakMatch MaxBy(std::vector<model::spectrum::Peak> peaks, 
        std::function<bool(const std::string&)> glycanFilter)
    {
        PeakMatch best;
        for(const auto& it : matches_)
        {
            best[it.first] = std::unordered_map<std::string, std::unordered_set<int>>();
            double best_score = 0;
            for(const auto& g: it.second)
            {
                if (!glycanFilter(g.first)) 
                    continue;
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
        return best;
    }

    PeakMatch MaxByHybrid(std::vector<model::spectrum::Peak> peaks)
    {
        PeakMatch best;
        for(const auto& it : matches_)
        {
            best[it.first] = std::unordered_map<std::string, std::unordered_set<int>>();
            std::unordered_map<std::string, std::vector<std::string>> mannose_part; // mannose -> id
            for(const auto& g: it.second)
            {
                if (!IsHybrid(g.first))
                    continue;
                std::string mannose = g.first.substr(8, 4);
                if (mannose_part.find(mannose) == mannose_part.end())
                {
                    mannose_part[mannose] = std::vector<std::string>();
                }
                mannose_part[mannose].push_back(g.first);
            }

            // tricky to avoid cross max on mannose part
            for(const auto& m: mannose_part)
            {
                std::unordered_map<std::string, std::unordered_set<int>> sub_best;
                double best_score = 0;
                for(const auto& g_id: m.second)
                {
                    double score = SearchHelper::ComputePeakScore(peaks, it.second.find(g_id)->second);
                    if (best_score < score)
                    {
                        best_score = score;
                        sub_best.clear();
                        sub_best[g_id] = it.second.find(g_id)->second;
                    }
                    else if (best_score == score)
                    {
                        sub_best[g_id] = it.second.find(g_id)->second;
                    }
                }
                for(const auto& g : sub_best)
                {
                    best[it.first][g.first] = g.second;
                }
            }
        }
        return best;
    }


    void Merge(PeakMatch& to, const PeakMatch& from)
    {
        for(const auto& it : from)
        {
            if(to.find(it.first) == to.end())
            {
                to[it.first] =
                    std::unordered_map<std::string, std::unordered_set<int>>();
            }
            
            for(const auto& g: it.second)
            {
                to[it.first][g.first] = g.second;
            }
        }
    }

    void Max(std::vector<model::spectrum::Peak> peaks)
    {
        PeakMatch best;
        Merge(best, MaxBy(peaks, IsComplex));
        Merge(best, MaxBy(peaks, IsHighMannose));
        Merge(best, MaxByHybrid(peaks));
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