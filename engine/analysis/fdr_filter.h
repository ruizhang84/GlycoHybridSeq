#ifndef ENGINE_ANALYSIS_FDR_FILTER_H
#define ENGINE_ANALYSIS_FDR_FILTER_H

#include <vector>
#include <algorithm>
#include <unordered_map>
#include "search_result.h"

namespace engine {
namespace analysis {

class FDRFilter
{
public:
    FDRFilter(double fdr): fdr_(fdr), cutoff_(-1){}

    void Init()
    {
         // init
        cutoff_ = -1;
        if (decoy_.size() == 0 || target_.size() == 0 || 
            (decoy_.size() * 1.0 / (decoy_.size() + target_.size()) < fdr_))   //trivial case
        {
            return;
        }

        std::vector<double> scores;
        for(const auto& it : target_)
        {
            scores.push_back(it.Score());
        }
        for(const auto& it : decoy_)
        {
            scores.push_back(it.Score());
        }
        std::sort(scores.begin(), scores.end());
        std::sort(target_.begin(), target_.end(), ScoreLess);
        std::sort(decoy_.begin(), decoy_.end(), ScoreLess);

        // compare and compute
        int i = 0, j = 0, k = 0;
        int target_size = (int) target_.size();
        int decoy_size = (int) decoy_.size();
        while (k < target_size + decoy_size)
        {
            double score = scores[k];
            while (i < target_size &&
                target_[i].Score() < score)
            {
                i++;
            }
            // decoy score is no less than targets
            while (j < decoy_size && 
                decoy_[j].Score() < score)
            {
                j++;
            }
            // compute fdr rate
            double rate = (decoy_size - j) * 1.0 / (target_size + decoy_size - i - j + 1);
            if (rate <= fdr_)
            {
                cutoff_ = score;
                return;
            }
            else
            {
                k++;
            }
            
        }
        // set max
        cutoff_ = INT64_MAX;
    }

    std::vector<SearchResult> Filter()
    {
        std::vector<SearchResult> res;
        if (cutoff_ < 0) return target_;

        for(const auto& it : target_)
        {
            if (it.Score() >= cutoff_)
            {
                res.push_back(it);
            }
        }
        if (!res.empty())
            std::sort(res.begin(), res.end(), ScanOrder);
        return res;
    }

    void set_data(std::vector<SearchResult>& targets, 
        std::vector<SearchResult>& decoys) 
    { 
        // acquire the best score of the scan
        std::unordered_map<int, double> score_map;
        for(auto& it: targets)
        {
            int scan = it.Scan();
            if (score_map.find(scan) == score_map.end())
            {
                score_map[scan] = it.Score();
            }
            else
            {
                if (score_map[scan] < it.Score())
                {
                    score_map[scan] = it.Score();
                }
            } 
        }
        // when target and decoy are in the same spectrum,
        // the one with higher score is picked.
        for(auto& it: decoys)
        {
            int scan = it.Scan();
            if (score_map.find(scan) == score_map.end())
            {
                score_map[scan] = it.Score();
            }
            else
            {
                if (score_map[scan] < it.Score())
                {
                    score_map[scan] = it.Score();
                }
            } 
        }
        for(const auto& it: targets)
        {
            int scan = it.Scan();
            if (score_map[scan] > it.Score())
                continue;
            target_.push_back(it);
        }
        for(const auto& it: decoys)
        {
            int scan = it.Scan();
            if (score_map[scan] > it.Score())
                continue;
            decoys.push_back(it);
        }
    }


    
    std::vector<SearchResult>& Target() { return target_; }
    std::vector<SearchResult>& Decoy() { return decoy_; }
    double Cutoff() const { return cutoff_; }
   
    void set_cutoff(double cutoff) { cutoff_ = cutoff; }

protected:
    static bool ScoreLess(const SearchResult& r1, SearchResult& r2) 
        { return r1.Score() < r2.Score(); }

    static bool ScanOrder(const SearchResult& r1, const SearchResult& r2)
        { return r1.Scan() < r2.Scan(); }

    double fdr_;
    double cutoff_;
    std::vector<SearchResult> target_;
    std::vector<SearchResult> decoy_;

};

} // namespace analysis
} // namespace engine 

#endif