#ifndef ENGINE_SEARCH_GLYCAN_H_
#define ENGINE_SEARCH_GLYCAN_H_

#include <string>
#include <queue> 
#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "../../algorithm/search/search.h"
#include "../../model/glycan/glycan.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/ion.h"
#include "../../util/mass/spectrum.h"
#include "search_glycan_helper.h"


namespace engine{
namespace search{

class GlycanSearch
{
public:
    GlycanSearch(std::unique_ptr<algorithm::search::ISearch<int>> searcher,
        const std::unordered_map<std::string, std::unique_ptr<model::glycan::Glycan>>& glycans_map):
        searcher_(std::move(searcher)), glycans_map_(glycans_map){}

    // peptide seq, glycan*
    std::unordered_map<std::string, std::vector<int>> Search(
        const std::vector<model::spectrum::Peak>& peaks, int max_charge, 
        const std::unordered_map<std::string, std::vector<model::glycan::Glycan*>>& candidates)
    {
        // init search engine
        InitSearch(peaks, max_charge);

        // init peak nodes
        auto queue = InitPriorityQueue(candidates);

        // dp
        auto dp_results = DynamicProgramming(peaks, queue);

        // filter results
        std::unordered_map<std::string, std::vector<int>> results;
        for(const auto& it : dp_results)
        {
            std::pair<std::string, std::string> glypeptide = SearchHelper::ExtractGlycoSequence(it.first);
            model::glycan::Glycan* glycan = glycans_map_.find(glypeptide.first)->second.get();
            for(auto c : candidates.find(glypeptide.first)->second)
            {
                if (SearchHelper::Satisify(c, glycan))
                {
                    results[glypeptide.first + " |" + c->ID()] = it.second;
                }
            }
        }
        
        return results;
    }

    double PeptideMass(const std::string seq)
    {
        if (peptide_mass_.find(seq) == peptide_mass_.end())
        {
            peptide_mass_[seq] = util::mass::PeptideMass::Compute(seq);
        }
        return peptide_mass_[seq];
    }
    

protected:
    std::unordered_map<std::string, std::vector<int>> DynamicProgramming(
        const std::vector<model::spectrum::Peak>& peaks,
        std::priority_queue<PeakNode, std::vector<PeakNode>, PeakNodeComparison>& queue)
    {

        std::unordered_map<std::string, std::vector<int>> results; // glycopeptides, matched peaks
        while (queue.size() > 0)
        {
            // get node
            PeakNode node = queue.top();
            queue.pop();
            // out of mass range
            if (node.Mass() > max_mass_)
            {
                // record a result of node
                for(const auto& it : node.Matches())
                {
                    results[it.first] = std::vector<int>(it.second.begin(), it.second.end());
                }
                continue;
            }

            // match peaks
            double target = node.Mass();
            std::vector<int> matched = searcher_->Search(target);
            // max if matched a peak
            if (matched.size() > 0)
            {
                double best_score = 0;
                std::unordered_map<std::string, std::unordered_set<int>> best_matches;
                for(const auto& it : node.Matches())
                {
                    double score = SearchHelper::PeakScore(peaks, it.second);
                    if (best_score < score)
                    {
                        best_score = score;
                        best_matches.clear();
                        best_matches[it.first] = it.second;
                    }
                    else if (best_score == score)
                    {
                        best_matches[it.first] = it.second;
                    }
                }
                node.set_matches(best_matches);
            }
            // update matches
            if (matched.size() > 0)
            {
                for(auto& it : node.Matches())
                {
                    it.second.insert(matched.begin(), matched.end());
                }
                node.set_miss(0);
            }

            // extending queue
            for(const auto& it : node.Matches())
            {
                std::pair<std::string, std::string> glypeptide = SearchHelper::ExtractGlycoSequence(it.first);
                model::glycan::Glycan* glycan = glycans_map_.find(glypeptide.first)->second.get();

                // record a result of glypeptide
                if (glycan->Children().size() == 0)
                {
                    results[it.first] = std::vector<int>(it.second.begin(), it.second.end());
                }

                for(const auto& g : glycan->Children())
                {
                    double mass = g->Mass() + PeptideMass(glypeptide.second);
                    if (peak_nodes_map_.find(mass) == peak_nodes_map_.end())
                    {
                        PeakNode next;
                        // set mass
                        next.set_mass(mass);
                        // set matches
                        next.set_matches(std::unordered_map<std::string, std::unordered_set<int>>
                            ({ {SearchHelper::MakeKeyGlycoSequence(g->ID(), glypeptide.second), it.second} }));
                        // set missing
                        next.set_miss(node.Missing() + 1);
                        // add node 
                        peak_nodes_map_[mass] = next;
                        // enqueue
                        queue.push(node);
                    }
                    else
                    {
                        PeakNode next = peak_nodes_map_[mass];
                        // set missing
                        next.set_miss(std::min(next.Missing(), node.Missing() + 1));
                        // set matches
                        next.set_matches(std::unordered_map<std::string, std::unordered_set<int>>
                            ({ {SearchHelper::MakeKeyGlycoSequence(g->ID(), glypeptide.second), it.second} }));
                    }
                }
            }
        }
        return results;
    }

    std::priority_queue<PeakNode, std::vector<PeakNode>, PeakNodeComparison> InitPriorityQueue(
        const std::unordered_map<std::string, std::vector<model::glycan::Glycan*>>& candidate)
    {
        std::priority_queue<PeakNode, std::vector<PeakNode>, PeakNodeComparison> queue;
        peak_nodes_map_.clear();
        for(const auto& it : candidate)
        {
            // Y1 mass
            double mass = util::mass::PeptideMass::Compute(it.first) + util::mass::GlycanMass::kHexNAc;
            // node matches
            std::string glycoseq = SearchHelper::MakeKeyGlycoSequence(kY1, it.first);
            if (peak_nodes_map_.find(mass) == peak_nodes_map_.end())
            {
                PeakNode node;
                // set mass
                node.set_mass(mass);
                // set matches
                node.set_matches(std::unordered_map<std::string, std::unordered_set<int>>
                    ({ {glycoseq, std::unordered_set<int>()} }));
                // add node 
                peak_nodes_map_[mass] = node;
                // enqueue
                queue.push(node);
            }
            else
            {
                // update glycopeptide match
                peak_nodes_map_[mass].Add(glycoseq, std::vector<int>());
            }
        }
        return queue;
    }


    void InitSearch(const std::vector<model::spectrum::Peak>& peaks, int max_charge)
    {
        std::vector<std::shared_ptr<algorithm::search::Point<int>>> peak_points; 
        max_mass_ = 2000;
        for(int i = 0; i < (int) peaks.size(); i++)
        {
            const model::spectrum::Peak peak = peaks[i];
            for(int charge = 1; charge <= max_charge; charge++)
            {
                double mass = util::mass::SpectrumMass::Compute(peak.MZ(), charge);
                std::shared_ptr<algorithm::search::Point<int>> point 
                    = std::make_shared<algorithm::search::Point<int>>(mass, i);
                peak_points.push_back(std::move(point));
                // update max mass range
                max_mass_ = std::max(max_mass_, mass);
            }
        }
        searcher_->Init(peak_points);
    }

    std::unique_ptr<algorithm::search::ISearch<int>> searcher_;
    const std::unordered_map<std::string, std::unique_ptr<model::glycan::Glycan>>& glycans_map_;
    double max_mass_;
    std::unordered_map<double, PeakNode> peak_nodes_map_;
    const std::string kY1 = "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ";
    std::unordered_map<std::string, double> peptide_mass_;
};



} // namespace engine
} // namespace search

#endif