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
#include "../../model/glycan/nglycan_complex.h"
#include "../../model/glycan/nglycan_hybrid.h"
#include "../../model/glycan/highmannose.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/ion.h"
#include "../../util/mass/spectrum.h"
#include "../../util/mass/glycan.h"
#include "search_glycan_helper.h"


namespace engine{
namespace search{

class GlycanSearch
{
public:
    GlycanSearch(std::unique_ptr<algorithm::search::ISearch<int>> searcher,
        const std::unordered_map<std::string, std::unique_ptr<model::glycan::Glycan>>& glycans_map,
        bool complex=true, bool hybrid=false, bool highmannose=false):
        searcher_(std::move(searcher)), glycans_map_(glycans_map), complex_(complex), 
        hybrid_(hybrid), highmannose_(highmannose){}

    // peptide seq, glycan*
    std::unordered_map<std::string, std::unordered_set<int>> Search(
        const std::vector<model::spectrum::Peak>& peaks, int max_charge, 
        const std::unordered_map<std::string, std::vector<model::glycan::Glycan*>>& candidates)
    {
        // init search engine
        InitSearch(peaks, max_charge);

        // init peak nodes
        std::unordered_map<double, std::unique_ptr<PeakNode>> peak_nodes_map;
        std::priority_queue<PeakNode*, std::vector<PeakNode*>, PeakNodeComparison> queue;
        InitPriorityQueue(candidates, peak_nodes_map, queue);

        // dp
        auto dp_results = DynamicProgramming(peaks, peak_nodes_map, queue);

        // filter results
        std::unordered_map<std::string, std::unordered_set<int>> results;
        for(const auto& node : dp_results)
        {
            for(const auto& it : node->Matches())
            {
                // std::string peptide_seq = it.first;
                // std::unordered_map<std::string, std::unordered_set<int>> matched = it.second;
                std::vector<model::glycan::Glycan*> candidate_glycan = candidates.find(it.first)->second;
                for(const auto& g : it.second)
                {
                    for(const auto& glycan : candidate_glycan)
                    {
                        if (Satisify(g.first, glycan))
                        {
                            std::string key = SearchHelper::MakeKeyGlycoSequence(glycan->ID(), it.first);
                            if(results.find(key) == results.end())
                            {
                                results[key] = std::unordered_set<int>();
                            }
                            results[key].insert(g.second.begin(), g.second.end());
                        }
                    }

                }
            }
        }


        return results;
    }

    bool Satisify(const std::string& identified_glycan_id, const model::glycan::Glycan* glycan) const
    {
        const std::vector<int> identified_glycan_table = 
            glycans_map_.find(identified_glycan_id)->second->TableConst();
        const std::vector<int> candidate_glycan_table = glycan->TableConst();
        if (candidate_glycan_table.size() != identified_glycan_table.size())
            return false;
        for(int i = 0; i < (int) identified_glycan_table.size(); i++)
        {
            if (candidate_glycan_table[i] < identified_glycan_table[i])
                return false;
        }

        // check terminal 
        if ((int) identified_glycan_table.size() == 24)
        {
            for(int i = 0; i < 4; i++)
            {
                if ( (identified_glycan_table[12 + i] > 0 || identified_glycan_table[16 + i] > 0
                        ) && identified_glycan_table[4 + i] != candidate_glycan_table[4 + i])
                        return false;
            }
            
        }
        else if ((int) identified_glycan_table.size() == 16)
        {
            for (int i = 0; i < 2; i++)
            {
                if ((identified_glycan_table[10 + i] > 0 || identified_glycan_table[12 + i] > 0
                        ) && identified_glycan_table[6 + i] != candidate_glycan_table[6 + i])
                    return false;
            }

        }
        return true;
    }

    double ComputePeptideMass(const std::string seq)
    {
        if (peptide_mass_.find(seq) == peptide_mass_.end())
        {
            peptide_mass_[seq] = util::mass::PeptideMass::Compute(seq);
        }
        return peptide_mass_[seq];
    }
    

protected:
    std::vector<PeakNode*> DynamicProgramming(
        const std::vector<model::spectrum::Peak>& peaks,
        std::unordered_map<double, std::unique_ptr<PeakNode>>& peak_nodes_map,
        std::priority_queue<PeakNode*, std::vector<PeakNode*>, PeakNodeComparison>& queue)
    {
        std::vector<PeakNode*> matched_nodes;
        while (queue.size() > 0)
        {
            // get node
            PeakNode* node = queue.top();
            queue.pop();

            // // match peaks
            double target = node->Mass();
            std::vector<int> matched = searcher_->Search(target);

            // max if matched a peak
            node->Max(peaks);

            // update matches
            if (matched.size() > 0)
            {
                node->Add(matched);
                node->set_miss(0);
                matched_nodes.push_back(node);
            }
                
            if (node->Missing() > kMissing)
                continue;

            // extending queue
            for(const auto& it : node->Matches())
            {
                std::string peptide = it.first;
                for(const auto& gt : it.second)
                {
                    std::string glycan_id = gt.first;
                    model::glycan::Glycan* glycan = glycans_map_.find(glycan_id)->second.get();

                    std::vector<int> peak_indexes(gt.second.begin(), gt.second.end());
                    for(const auto& g : glycan->Children())
                    {
                        double mass = g->Mass() + util::mass::PeptideMass::Compute(peptide);
                        if (peak_nodes_map.find(mass) == peak_nodes_map.end())
                        {
                            std::unique_ptr<PeakNode> next = 
                                std::make_unique<PeakNode>();
                            // set mass
                            next->set_mass(mass);
                            // set matches
                            next->Add(peptide, g->ID(), peak_indexes);
                            // set missing
                            next->set_miss(node->Missing() + 1);
                            // add node 
                            peak_nodes_map[mass] = std::move(next);
                            // enqueue
                            queue.push(peak_nodes_map[mass].get());
                        }
                        else
                        {
                            // std::cout << "here" << std::endl;
                            PeakNode* next = peak_nodes_map[mass].get();
                            // set missing
                            next->set_miss(std::min(next->Missing(), node->Missing()+1));
                            // set matches
                            next->Add(peptide, g->ID(),  peak_indexes);
                        }
                    }
                }
            }
        }
        return matched_nodes;
    }

    void InitPriorityQueue(
        const std::unordered_map<std::string, std::vector<model::glycan::Glycan*>>& candidate, 
        std::unordered_map<double, std::unique_ptr<PeakNode>>& peak_nodes_map,
        std::priority_queue<PeakNode*, std::vector<PeakNode*>, PeakNodeComparison>& queue)
    {
        for(const auto& it : candidate)
        {
            // Y1 mass
            double mass = ComputePeptideMass(it.first) + util::mass::GlycanMass::kHexNAc;
            // node matches
            if (peak_nodes_map.find(mass) == peak_nodes_map.end())
            {
                std::unique_ptr<PeakNode> node = 
                    std::make_unique<PeakNode>();
                // set mass
                node->set_mass(mass);
                // set matches
                if (complex_)
                    node->Add(it.first, kY1, std::vector<int>());
                if (hybrid_)
                    node->Add(it.first, kY1_mannose, std::vector<int>());
                if (highmannose_)
                    node->Add(it.first, kY1_hybrid, std::vector<int>());
                // add node 
                peak_nodes_map[mass] = std::move(node);
                // enqueue
                queue.push(peak_nodes_map[mass].get());
            }
            else
            {
                // update glycopeptide match
                if (complex_)
                    peak_nodes_map[mass]->Add(it.first, kY1, std::vector<int>());
                if (hybrid_)
                    peak_nodes_map[mass]->Add(it.first, kY1_mannose, std::vector<int>());
                if (highmannose_)
                    peak_nodes_map[mass]->Add(it.first, kY1_hybrid, std::vector<int>());
            }
        }
    }


    void InitSearch(const std::vector<model::spectrum::Peak>& peaks, int max_charge)
    {
        std::vector<std::shared_ptr<algorithm::search::Point<int>>> peak_points; 
        for(int i = 0; i < (int) peaks.size(); i++)
        {
            const model::spectrum::Peak peak = peaks[i];
            for(int charge = 1; charge <= max_charge; charge++)
            {
                double mass = util::mass::SpectrumMass::Compute(peak.MZ(), charge);
                std::shared_ptr<algorithm::search::Point<int>> point 
                    = std::make_shared<algorithm::search::Point<int>>(mass, i);
                peak_points.push_back(std::move(point));
            }
        }
        searcher_->Init(peak_points);
    }

    std::unique_ptr<algorithm::search::ISearch<int>> searcher_;
    const std::unordered_map<std::string, std::unique_ptr<model::glycan::Glycan>>& glycans_map_;
    bool complex_;
    bool hybrid_;
    bool highmannose_;
    const std::string kY1 = "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ";
    const std::string kY1_hybrid = "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ";
    const std::string kY1_mannose = "1 0 0 0 0 0 ";
    const int kMissing = 4;
    std::unordered_map<std::string, double> peptide_mass_;
};



} // namespace engine
} // namespace search

#endif