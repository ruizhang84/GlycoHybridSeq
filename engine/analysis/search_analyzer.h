#ifndef ENGINE_ANALYSIS_SEARCH_ANALYZER_H_
#define ENGINE_ANALYSIS_SEARCH_ANALYZER_H_


#include "../../model/spectrum/peak.h"
#include "../../model/glycan/glycan.h"
#include "search_result.h"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath> 
#include <memory>
#include <algorithm>



namespace engine{
namespace analysis{

class SearchAnalyzer
{
public:
    SearchAnalyzer() = default;

    std::vector<SearchResult> Analyze(
        int scan,
        const std::vector<model::spectrum::Peak>& peaks,
        const std::unordered_map<std::string, std::unordered_set<int>>& peptide_results,
        const std::unordered_map<std::string, std::unordered_set<int>>& glycan_results)
    {
        std::vector<SearchResult> results;
        // preprocess to extract peptide
        std::unordered_map<std::string, std::vector<std::string>> peptides_map;
        std::unordered_map<std::string, std::vector<std::string>> glycans_map;
        for (const auto& it : peptide_results)
        {
            std::string peptide = it.first.substr(0, it.first.find("|"));
            
            if (peptides_map.find(peptide) == peptides_map.end())
            {
                peptides_map[peptide] = std::vector<std::string>();
            }
            peptides_map[peptide].push_back(it.first);
        }
        for (const auto& it : glycan_results)
        {
            std::string peptide = it.first.substr(0, it.first.find("|"));
            if (glycans_map.find(peptide) == glycans_map.end())
            {
                glycans_map[peptide] = std::vector<std::string>();
            }
            glycans_map[peptide].push_back(it.first);
        }
        
        // analyze the best matches
        double best_score = 0;
        for(const auto& it : peptides_map)
        {
            std::string peptide = it.first;
            for(const auto& p : it.second)
            {
                
                for(const auto& g : glycans_map[peptide])
                {
                    // get index
                    std::unordered_set<int> peptides_index = peptide_results.find(p)->second;
                    std::unordered_set<int> glycans_index = glycan_results.find(g)->second;
                    // compute score
                    double score = ComputePeakScore(peaks, peptides_index, glycans_index);
                    // create results if higher score
                    if (score > best_score)
                    {
                        best_score = score;
                        results.clear();
                    }
                    if (score == best_score)
                    {
                        int pos = std::stoi(p.substr(peptide.length()+1));
                        std::string glycan = g.substr(peptide.length()+1);
                        SearchResult r;
                        r.set_glycan(glycan);
                        r.set_peptide(peptide);
                        r.set_scan(scan);
                        r.set_site(pos);
                        r.set_score(score);
                        results.push_back(r);
                    }
                }
            }
        }
 
        return results;
    }

    double ComputePeakScore(const std::vector<model::spectrum::Peak>& peaks, 
        const std::unordered_set<int>& peptides_index, 
        const std::unordered_set<int>& glycans_index) const
    {
        double sum = 0;
        std::vector<int> peak_index;
        peak_index.insert(peak_index.end(), peptides_index.begin(), peptides_index.end());
        peak_index.insert(peak_index.end(), glycans_index.begin(), glycans_index.end());

        double peptide_score = 0;
        double glycan_score = 0;
        for(int index : peptides_index)
        {
            peptide_score += log(peaks[index].Intensity()) ;
        }
        for(int index : glycans_index)
        {
            glycan_score += log(peaks[index].Intensity()) ;
        }

        for(const auto& it : peaks)
        {
            sum += log(it.Intensity());
        }

        return sqrt(peptide_score * glycan_score) / sum;
    }

};



} // namespace analysis
} // namespace search

#endif