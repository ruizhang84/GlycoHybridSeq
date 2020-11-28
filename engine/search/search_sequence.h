#ifndef ENGINE_SEARCH_SEQUENCE_H_
#define ENGINE_SEARCH_SEQUENCE_H_

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


class SequenceSearch
{
public:
    SequenceSearch(std::unique_ptr<algorithm::search::ISearch<std::string>> searcher):
        searcher_(std::move(searcher)){}

    // peptide seq, glycan*
    std::unordered_map<std::string, std::vector<model::spectrum::Peak>> Search(
        const std::vector<model::spectrum::Peak>& peaks, int max_charge, 
        const std::unordered_map<std::string, std::vector<model::glycan::Glycan*>>& candidate)
    {
        InitSearch(candidate);
        std::unordered_map<std::string, std::unordered_set<int>> results;

        // search peaks
        for(int i = 0; i < (int) peaks.size(); i++)
        {
            const model::spectrum::Peak peak = peaks[i];
            for(int charge = 1; charge < max_charge; charge++)
            {
                double target = util::mass::SpectrumMass::Compute(peak.MZ(), charge);
                std::vector<std::string> glyco_sequences = searcher_->Search(target);
                for(std::string seq: glyco_sequences)
                {
                    if(results.find(seq) == results.end())
                    {
                        results[seq] = std::unordered_set<int>();
                    }
                    results[seq].insert(i);
                }
            }
        }

        // convert to peak results
        std::unordered_map<std::string, std::vector<model::spectrum::Peak>> peak_results;
        for(const auto& it : results)
        {
            peak_results[it.first] = std::vector<model::spectrum::Peak>();
            for(int index : it.second)
            {
                peak_results[it.first].push_back(peaks[index]);
            }
        }
        return peak_results;
    }
    

protected:
    double ComputeGlycanMass(const std::vector<model::glycan::Glycan*>& glycans)
    {
        double sums = 0.0;
        int glycan_count = (int) glycans.size();
        for(const auto& glycan : glycans)
        {
            sums += glycan->Mass();
        }
        return sums * 1.0 / glycan_count;
    }

    void InitSearch(
        const std::unordered_map<std::string, std::vector<model::glycan::Glycan*>>& candidate)
    {
        std::vector<std::shared_ptr<algorithm::search::Point<std::string>>> peptides_points;
        for(const auto& it : candidate)
        {
            std::string peptide = it.first;
            // get glycan mass
            double glycan_mean_mass = ComputeGlycanMass(it.second);

            // create points
            for (int pos : engine::protein::ProteinPTM::FindNGlycanSite(peptide))
            {
                std::string table_key = SearchHelper::MakeKeySequence(peptide, pos);
                //init table
                if (mass_table_.find(table_key) == mass_table_.end())
                {
                    
                    std::vector<double> mass_list = SearchHelper::ComputeNonePTMPeptideMass(peptide, pos);
                    mass_table_[table_key] = std::vector<double>(mass_list);
                    mass_list = SearchHelper::ComputePTMPeptideMass(peptide, pos);
                    ptm_mass_table_[table_key] = std::vector<double>(mass_list);
                }
                // retreive table
                for(double mass : mass_table_[table_key])
                {
                    std::shared_ptr<algorithm::search::Point<std::string>> point = 
                        std::make_shared<algorithm::search::Point<std::string>>(mass, table_key);
                    peptides_points.push_back(std::move(point));
                }
                for(double mass : ptm_mass_table_[table_key])
                {
                    std::shared_ptr<algorithm::search::Point<std::string>> point = 
                        std::make_shared<algorithm::search::Point<std::string>>(mass + glycan_mean_mass, table_key);
                    peptides_points.push_back(std::move(point));
                }
            }
        } 

        searcher_->Init(peptides_points);
    } 


    std::unique_ptr<algorithm::search::ISearch<std::string>> searcher_;
    // mass table to store the computed ptm / non-ptm results
    std::unordered_map<std::string, std::vector<double>> ptm_mass_table_;
    std::unordered_map<std::string, std::vector<double>> mass_table_;
};



} // namespace engine
} // namespace search

#endif