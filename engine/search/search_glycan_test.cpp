#define BOOST_TEST_MODULE SearchTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <iomanip>
#include <queue> 
#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../protein/protein_digest.h"
#include "../protein/protein_ptm.h"
#include "../glycan/glycan_builder.h"
#include "../../algorithm/search/bucket_search.h"
#include "../../algorithm/search/binary_search.h"
#include "search_helper.h"
#include "precursor_match.h"
#include "search_glycan_helper.h"
// #include "search_glycan.h"

#include <chrono> 

namespace engine{
namespace search {

bool Satisify(const std::unordered_map<std::string, std::unique_ptr<model::glycan::Glycan>>& glycans_map_,
    const std::string& identified_glycan_id, const model::glycan::Glycan* glycan)
{
    const std::vector<int> identified_glycan_table = 
        glycans_map_.find(identified_glycan_id)->second->TableConst();
    const std::vector<int> candidate_glycan_table = glycan->TableConst();
    for(int i = 0; i < (int) identified_glycan_table.size(); i++)
    {
        if (candidate_glycan_table[i] < identified_glycan_table[i])
            return false;
    }
    return true;
}


BOOST_AUTO_TEST_CASE( search_engine_test ) 
{
    // read spectrum
    std::string path = "/home/ruiz/Documents/GlycoCrushSeq/data/ZC_20171218_C22_R1.mgf";
    std::unique_ptr<util::io::SpectrumParser> parser = std::make_unique<util::io::MGFParser>();
    util::io::SpectrumReader spectrum_reader(path, std::move(parser));


    int start_scan = spectrum_reader.GetFirstScan();
    int last_scan = spectrum_reader.GetLastScan();
    BOOST_CHECK(start_scan < last_scan);

    // process spectrum by normalization
    model::spectrum::Spectrum spec = spectrum_reader.GetSpectrum(start_scan);

    // read fasta and build peptides
    // util::io::FASTAReader fasta_reader("/home/ruiz/Documents/Glycoseq-Cpp/data/haptoglobin.fasta");
     util::io::FASTAReader fasta_reader("/home/ruiz/Documents/Glycoseq-Cpp/data/titin.fasta");
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
 
    engine::protein::Digestion digest;
    digest.SetProtease(engine::protein::Proteases::Trypsin);
    std::unordered_set<std::string> seqs = digest.Sequences(proteins.front().Sequence(),
         engine::protein::ProteinPTM::ContainsNGlycanSite);
    digest.SetProtease(engine::protein::Proteases::GluC);
    std::vector<std::string> peptides;
    for(auto& it : seqs)
    {
        std::unordered_set<std::string> seq = digest.Sequences(it,
         engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(peptides.end(), seq.begin(), seq.end());
    }
    // BOOST_CHECK(std::find(peptides.begin(), peptides.end(), "VVLHPNYSQVD") != peptides.end());


    // // build glycans
    int hexNAc = 12, hex = 12, Fuc = 5, NeuAc = 4, NeuGc = 0;
    std::unique_ptr<engine::glycan::GlycanBuilder> builder =
        std::make_unique<engine::glycan::GlycanBuilder>(hexNAc, hex, Fuc, NeuAc, NeuGc);
    builder->Build();


    // spectrum matching
    int special_scan = 6765;
    double ms1_tol = 10;
    model::spectrum::ToleranceBy ms1_by = model::spectrum::ToleranceBy::PPM;
    std::unique_ptr<algorithm::search::ISearch<std::string>> searcher =
        std::make_unique<algorithm::search::BucketSearch<std::string>>(ms1_by, ms1_tol);

    PrecursorMatcher precursor_runner(std::move(searcher));
    precursor_runner.Init(peptides, builder->GlycanMapsRef());
    auto special_spec = spectrum_reader.GetSpectrum(special_scan);

    auto results = precursor_runner.Match(special_spec.PrecursorMZ(), special_spec.PrecursorCharge());   
    // std::cout << special_spec.Scan() << " : " << std::endl;
    // for(auto it : results)
    // {
    //     std::cout << it.first << std::endl;
    //     for(auto g: it.second)
    //     {
    //         std::cout << g->Name() << std::endl;
    //     }
    // }

    // search glycans
    auto start = std::chrono::high_resolution_clock::now(); 
    double ms2_tol = 0.01;
    model::spectrum::ToleranceBy ms2_by = model::spectrum::ToleranceBy::Dalton; 
    std::unique_ptr<algorithm::search::ISearch<int>> more_searcher =
        std::make_unique<algorithm::search::BucketSearch<int>>(ms2_by, ms2_tol);

    // GlycanSearch spectrum_runner(std::move(more_searcher), builder->GlycanMapsRef());
    // auto peptide_results = spectrum_runner.Search(special_spec.Peaks(), special_spec.PrecursorCharge(), results);
    
    const std::unordered_map<std::string, std::unique_ptr<model::glycan::Glycan>>& glycans_map_ = builder->GlycanMapsRef();
    const std::string kY1 = "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ";
    std::unordered_map<std::string, double> peptide_mass_;

    // init search
    std::vector< model::spectrum::Peak> peaks = special_spec.Peaks();
    std::vector<std::shared_ptr<algorithm::search::Point<int>>> peak_points; 
    for(int i = 0; i < (int) peaks.size(); i++)
    {
        const model::spectrum::Peak peak = peaks[i];
        for(int charge = 1; charge <= special_spec.PrecursorCharge(); charge++)
        {
            double mass = util::mass::SpectrumMass::Compute(peak.MZ(), charge);
            std::shared_ptr<algorithm::search::Point<int>> point 
                = std::make_shared<algorithm::search::Point<int>>(mass, i);
            peak_points.push_back(std::move(point));
        }
    }
    more_searcher->Init(peak_points);
    // std::cout << peak_points.size() << " " << special_spec.PrecursorCharge() << " " << peaks.size()  << std::endl;

    // init priority queue
    std::unordered_map<double, std::unique_ptr<PeakNode>> peak_nodes_map;
    std::priority_queue<PeakNode*, std::vector<PeakNode*>, PeakNodeComparison> queue;
    std::vector<PeakNode*> matched_nodes;

    for(const auto& it : results)
    {
        // Y1 mass
        double mass = util::mass::PeptideMass::Compute(it.first) + util::mass::GlycanMass::kHexNAc;
        // node matches
        if (peak_nodes_map.find(mass) == peak_nodes_map.end())
        {
            std::unique_ptr<PeakNode> node = 
                std::make_unique<PeakNode>();
            // set mass
            node->set_mass(mass);
            // set matches
            node->Add(it.first, kY1, std::vector<int>());
            // add node 
            peak_nodes_map[mass] = std::move(node);
            // enqueue
            queue.push(peak_nodes_map[mass].get());
        }
        else
        {
            // update glycopeptide match
            peak_nodes_map[mass]->Add(it.first, kY1, std::vector<int>());
        }
    }

    // while (queue.size() > 0)
    // {
    //     PeakNode node = queue.top();
    //     queue.pop();
    //     std::cout << node.Mass() << " :" << std::endl;
    //     for(const auto& it : node.Matches())
    //     {
    //         std::cout << it.first << std::endl;
    //         for(const auto& item : it.second)
    //         {
    //             std::cout << item.first << std::endl;
    //         }
    //     }
    // }

    // for(auto& it : peak_nodes_map)
    // {
    //     std::cout << it.first<< " : " << std::endl;
    //     for(const auto& j : it.second.Matches())
    //     {
    //         std::cout << j.first << std::endl;
    //     }
    // }
    
    // dynamic programming

        while (queue.size() > 0)
        {
            // get node
            PeakNode* node = queue.top();
            queue.pop();

            // // match peaks
            double target = node->Mass();
            std::vector<int> matched = more_searcher->Search(target);
            // std::cout << matched.size() << std::endl;

            // max if matched a peak
            node->Max(peaks);
            
            // update matches
            if (matched.size() > 0)
            {
                // for(const auto& it : node.Matches())
                // {
                //     std::cout << it.first << std::endl;
                //     for(const auto& item : it.second)
                //     {
                //         std::cout << item.first << 
                //             SearchHelper::ComputePeakScore(peaks, item.second) << std::endl;
                //     }
                // }
                // std::cout << std::endl;
                

                node->Add(matched);
                node->set_miss(0);

                // for(const auto& it : node.Matches())
                // {
                //     std::cout << it.first << std::endl;
                //     for(const auto& item : it.second)
                //     {
                //         std::cout << item.first << 
                //             SearchHelper::ComputePeakScore(peaks, item.second) << std::endl;
                //     }
                // }
                matched_nodes.push_back(node);

            }

            if (node->Missing() > 5)
                continue;

            // extending queue
            for(const auto& it : node->Matches())
            {
                std::string peptide = it.first;
                for(const auto& gt : it.second)
                {
                    std::string glycan_id = gt.first;
                    std::vector<int> peak_indexes(gt.second.begin(), gt.second.end());

                    model::glycan::Glycan* glycan = glycans_map_.find(glycan_id)->second.get();
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

            // if (queue.size() > 8)
            // {
            //     while (queue.size() > 0)
            //     {
            //         PeakNode node = queue.top();
            //         queue.pop();
            //         std::cout << node.Mass() << " :" << std::endl;
            //         for(const auto& it : node.Matches())
            //         {
            //             std::cout << it.first << std::endl;
            //             for(const auto& item : it.second)
            //             {
            //                 std::cout << item.first << " " << SearchHelper::ComputePeakScore(peaks, item.second) << std::endl;
            //             }
            //         }
            //     }
            // }
        }
        // std::cout << peak_nodes_map.size() << std::endl;
        // for(auto& it : peak_nodes_map)
        // {
        //     std::cout << it.first<< " : " << std::endl;
        //     for(const auto& j : it.second->Matches())
        //     {
        //         std::cout << j.first << std::endl;
        //         for(const auto& item : j.second)
        //         {
        //             std::cout << glycans_map_.find(item.first)->second->Name() <<  item.first << " " << SearchHelper::ComputePeakScore(peaks, item.second) << std::endl;
        //         }
        //     }
        //     std::cout << std::endl;
        // }

        // for(const auto& it : matched_nodes)
        // {
        //     for(const auto& j : it->Matches())
        //     {
        //         std::cout << j.first << std::endl;
        //         for(const auto& item : j.second)
        //         {
        //             std::cout << glycans_map_.find(item.first)->second->Name() <<  item.first << " " << SearchHelper::ComputePeakScore(peaks, item.second) << std::endl;
        //         }
        //     }
        //     std::cout << std::endl;
        // }

        std::unordered_map<std::string, std::unordered_set<int>> more_results;
        for(const auto& node : matched_nodes)
        {
            for(const auto& it : node->Matches())
            {
                // std::string peptide_seq = it.first;
                // std::unordered_map<std::string, std::unordered_set<int>> matched = it.second;
                std::vector<model::glycan::Glycan*> candidate_glycan = results.find(it.first)->second;
                for(const auto& g : it.second)
                {
                    for(const auto& glycan : candidate_glycan)
                    {
                        if (Satisify(glycans_map_, g.first, glycan))
                        {
                            std::string key = SearchHelper::MakeKeyGlycoSequence(glycan->ID(), it.first);
                            if(more_results.find(key) == more_results.end())
                            {
                                more_results[key] = std::unordered_set<int>();
                            }
                            more_results[key].insert(g.second.begin(), g.second.end());
                        }
                    }

                }
            }
        }

        for(const auto& it : more_results)
        {

            std::cout << it.first << " :"  << SearchHelper::ComputePeakScore(peaks, it.second) << std::endl;
            
        }



    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    
}



} // namespace search
} // namespace engine