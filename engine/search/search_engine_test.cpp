#define BOOST_TEST_MODULE SearchTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <iomanip>

#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../protein/protein_digest.h"
#include "../protein/protein_ptm.h"
#include "../glycan/glycan_builder.h"
#include "../../algorithm/search/bucket_search.h"
#include "../../algorithm/search/binary_search.h"
#include "search_sequence.h"
#include "search_glycan.h"
#include "search_helper.h"
#include "precursor_match.h"
#include "../analysis/search_result.h"

#include <chrono>

namespace engine{
namespace search {

BOOST_AUTO_TEST_CASE( search_engine_test )
{
    // read spectrum
    std::string path = "/home/ruiz/Documents/GlycoCrushSeq/data/ZC_20171218_H95_R1.mgf";
    std::unique_ptr<util::io::SpectrumParser> parser = std::make_unique<util::io::MGFParser>();
    util::io::SpectrumReader spectrum_reader(path, std::move(parser));

    int start_scan = spectrum_reader.GetFirstScan();
    int last_scan = spectrum_reader.GetLastScan();
    BOOST_CHECK(start_scan < last_scan);

    // process spectrum by normalization
    model::spectrum::Spectrum spec = spectrum_reader.GetSpectrum(start_scan);

    // read fasta and build peptides
    util::io::FASTAReader fasta_reader("/home/ruiz/Documents/GlycoCrushSeq/data/haptoglobin.fasta");
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();

    // for(auto& p: proteins)
    // {
    //     std::string seq = p.Sequence();
    //     std::reverse(seq.begin(), seq.end());
    //     p.set_sequence(seq);
    // }

    engine::protein::Digestion digest;
    digest.SetProtease(engine::protein::Proteases::Trypsin);
    std::unordered_set<std::string> seqs = digest.Sequences(proteins.front().Sequence(),
         engine::protein::ProteinPTM::ContainsNGlycanSite);
    digest.SetProtease(engine::protein::Proteases::GluC);
    std::vector<std::string> peptides;
    peptides.insert(peptides.end(), seqs.begin(), seqs.end());

    for(auto& it : seqs)
    {
        std::unordered_set<std::string> seq = digest.Sequences(it,
         engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(peptides.end(), seq.begin(), seq.end());
    }

    // std::unordered_set<std::string> seen;
    // int count = 0;
    // for(auto it : peptides)
    // {
    //     if (seen.find(it) == seen.end())
    //     {
    //         seen.insert(it);
    //         count++;
    //     }
            
    // }
    // std::cout << count << std::endl;
    // BOOST_CHECK(std::find(peptides.begin(), peptides.end(), "VVLHPNYSQVD") != peptides.end());
    // BOOST_CHECK(std::find(peptides.begin(), peptides.end(), "KDNLTYVGDGETR") != peptides.end());

    // // build glycans
    int hexNAc = 12, hex = 12, Fuc = 5, NeuAc = 4, NeuGc = 0;
    std::unique_ptr<engine::glycan::GlycanBuilder> builder =
        std::make_unique<engine::glycan::GlycanBuilder>(hexNAc, hex, Fuc, NeuAc, NeuGc);
    builder->Build();
    // std::cout << builder->GlycanMapsRef().size() << std::endl;


    // spectrum matching
    int special_scan = 8019;
    double ms1_tol = 10;
    model::spectrum::ToleranceBy ms1_by = model::spectrum::ToleranceBy::PPM;
    std::unique_ptr<algorithm::search::ISearch<std::string>> searcher =
        std::make_unique<algorithm::search::BucketSearch<std::string>>(ms1_by, ms1_tol);

    PrecursorMatcher precursor_runner(std::move(searcher));
    precursor_runner.Init(peptides, builder->GlycanMapsRef());
    auto special_spec = spectrum_reader.GetSpectrum(special_scan);

    auto results = precursor_runner.Match(special_spec.PrecursorMZ(), special_spec.PrecursorCharge());
    std::cout << "scan " << special_spec.Scan() << " : " << std::endl;
    for(auto it : results)
    {
        std::cout << it.first << std::endl;
        for(auto g: it.second)
        {
            std::cout << g->Name() << std::endl;
        }
    }

    // search peptide
    auto start = std::chrono::high_resolution_clock::now();
    double ms2_tol = 0.01;
    model::spectrum::ToleranceBy ms2_by = model::spectrum::ToleranceBy::Dalton;
    std::unique_ptr<algorithm::search::ISearch<std::string>> more_searcher =
        std::make_unique<algorithm::search::BucketSearch<std::string>>(ms2_by, ms2_tol);

    SequenceSearch spectrum_runner(std::move(more_searcher));
    auto peptide_results = spectrum_runner.Search(special_spec.Peaks(), special_spec.PrecursorCharge(), results);

    // search glycan
    std::unique_ptr<algorithm::search::ISearch<int>> extra_searcher =
        std::make_unique<algorithm::search::BucketSearch<int>>(ms2_by, ms2_tol);
    GlycanSearch spectrum_searcher(std::move(extra_searcher), builder->GlycanMapsRef());
    auto glycan_results = spectrum_searcher.Search(special_spec.Peaks(), special_spec.PrecursorCharge(), results);

    

    for(const auto& it : glycan_results)
    {
        std::cout << it.first << " :" << it.second.size() << " " << SearchHelper::ComputePeakScore(special_spec.Peaks(), it.second) << std::endl;

    }
    for(const auto& it : peptide_results)
    {
        std::cout << it.first << " :"  << it.second.size() << " " << SearchHelper::ComputePeakScore(special_spec.Peaks(), it.second) << std::endl;

    }

    // engine::analysis::SearchAnalyzer analyzer;
    // auto searched = analyzer.Analyze(special_scan, special_spec.Peaks(), peptide_results, glycan_results);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << duration.count() << std::endl;
    // for(const auto& r : searched)
    // {
    //     std::cout << r.Sequence() << " | " << r.Glycan() << " | " << r.Scan() << " " << r.Score() << std::endl;
    // }

   
}



} // namespace search
} // namespace engine