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
#include "precursor_match.h"

#include <chrono> 

namespace engine{
namespace search {

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
    int special_scan = 6879;
    double ms1_tol = 10;
    model::spectrum::ToleranceBy ms1_by = model::spectrum::ToleranceBy::PPM;
    std::unique_ptr<algorithm::search::ISearch<std::string>> searcher =
        std::make_unique<algorithm::search::BucketSearch<std::string>>(ms1_by, ms1_tol);

    PrecursorMatcher precursor_runner(std::move(searcher));
    precursor_runner.Init(peptides, builder->GlycanMapsRef());
    auto special_spec = spectrum_reader.GetSpectrum(special_scan);

    auto start = std::chrono::high_resolution_clock::now(); 

    auto results = precursor_runner.Match(special_spec.PrecursorMZ(), special_spec.PrecursorCharge());   

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    std::cout << special_spec.Scan() << " : " << std::endl;
    for(auto it : results)
    {
        std::cout << it.first << std::endl;
        for(auto g: it.second)
        {
            std::cout << g->Name() << std::endl;
        }
    }


   

}

} // namespace search
} // namespace engine