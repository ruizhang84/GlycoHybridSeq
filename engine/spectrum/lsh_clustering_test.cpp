#define BOOST_TEST_MODULE BinPackingTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "../../model/spectrum/spectrum.h"
#include "lsh_clustering.h"
#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "spectrum_sim.h"
#include <chrono>


namespace engine {
namespace spectrum {

BOOST_AUTO_TEST_CASE( cluster_data ) 
{
    // read spectrum
    std::string path = "/home/ruiz/Documents/GlycoCrushSeq/data/MGF2/ZC_20171218_C16_R1.mgf";
    std::unique_ptr<util::io::SpectrumParser> parser = std::make_unique<util::io::MGFParser>();
    util::io::SpectrumReader spectrum_reader(path, std::move(parser));

    int start_scan = spectrum_reader.GetFirstScan();
    int last_scan = spectrum_reader.GetLastScan();
    BOOST_CHECK(start_scan < last_scan);

    // process spectrum by normalization
    model::spectrum::Spectrum spec = spectrum_reader.GetSpectrum(start_scan);

    // clustering
    auto start = std::chrono::high_resolution_clock::now();

    LSHClustering cluster(model::spectrum::ToleranceBy::Dalton, 0.01);
    std::vector<EmbededSpectrum> embeds = cluster.Embed(spectrum_reader.GetSpectrum());
    std::unordered_map<int, std::vector<EmbededSpectrum>> hash_table =  cluster.Hashing(embeds);

    SpectrumSim sim;

    // for(const auto& it : hash_table)
    // {
    //     if (it.second.size() > 10)
    //     {
           
    //         for(const auto& s : it.second)
    //         {
    //             for(const auto& ss : it.second)
    //             {
    //                 double cos = sim.ComputeCosine(s.Embed(), ss.Embed());
    //                 std::cout << cos << std::endl;
    //             }
    //         }
    //     }
    // }


    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << duration.count() << std::endl;

    int count = 0;
    for(auto& it : hash_table)
    {
        if (it.second.size() > 2)
        {
            std::cout << it.first << " : " << it.second.front().Scan() << ", " << it.second.back().Scan() << " , " 
                << sim.ComputeCosine(it.second.front().Embed(), it.second.back().Embed())  << std::endl;
            count += (int) it.second.size();
        }
            
    }
    std::cout << count << std::endl;
}


} // namespace spectrum
} // namespace engine
