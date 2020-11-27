#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include "bin_search.h"
#include "bucket_search.h"
#include "../../util/io/mgf_parser.h"
#include <unordered_map>


using namespace std;
namespace algorithm {
namespace search {

std::shared_ptr<Point<double>> CreatePoint(double value)
{
    return make_shared<Point<double>>(value, value);
}

BOOST_AUTO_TEST_CASE( Algorithm_test ) 
{
    BinSearch searcher(model::spectrum::ToleranceBy::Dalton, 20);
    std::vector<double> box; 

    for(int i=2; i<100; i++)
    {
        box.push_back(i);
    }

    searcher.Init(box);
    std::cout << "Init Done!" << std::endl;

    BOOST_CHECK(searcher.Match(50));
}


BOOST_AUTO_TEST_CASE( spectrum_search_test ) 
{

    std::unique_ptr<util::io::SpectrumParser> parser = std::make_unique<util::io::MGFParser>();
    util::io::SpectrumReader spectrum_reader("/home/ruiz/Documents/GlycoCrushSeq/data/ZC_20171218_C22_R1.mgf", std::move(parser));

    model::spectrum::Spectrum spectrum = spectrum_reader.GetSpectrum(64);
    std::vector<std::shared_ptr<Point<model::spectrum::Peak>>> mz_points; 

    for(const auto& it : spectrum.Peaks())
    {
        std::shared_ptr<Point<model::spectrum::Peak>> p = 
            std::make_shared<Point<model::spectrum::Peak>>(it.MZ(), it);
        mz_points.push_back(std::move(p));
    }


    auto start = std::chrono::high_resolution_clock::now(); 
    BucketSearch<model::spectrum::Peak> searcher(model::spectrum::ToleranceBy::PPM, 200);
    searcher.Init(mz_points);
   
    std::vector<model::spectrum::Peak> res = searcher.Search(662.0826);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    for(auto& it : res)
    {
        std::cout << it.MZ() << " " 
            << util::mass::SpectrumMass::ComputePPM(662.0826, it.MZ()) << std::endl;
    }

    BOOST_CHECK(searcher.Match(113.0671));

}


} // namespace algorithm
} // namespace search 