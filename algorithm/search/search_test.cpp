#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <chrono>
#include <iostream>
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
    BucketSearch<double> searcher(model::spectrum::ToleranceBy::Dalton, 20, 1, 1000);
    std::vector<std::shared_ptr<Point<double>>> box; 

    for(int i=2; i<100; i++)
    {
        box.push_back(CreatePoint(i));
    }

    searcher.Init(box);
    std::cout << "Init Done!" << std::endl;

    std::vector<double> res = searcher.Search(40);
    BOOST_CHECK(searcher.Match(50));
    BOOST_CHECK(res.size() == 39);
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

    BucketSearch<model::spectrum::Peak> searcher(model::spectrum::ToleranceBy::PPM, 200, 100, 2000);
    searcher.Init(mz_points);


    auto start = std::chrono::high_resolution_clock::now(); 
    std::vector<model::spectrum::Peak> res = searcher.Search(1089.1);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    for(auto& it : res)
    {
        std::cout << it.MZ() << " " 
            << util::mass::SpectrumMass::ComputePPM(1089.1, it.MZ()) << std::endl;
    }

    BOOST_CHECK(searcher.Match(113.0661));

}


} // namespace algorithm
} // namespace search 