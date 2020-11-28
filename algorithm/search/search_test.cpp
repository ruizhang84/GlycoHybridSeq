#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include "bucket_search.h"
#include "binary_search.h"
#include "search.h"
#include "../../util/io/mgf_parser.h"
#include "../../util/mass/spectrum.h"
#include <unordered_map>


using namespace std;
namespace algorithm {
namespace search {

std::shared_ptr<Point<double>> CreatePoint(double value)
{
    return make_shared<Point<double>>(value, value);
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

    std::shared_ptr<Point<model::spectrum::Peak>> p = 
            std::make_shared<Point<model::spectrum::Peak>>(662.0826, model::spectrum::Peak(662.0826, 2));
    searcher.Add(p);
   
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

BOOST_AUTO_TEST_CASE( spectrum_search_test2 ) 
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
    BinarySearch<model::spectrum::Peak> searcher(model::spectrum::ToleranceBy::PPM, 200);
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

BOOST_AUTO_TEST_CASE( spectrum_search_test3 ) 
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

    std::unique_ptr<BinarySearch<model::spectrum::Peak>> searcher = 
        std::make_unique< BinarySearch<model::spectrum::Peak>>(model::spectrum::ToleranceBy::PPM, 200);

    searcher->Init(mz_points);

    std::vector<model::spectrum::Peak> res = searcher->Search(662.0826);
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    for(auto& it : res)
    {
        std::cout << it.MZ() << " " 
            << util::mass::SpectrumMass::ComputePPM(662.0826, it.MZ()) << std::endl;
    }

    BOOST_CHECK(searcher->Match(113.0671));

}


} // namespace algorithm
} // namespace search 