#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <iostream>
#include "bucket_search.h"
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


} // namespace algorithm
} // namespace search 