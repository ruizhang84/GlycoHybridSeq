#define BOOST_TEST_MODULE SearchTest
#include <boost/test/unit_test.hpp>
#include "multi_comparison.h"


namespace engine{
namespace analysis{

BOOST_AUTO_TEST_CASE( search_engine_test ) 
{
    std::vector<double> v {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    MultiComparison tester(0.01);
    // BOOST_CHECK(tester.mean(v) == 5);
    // BOOST_CHECK_CLOSE(tester.stdv(v), 2.738613, 0.001);
    // BOOST_CHECK_CLOSE(tester.pValue(v, 2.5), 0.8193448, 0.001);

}


} // namespace analysis
} // namespace engine