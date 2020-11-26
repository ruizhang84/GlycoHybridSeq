#define BOOST_TEST_MODULE BinPackingTest
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "binpacking.h"

namespace engine {
namespace spectrum {

BOOST_AUTO_TEST_CASE( read_training_data ) 
{
    std::vector<model::spectrum::Peak> peaks;
    peaks.push_back(model::spectrum::Peak(1, 1));
    peaks.push_back(model::spectrum::Peak(1, 2));
    peaks.push_back(model::spectrum::Peak(2, 3));
    peaks.push_back(model::spectrum::Peak(3, 1));
    peaks.push_back(model::spectrum::Peak(4, 4));
    peaks.push_back(model::spectrum::Peak(5, 1));
    peaks.push_back(model::spectrum::Peak(6, 1));
    peaks.push_back(model::spectrum::Peak(7, 8));
    peaks.push_back(model::spectrum::Peak(8, 1));
    peaks.push_back(model::spectrum::Peak(9, 1));
    peaks.push_back(model::spectrum::Peak(10, 12));
    peaks.push_back(model::spectrum::Peak(11, 1));

    model::spectrum::Spectrum spec;
    spec.set_peaks(peaks);
    BinPacking packer(model::spectrum::ToleranceBy::PPM, 2, 0.000001, 20);
    std::vector<double> bins = packer.Packing(spec);
    for(int i = 0; i < (int) bins.size(); i++)
    {
        if (bins[i] > 0)
        {
            std::cout << i << ": " << bins[i] << std::endl;
        }
    }

}


} // namespace spectrum
} // namespace engine
