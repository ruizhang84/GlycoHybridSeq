#define BOOST_TEST_MODULE BinPackingTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "../../model/spectrum/spectrum.h"
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
    BinPacking packer(model::spectrum::ToleranceBy::Dalton, 2, 0.1, 20);
    std::unordered_map<int, model::spectrum::Peak> bins = packer.Packing(spec);
    for(const auto& it : bins)
    {
        std::cout << it.first << ": " << it.second.MZ() << ", " << it.second.Intensity() << std::endl;

    }

}


} // namespace spectrum
} // namespace engine
