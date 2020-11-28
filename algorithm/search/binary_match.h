#ifndef ALGORITHM_BINARY_MATCH_H_
#define ALGORITHM_BINARY_MATCH_H_

#include <vector>
#include <cstdlib> 
#include <algorithm>
#include "search.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/spectrum.h"

namespace algorithm {
namespace search {

// enhanced binary search
class BinaryMatch
{
public:
    BinaryMatch(model::spectrum::ToleranceBy type, double tol):
        type_(type), tolerance_(tol) {};

    virtual void Init(std::vector<double> inputs, bool sorted=false) 
    {
        if (!inputs.empty() && !sorted)
            std::sort(inputs.begin(), inputs.end());
        data_ = inputs;
    }

    double Tolerance() const { return tolerance_; }
    model::spectrum::ToleranceBy ToleranceType() const { return type_; }
    void set_tolerance(double tol) { tolerance_ = tol; }
    void set_tolerance_by(model::spectrum::ToleranceBy type) { type_ = type; }

    virtual bool ISMatch(const double p, const double target)
    {
        switch (type_)
        {
        case model::spectrum::ToleranceBy::PPM:
            return util::mass::SpectrumMass::ComputePPM(p, target) < tolerance_;
        case model::spectrum::ToleranceBy::Dalton:
            return std::abs(p - target) < tolerance_;
        default:
            break;
        }
        return false;
    }


    virtual bool Match(const double target)
    {
        if (data_.empty()) 
            return false;

        int start = 0, end = data_.size()-1;
        while (start <= end)
        {
            int mid = (end - start) / 2 + start;
            if (ISMatch(data_[mid], target))
                return true;
            else if (data_[mid] < target)
                start = mid + 1;
            else
                end = mid - 1;
        }
        return false;
    }

protected:
    double tolerance_; 
    model::spectrum::ToleranceBy type_;
    std::vector<double> data_;
    double lower_;
    double upper_;
};

} // namespace algorithm
} // namespace search 

#endif
