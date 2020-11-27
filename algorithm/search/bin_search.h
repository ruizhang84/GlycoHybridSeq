#ifndef ALGORITHM_BIN_SEARCH_H
#define ALGORITHM_BIN_SEARCH_H

#include <algorithm> 
#include <iostream>
#include <cmath>
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/spectrum.h"

namespace algorithm {
namespace search {

class BinSearch
{
public:
    BinSearch(model::spectrum::ToleranceBy type, double tol):
        type_(type), tolerance_(tol){}

    double Tolerance() const { return tolerance_; }
    model::spectrum::ToleranceBy ToleranceType() const { return type_; }
    void set_tolerance(double tol) { tolerance_ = tol; }
    void set_tolerance_by(model::spectrum::ToleranceBy type) { type_ = type; }
    
    void Init(std::vector<double> inputs, bool sorted=false) 
    { 
        if (!sorted)
        {
            std::sort(inputs.begin(), inputs.end());
        }

        lower_ = *min_element(inputs.begin(), inputs.end()) - 1;
        upper_ = *max_element(inputs.begin(), inputs.end());

        if (type_ == model::spectrum::ToleranceBy::PPM)
            return PPMInit(inputs);
        return DaltonInit(inputs);
    }

    bool IsMatch(const double expect, const double observe)
    {
        if (type_ == model::spectrum::ToleranceBy::PPM)
        {
           return util::mass::SpectrumMass::ComputePPM(expect, observe) < tolerance_;
        }
        return std::abs(expect - observe) < tolerance_;
    }

    bool Match(double expect)
    {
        int index = Index(expect);
        if (data_[index].size() > 0)
            return true;

        int size = (int) data_.size();
        if (index < size-1 && data_[index+1].size() > 0)
        {
            if (IsMatch(expect, data_[index+1][0]))
                return true;
        }

        if (index > 0 && data_[index-1].size() > 0)
        {
            int i = (int) data_[index-1].size();
            if(IsMatch(expect, data_[index-1][i-1]))
                return true;
        } 

        return false;
    }

    int Index(double expect)
    {
        if (type_ == model::spectrum::ToleranceBy::Dalton)
            return floor((expect - lower_) / tolerance_);
        return  floor(log(expect / lower_) / log1p(tolerance_ / 1000000));
    }

protected:
    void DaltonInit(std::vector<double> inputs)
    {
        // allocate vector
        data_.clear();

        // fill the bucket
        int size = ceil((upper_ - lower_ + 1.0) / tolerance_);
        data_.assign(size, std::vector<double>());

        // assign the value
        for(const auto& it : inputs)
        {
            if (it < lower_ || it > upper_)
                continue;
            int index = floor((it - lower_) / tolerance_);
            data_[index].push_back(it);
        }
    }

    void PPMInit(std::vector<double> inputs)
    {
        // allocate vector
        data_.clear();

        // fill the bucket
        double ratio = tolerance_ / 1000000;
        int size = ceil(log(upper_ * 1.0 / lower_) / log1p(ratio));
        data_.assign(size, std::vector<double>());

        // assign the value
        for(const auto& it : inputs)
        {
            if (it < lower_ || it > upper_)
                continue;
            int index = floor(log(it / lower_) / log1p(ratio));
            data_[index].push_back(it);
        }
    }

    std::vector<std::vector<double>> data_;
    model::spectrum::ToleranceBy type_;
    double tolerance_;
    double lower_;
    double upper_;
};

} // namespace algorithm
} // namespace search 

#endif
