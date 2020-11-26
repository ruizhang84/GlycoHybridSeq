#ifndef ALGORITHM_BUCKET_SEARCH_H
#define ALGORITHM_BUCKET_SEARCH_H

#include <algorithm> 
#include <iostream>
#include <cmath>
#include "point.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/spectrum.h"

namespace algorithm {
namespace search {

template <class T>
class BucketSearch
{
typedef std::vector<std::shared_ptr<Point<T>>> Points;
typedef std::vector<std::vector<std::shared_ptr<Point<T>>>> Bucket;

public:
    BucketSearch(model::spectrum::ToleranceBy type, double tol, double lower, double upper):
        type_(type), tolerance_(tol), lower_(lower), upper_(upper){}
    
    void Init(Points inputs) 
    { 
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

    std::vector<T> Search(double expect)
    {
        std::vector<T> result;
        int index = Index(expect);
        if (index < 0 || index >= (int) data_.size())
            return result;
        
        for(const auto& it : data_[index])
        {
            result.push_back(it->Content());
        }

        for (int i : Neighbor(index))
        {
            for(const auto& it : data_[i])
            {
                const double observe = it->Value();
                if (IsMatch(expect, observe))
                {
                    result.push_back(it->Content());
                }
            }
        }
        return result;
    }

    bool Match(double expect)
    {
        int index = Index(expect);
        if (data_[index].size() > 0)
            return true;

        for(int i : Neighbor(index))
        {
            for(const auto& it : data_[i])
            {
                const double val = it->Value();
                if (IsMatch(expect, val))
                    return true;
            }
        }

        return false;
    }

    int Index(double expect)
    {
        if (type_ == model::spectrum::ToleranceBy::Dalton)
            return floor((expect - lower_) / tolerance_);
        return  ceil(log(expect / lower_) / log1p(tolerance_ / 1000000));
    }

protected:
    std::vector<int> Neighbor(int index)
    {
        std::vector<int> neighbors;
        if (index > 0)
            neighbors.push_back(index-1);
        if (index < (int) data_.size()-1)
            neighbors.push_back(index+1);
        return neighbors;
    }

    void DaltonInit(Points inputs)
    {
        // allocate vector
        data_.clear();
        std::sort(inputs.begin(), inputs.end(), Comp);

        // fill the bucket
        int size = ceil((upper_ - lower_ + 1.0) / tolerance_);
        data_.assign(size, std::vector<std::shared_ptr<Point<T>>>());

        // assign the value
        for(const auto& it : inputs)
        {
            double val = it->Value();
            if (val < lower_ || val > upper_)
                continue;
            int index = floor((val - lower_) / tolerance_);
            data_[index].push_back(it);
        }
    }

    void PPMInit(Points inputs)
    {
        // allocate vector
        data_.clear();
        std::sort(inputs.begin(), inputs.end(), Comp);

        // fill the bucket
        double ratio = tolerance_ / 1000000;
        int size = ceil(log(upper_ * 1.0 / lower_) / log1p(ratio));
        data_.assign(size, std::vector<std::shared_ptr<Point<T>>>());

        // assign the value
        for(const auto& it : inputs)
        {
            double val = it->Value();
            if (val < lower_ || val > upper_)
                continue;
            int index = floor(log(val / lower_) / log1p(ratio));
            data_[index].push_back(it);
        }
    }

    static bool Comp(const std::shared_ptr<Point<T>>& p1, const std::shared_ptr<Point<T>>& p2)
        { return p1->Value() < p2->Value(); }


    Bucket data_;
    model::spectrum::ToleranceBy type_;
    double tolerance_;
    double lower_;
    double upper_;
};

} // namespace algorithm
} // namespace search 

#endif
