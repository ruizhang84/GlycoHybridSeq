#ifndef ALGORITHM_BUCKET_SEARCH_H_
#define ALGORITHM_BUCKET_SEARCH_H_

#include <algorithm> 
#include <iostream>
#include <cmath>
#include "point.h"
#include "../../model/spectrum/spectrum.h"
#include "search.h"

namespace algorithm {
namespace search {

template <class T>
class BucketSearch : public ISearch<T>
{
typedef std::vector<std::shared_ptr<Point<T>>> Points;
typedef std::vector<std::vector<std::shared_ptr<Point<T>>>> Bucket;

public:
    BucketSearch(model::spectrum::ToleranceBy type, double tol):
        type_(type), tolerance_(tol){}
    ~BucketSearch(){}
    
    void Init(Points inputs, bool sorted=false) override
    {
        lower_ = INT_MAX;
        upper_ = 0;
        for(const auto& it : inputs)
        {
            double val = it->Value();
            lower_ = lower_ < val ? lower_ : val;
            upper_ = upper_ > val ? upper_ : val;
        }
        lower_--;

        if (type_ == model::spectrum::ToleranceBy::PPM)
            return PPMInit(inputs);
        return DaltonInit(inputs);
    }

    void Add(std::shared_ptr<Point<T>> point)
    {
        double expect = point->Value();
        int index = Index(expect);
        if (index >= 0 && index < (int) data_.size())
            data_[index].push_back(point);
    }

    bool IsMatch(double expect, double observe, double base)
    {
        if (type_ == model::spectrum::ToleranceBy::PPM)
        {
           return std::abs(expect - observe) / base * 1000000.0 < tolerance_;
        }
        return std::abs(expect - observe) < tolerance_;
    }

    // in case when compute delta of two value,
    // to handle ppm correctly, (m1-m2)/base
    std::vector<T> Search(double expect, double base) override
    {
        std::vector<T> result;
        int index = Index(expect);
        int size = (int) data_.size();
        if (index < 0 || index >= size)
            return result;
        
        for(const auto& it : data_[index])
        {            
            result.push_back(it->Content());
        }


        if (index < size - 1)
        {
            const auto& it = data_[index+1];
            for(int i = 0; i < (int) it.size(); i++)
            {
                if (IsMatch(expect, it[i]->Value(), base))
                {    
                    result.push_back(it[i]->Content());
                }
            }
        }


        if (index > 0)
        {
            const auto& it = data_[index-1];
            int bin_size = (int) it.size();
            for(int i = bin_size-1; i >= 0; i--)
            {
                if (IsMatch(expect, it[i]->Value(), base))
                {
                    result.push_back(it[i]->Content());
                }
            }
        }
       
        return result;
    }
    std::vector<T> Search(double expect) override
    {
        return Search(expect, expect);
    }

    // base to handle ppm
    bool Match(double expect, double base) override
    {
        int index = Index(expect);

        if (index < 0 || index >= (int) data_.size())
            return false;

        if (data_[index].size() > 0)
            return true;

        int size = (int) data_.size();
        if (index < size - 1)
        {
            const auto& it = data_[index+1];
            for(int i = 0; i < (int) it.size(); i++)
            {
                if (IsMatch(expect, it[i]->Value(), base))
                {    
                    return true;
                }
            }
        }

        if (index > 0)
        {
            const auto& it = data_[index-1];
            int bin_size = (int) it.size();
            for(int i = bin_size-1; i >= 0; i--)
            {
                if (IsMatch(expect, it[i]->Value(), base))
                {
                    return true;
                }
            }
        }

        return false;
    }
    bool Match(double expect) override
    {
        return Match(expect, expect);
    }

    int Index(double expect)
    {
        if (type_ == model::spectrum::ToleranceBy::Dalton)
            return floor((expect - lower_) / tolerance_);
        double ratio = 1.0/(1.0 - tolerance_ / 1000000);
        return  floor(log(expect * 1.0/ lower_) / log(ratio));
    }

protected:
    void DaltonInit(Points inputs)
    {
        // allocate vector
        data_.clear();

        // fill the bucket
        int size = ceil((upper_ - lower_ + 1.0) / tolerance_);
        size = std::max(size, 0);
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

        // fill the bucket
        double ratio = 1.0/(1.0 - tolerance_ / 1000000);
        int size = ceil(log(upper_ / lower_) / log(ratio));
        size = std::max(size, 0);
        data_.assign(size, std::vector<std::shared_ptr<Point<T>>>());

        // assign the value
        for(const auto& it : inputs)
        {
            double val = it->Value();
            if (val < lower_ || val > upper_)
                continue;
            int index = floor(log(val / lower_) / log(ratio));    
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
