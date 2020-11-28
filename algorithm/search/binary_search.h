#ifndef ALGORITHM_BINARY_SEARCH_H_
#define ALGORITHM_BINARY_SEARCH_H_

#include <vector>
#include <memory>
#include <cstdlib>
#include <algorithm> 
#include "point.h"
#include "../../model/spectrum/spectrum.h"
#include "search.h"

namespace algorithm {
namespace search {

template <class T>
class BinarySearch : public ISearch<T>
{
typedef std::vector<std::shared_ptr<Point<T>>> Points;
public:
    BinarySearch(model::spectrum::ToleranceBy type, double tol):
        type_(type), tolerance_(tol){}

    ~BinarySearch(){}

    double Tolerance() const { return tolerance_; }
    model::spectrum::ToleranceBy ToleranceType() const { return type_; }
    void set_tolerance(double tol) { tolerance_ = tol; }
    void set_tolerance_by(model::spectrum::ToleranceBy type) { type_ = type; }

    void Init(Points inputs, bool sorted=false) override
    {
        if (!inputs.empty() && !sorted)
        {
            std::sort(inputs.begin(), inputs.end(), BasicComp);
        }
        data_ = inputs;
    }

    bool IsMatch(double expect, double observe, double base)
    {
        if (type_ == model::spectrum::ToleranceBy::PPM)
        {
           return std::abs(expect - observe) / base * 1000000.0 < tolerance_;
        }
        return std::abs(expect - observe) < tolerance_;
    }
    
    std::vector<T> Search(double target, double base) override
    {
        std::vector<T> result;
        if (data_.empty()) 
            return result;

        int start = 0, end = data_.size()-1;
        while (start <= end)
        {
            int mid = (end - start) / 2 + start;
            if (IsMatch(target, data_[mid]->Value(), base))
            {
                for(int left = mid; left >= 0 && IsMatch(target, data_[left]->Value(), base); left--)
                {
                    result.push_back(data_[left]->Content());
                }

                for (int right = mid+1; right < (int) data_.size() && IsMatch(target, data_[right]->Value(), base); right++)
                {
                    result.push_back(data_[right]->Content());
                }
                break;
            }
            else if (data_[mid]->Value() < target)
                start = mid + 1;
            else
                end = mid - 1;
        }
        return result;
    }
    std::vector<T> Search(double target) override
    {
        return Search(target, target);
    }

    bool Match(double target, double base) override
    {
        if (data_.empty()) 
            return false;

        int start = 0, end = data_.size()-1;
        while (start <= end)
        {
            int mid = (end - start) / 2 + start;
            if (IsMatch(target, data_[mid]->Value(), base))
                return true;
            else if (data_[mid]->Value() < target)
                start = mid + 1;
            else
                end = mid - 1;
        }
        return false;
    }

    bool Match(double target) override
    {
        return Match(target, target);
    }

protected:
    static bool BasicComp(const std::shared_ptr<Point<T>>& p1, const std::shared_ptr<Point<T>>& p2)
        { return p1->Value() < p2->Value(); }

    model::spectrum::ToleranceBy type_;
    double tolerance_;
    Points data_;
};

} // namespace algorithm
} // namespace search 

#endif
