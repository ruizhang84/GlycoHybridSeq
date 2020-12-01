#ifndef ALGORITHM_SEARCH_H_
#define ALGORITHM_SEARCH_H_

#include <algorithm> 
#include <iostream>
#include <cmath>
#include "point.h"
#include "../../model/spectrum/spectrum.h"

namespace algorithm {
namespace search {

template <class T>
class ISearch
{
typedef std::vector<std::shared_ptr<Point<T>>> Points;
public:
    virtual ~ISearch(){}
    virtual void Init(Points inputs, bool sorted=false){}
    virtual std::vector<T> Search(double expect, double base) {return std::vector<T>();}
    virtual std::vector<T> Search(double expect) {return std::vector<T>();}
    virtual bool Match(double expect, double base) {return false;}
    virtual bool Match(double expect) {return false;}

};

} // namespace algorithm
} // namespace search 

#endif
