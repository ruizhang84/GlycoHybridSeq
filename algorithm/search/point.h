#ifndef ALGORITHM_SEARCH_POINT_H_
#define ALGORITHM_SEARCH_POINT_H_

#include <vector>

namespace algorithm {
namespace search {

template <class T>
class Point
{
public:
    Point() = default;
    Point(double value,T content):
        value_(value), content_(content){}
    
    double Value() const { return value_; }
    void set_value(double value) { value_ = value; }
    T Content() { return content_; }
    void set_content(T content) { content_ = content; }
    bool operator<(const Point& other)
    {
        return value_ < other.value_;
    }

protected:
    double value_;
    T content_;
};

} // namespace algorithm
} // namespace search 

#endif
