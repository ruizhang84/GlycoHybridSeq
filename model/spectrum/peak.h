#ifndef MODEL_SPECTRUM_SPECTRUM_PEAK_H_
#define MODEL_SPECTRUM_SPECTRUM_PEAK_H_

namespace model {
namespace spectrum {

class Peak
{
public:
    Peak() = default;
    Peak(double mz, double intensity):
        mz_(mz), intensity_(intensity){}

    Peak(const Peak& other)
    {
        mz_ = other.mz_;
        intensity_ = other.intensity_;
    }
    Peak& operator=(const Peak& other)
    {
        mz_ = other.mz_;
        intensity_ = other.intensity_;
        return *this;
    }

    double MZ() const { return mz_; }
    void set_mz(double mz) { mz_ = mz; }
    double Intensity() const { return intensity_; }
    void set_intensity(double intensity) 
        { intensity_ = intensity; }
    
    bool operator<(const Peak& other) const
        { return mz_ < other.mz_; }

protected:
    double mz_;
    double intensity_;
};



} // namespace spectrum
} // namespace model


#endif