#ifndef UTIL_IO_SPECTRUM_READER_H_
#define UTIL_IO_SPECTRUM_READER_H_

#include <memory>
#include <string>
#include <vector>
#include "../../model/spectrum/spectrum.h"

namespace util {
namespace io {

using namespace model::spectrum;

class SpectrumParser
{
public:
    SpectrumParser() = default;
    virtual ~SpectrumParser(){}

    virtual double ParentMZ(int scan_num){ return 0; }
    virtual int ParentCharge(int scan_num){ return 0; }
    virtual int GetFirstScan(){ return 0; }
    virtual int GetLastScan(){ return 0; }

    virtual std::vector<Peak> Peaks(int scan_num){ return std::vector<Peak>(); }
    virtual std::string GetScanInfo(int scan_num){ return ""; }
    virtual double RTFromScanNum(int scan_num){ return 0; }
    virtual bool Exist(int scan_num){ return false; }
    virtual void Init(std::string path){ } 
};


class SpectrumReader
{
public:
    SpectrumReader(std::string path,
        std::unique_ptr<SpectrumParser> parser):
            path_(path), parser_(std::move(parser))
                { parser_->Init(path); }
    
    int GetFirstScan() { return parser_->GetFirstScan(); }
    int GetLastScan() { return parser_->GetLastScan(); }

    std::string GetScanInfo(int scan_num) 
    {   
        if (! parser_->Exist(scan_num))
            return "Not Exist!";
        return parser_->GetScanInfo(scan_num); 
    }

    double RTFromScanNum(int scan_num) 
    { 
        if (! parser_->Exist(scan_num))
            return -1;
        return parser_->RTFromScanNum(scan_num); 
    }
    
    Spectrum GetSpectrum(int scan_num)
    {
        Spectrum spectrum;
        if (parser_->Exist(scan_num))
        {
            std::vector<Peak> peaks = parser_->Peaks(scan_num);
            double mz = parser_->ParentMZ(scan_num);
            int charge = parser_->ParentCharge(scan_num);
            double retention = parser_->RTFromScanNum(scan_num);
            spectrum.set_peaks(peaks);
            spectrum.set_scan(scan_num);
            spectrum.set_retention(retention);
            spectrum.set_parent_mz(mz);
            spectrum.set_parent_charge(charge);
        }
        return spectrum;
    }
  
    std::vector<Spectrum> GetSpectrum(int start, int last)
    {
        std::vector<Spectrum> result;
        for (int scan_num = start; scan_num <= last; scan_num++)
        {
            if (parser_->Exist(scan_num))
                result.push_back(GetSpectrum(scan_num));
        }
        return result;
    }

    std::vector<Spectrum> GetSpectrum()
    {
        int start = GetFirstScan();
        int last = GetLastScan();
        return GetSpectrum(start, last);
    }

protected:
    std::string path_;
    std::unique_ptr<SpectrumParser> parser_;
    
};

} // namespace io
} // namespace util


#endif