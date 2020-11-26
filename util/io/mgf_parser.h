#ifndef UTIL_IO_MGF_PARSER_H_
#define UTIL_IO_MGF_PARSER_H_

#include <string>
#include <map> 
#include <fstream>
#include <regex>
#include "spectrum_reader.h"

namespace util {
namespace io {

class MGFParser : public SpectrumParser
{   
public:
    MGFParser(std::string path)
        { path_ = path; }
        
    void Init() override
    {
        MGFData data;
        int scan_num = -1;

        std::ifstream file(path_);
        std::string line;

        std::smatch result;
        // std::regex start("BEGIN\\s+IONS");
        // std::regex end("END\\s+IONS");
        // std::regex title("TITLE=(.*)");
        std::regex pepmass("PEPMASS=(\\d+\\.?\\d*)");
        std::regex charge("CHARGE=(\\d+)");
        std::regex rt_second("RTINSECONDS=(\\d+\\.?\\d*)");
        std::regex scan("SCANS=(\\d+)");
        std::regex mz_intensity("^(\\d+\\.?\\d*)\\s+(\\d+\\.?\\d*)");

        if (file.is_open()){
            while(std::getline(file, line)){
                if (std::regex_search(line, result, mz_intensity))
                {
                    Peak pk(std::stod(result[1]), std::stod(result[2]));
                    data.peaks.push_back(pk);
                }
                else if (line.substr(0, 10) == "BEGIN IONS")
                {
                    data = MGFData();
                    scan_num++;
                }
                else if (line.substr(0, 8) == "END IONS")
                {
                    data_set_.emplace(scan_num, data);
                } 
                else if (line.substr(0, 6) == "TITLE=")
                {
                    data.title = line.substr(6, line.length());
                }
                else if (std::regex_search(line, result, pepmass))
                {
                    data.pep_mass = std::stod(result[1]);
                }
                else if (std::regex_search(line, result, charge)){
                    data.charge = std::stoi(result[1]);
                }
                else if (std::regex_search(line, result, scan))
                {
                    scan_num = std::stoi(result[1]);
                    data.scans = scan_num;
                }
                else if (std::regex_search(line, result, rt_second))
                {
                    data.rt_seconds = std::stod(result[1]);
                }
            }
        }
    }

    double ParentMZ(int scan_num) override 
    { 
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.pep_mass;
        }
        return 0;
    }
    int ParentCharge(int scan_num) override
    { 
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.charge;
        }
        return 0;
    }
    int GetFirstScan() override 
    {
        auto it = data_set_.begin();
        if (it != data_set_.end())
        {
            return it->first;
        }
        return -1;
    }
    int GetLastScan() override
    {
        auto it = data_set_.rbegin();
        if (it != data_set_.rend())
        {
            return it->first;
        }
        return -1;
    }
    std::vector<Peak> Peaks(int scan_num) override
    {         
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            
            return it->second.peaks;
        }
        return std::vector<Peak>();
    }

    double RTFromScanNum(int scan_num) override
    {
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.rt_seconds;
        }
        return -1;
    }
    std::string GetScanInfo(int scan_num) override
    {
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.title;
        }
        return "";
    }
    bool Exist(int scan_num) override
    {
        return data_set_.find(scan_num) != data_set_.end();
    }
    
private:
    class MGFData
    {
    public:
        MGFData() = default;

        std::vector<Peak> peaks;
        double pep_mass;
        int charge;
        double rt_seconds;
        int scans;
        std::string title;
    };
    std::map<int, MGFData> data_set_;
};


} // namespace io
} // namespace util


#endif