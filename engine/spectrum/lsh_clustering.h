#ifndef ENGINE_SPECTRUM_LSH_CLUSTERING_H_
#define ENGINE_SPECTRUM_LSH_CLUSTERING_H_

#include "binpacking.h"
#include "spectrum_sim.h"
#include "../../model/spectrum/spectrum.h"
#include "../../algorithm/lsh/lsh.h"

namespace engine {
namespace spectrum {

class EmbededSpectrum
{
public:
    EmbededSpectrum(model::spectrum::Spectrum spectrum):
        spectrum_(spectrum){};

    model::spectrum::Spectrum Spectrum()
        { return spectrum_; }

    const std::unordered_map<int, model::spectrum::Peak> Embed() const
        { return embedded_; }
    void set_embed(std::unordered_map<int, model::spectrum::Peak> embedded)
        { embedded_ = embedded;}

    int Scan() { return spectrum_.Scan(); }
    void set_scan(int scan) { spectrum_.set_scan(scan); }
    int Retention() { return spectrum_.Retention(); }
    void set_retention(double retention) { spectrum_.set_retention(retention); }

    std::vector<model::spectrum::Peak> Peaks() 
        { return spectrum_.Peaks(); }
    const std::vector<model::spectrum::Peak>& Peaks() const 
        { return spectrum_.Peaks(); }
    void set_peaks(std::vector<model::spectrum::Peak>& peaks) 
        { spectrum_.set_peaks(peaks); }

    double PrecursorMZ() { return spectrum_.PrecursorMZ(); }
    double PrecursorCharge() { return spectrum_.PrecursorCharge(); }

    void set_parent_mz(double mz) { spectrum_.set_parent_mz(mz); }
    void set_parent_charge(int charge) { spectrum_.set_parent_charge(charge); }

protected:
    model::spectrum::Spectrum spectrum_;
    std::unordered_map<int, model::spectrum::Peak> embedded_;
};

class LSHClustering
{
public:
    LSHClustering(model::spectrum::ToleranceBy type, double tol, double lower=200.0, double upper=2000, int hash_num=16):
        binpacker_(BinPacking(type, tol, lower, upper)), 
        lsh_(algorithm::lsh::LSH(binpacker_.Size(), hash_num)){}

    std::vector<EmbededSpectrum> Embed(std::vector<model::spectrum::Spectrum> spectra)
    {
        std::vector<EmbededSpectrum> embeds;
        for(const auto& spectrum : spectra)
        {
            EmbededSpectrum embed(spectrum);
            embed.set_embed(binpacker_.Packing(spectrum));
            embeds.push_back(embed);
        }
        return embeds;
    }

    std::unordered_map<int, std::vector<EmbededSpectrum>> Hashing(std::vector<EmbededSpectrum> embeds)
    {
        std::unordered_map<int, std::vector<EmbededSpectrum>> hash_table;
        for(const auto& embed : embeds)
        {
            int key = lsh_.Key(embed.Embed());
            if (hash_table.find(key) == hash_table.end())
            {
                hash_table[key] = std::vector<EmbededSpectrum>();
            }
            hash_table[key].push_back(embed);
        }
        return hash_table;
    }





    std::vector<model::spectrum::Spectrum> filter(std::unordered_map<int, std::vector<EmbededSpectrum>> hash_table)
    {
        std::vector<model::spectrum::Spectrum> results;
        for(auto& it : hash_table)
        {
            if (it.second.size() > 3)
            {
                for(auto& s : it.second)
                {
                    results.push_back(s.Spectrum());
                }
            }
        }
        return results;
    }


protected:
    BinPacking binpacker_;
    algorithm::lsh::LSH lsh_;
    SpectrumSim sim_;
};


} // namespace spectrum
} // namespace engine


#endif