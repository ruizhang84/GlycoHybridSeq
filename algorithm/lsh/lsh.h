#ifndef UTIL_CALC_LSH_H
#define UTIL_CALC_LSH_H

#include <random>
#include <vector>
#include <numeric>
#include <unordered_map>
#include "../../model/spectrum/peak.h"

namespace algorithm {
namespace lsh {

class LSH
{
public:
    LSH(int size, int num): 
        weight_size_(size), hash_func_num_(num){ Init(); };
        
    void Init() 
    { 
        if (!weights_.empty())
            weights_.clear(); 
        GenHashFunc(); 
    }

    int Key(const std::unordered_map<int, model::spectrum::Peak>& embedded)
    {
        int key = 0;
        for (int i = 0; i < hash_func_num_; i++)
        {
            key <<= 1; 
            key |= RandomProjection(embedded, weights_[i]);
        }
        return key;
    }

    
    int WeightSize() { return weight_size_; }
    int HashFuncNum() { return hash_func_num_; }
    void set_weight_size(int size) { weight_size_ = size; }
    void set_hash_func_num(int num) { hash_func_num_ = num; } 

protected:
    int RandomProjection(const std::unordered_map<int, model::spectrum::Peak>& embedded, 
        const std::vector<double>& weight) const
    {
        double sum = 0.0;
        for(const auto& it : embedded)
        {
            sum += weight[it.first] * sqrt(it.second.Intensity());
        }
        return sum >=0 ? 1 : 0;
    }

    std::vector<double> GenNormalDistribWeight()
    {
        std::vector<double> weight;
        std::random_device rd{};
        std::mt19937 gen{rd()};

        std::normal_distribution<double> distribution(kMean, kSTD);

        for (int i = 0; i < weight_size_; i++) 
        {
            double number = distribution(gen);
            weight.push_back(number);
        }
        return weight;
    }

    void GenHashFunc()
    {
        for (int i = 0; i < hash_func_num_; i++)
        {
            weights_.push_back(GenNormalDistribWeight());
        }

    }

    std::vector<std::vector<double>> weights_;
    int weight_size_; 
    int hash_func_num_;
    const int kMean = 0;
    const int kSTD = 1; 
};


} // namespace io
} // namespace util


#endif