#ifndef ENGINE_ANALYSIS_MULTI_COMPARISON_H_
#define ENGINE_ANALYSIS_MULTI_COMPARISON_H_

#include <algorithm>
#include <numeric>
#include <map>
#include <boost/math/distributions/normal.hpp>
#include "search_result.h"

using boost::math::normal_distribution;

namespace engine{
namespace analysis{

class MultiComparison
{
public:
    MultiComparison(double fdr): fdr_(fdr){}

    std::vector<SearchResult> Tests(
        const std::vector<SearchResult>& targets,
        const std::vector<SearchResult>& decoys)
    {
        std::vector<SearchResult> results;

         // get score lists of each spectrum
        std::map<int, std::vector<double>> v;
        for(const auto& it : decoys)
        {
            int i = it.Scan();
            if (v.find(i) == v.end())
            {
                v[i] = std::vector<double>();
            }
            v[i].push_back(it.Score());
        }

        // compute p value
        std::map<int, SearchResult> s_results;
        std::map<int, double> p_v;
        for(const auto& it : targets)
        {
            if (v.find(it.Scan()) == v.end()) // no decoys find
            {
                results.push_back(it);
            }
            else
            {
                std::vector<double>& score_list = v[it.Scan()];
                double avg = mean(score_list);
                if (avg > it.Score()) continue;

                double p = pValue(score_list, it.Score());

                p_v[it.Scan()] = p;
                s_results[it.Scan()] = it;
            }
        }

        // sort
        std::vector<double> p_values;
        for(const auto& it : p_v)
        {
            p_values.push_back(it.second);
        }
        std::sort(p_values.begin(), p_values.end());
        // find significant 
        int size = (int) p_values.size();
        double max_p = 0;
        for(int i = 0; i < size; i++)
        {
            double p = p_values[i];
            double r_p = i * 1.0 / size * fdr_;
            if (p < r_p)
            {
                max_p = std::max(max_p, p);
            }
        }

        for(const auto& it : p_v)
        {
            double p = it.second;
            if (max_p  >= p)
            {
                results.push_back(s_results[it.first]);
            }
        }
        if (!results.empty())
            std::sort(results.begin(), results.end(), OrderByScan);
        return results;
    }

protected:
    double fdr_;

    static double mean(std::vector<double> v)
    {
        int size = (int) v.size();
        double sum = 0;
        for(const auto& it : v)
        {
            sum += it;
        }
        return sum * 1.0 / size;
    }

    static double stdv(std::vector<double> v)
    {
        
        double m = 0;
        double avg = mean(v);
        int size = (int) v.size();
        for(const auto& it : v)
        {
            m += (it - avg) * (it - avg);
        }
        return std::sqrt( m * 1.0 / (size - 1));
    }

    static double pValue(std::vector<double> v, double q)
    {
        double m = mean(v);
        double s = stdv(v);
        return 0.5 * erfc( (q-m) * 1.0 / (s * std::sqrt(2)) );
    }

    static bool OrderByScan(const SearchResult& r1, const SearchResult& r2)
        { return r1.Scan() < r2.Scan(); }
};

}
}

#endif