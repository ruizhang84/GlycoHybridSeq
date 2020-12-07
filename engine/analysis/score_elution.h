#ifndef ENGINE_SCORE_COELUTION_H
#define ENGINE_SCORE_COELUTION_H

#include <climits> 
#include "search_result.h"
#include "../../util/io/spectrum_reader.h"

namespace engine{
namespace analysis {

class CoElution
{
public:
    CoElution(double range=1.0): range_(range){}

    void Update(std::vector<engine::analysis::SearchResult>& results)
    {
        // compute coelution of peptide sequence
        std::unordered_map<std::string, int> total;
        double start = INT_MAX,  end = 0;
        for(const auto& r : results)
        {
            double time = r.Retention();
            start = std::min(start, time);
            end = std::max(end, time);
            if (total.find(r.Sequence()) == total.end())
            {
                total[r.Sequence()] = 0;
            }
            total[r.Sequence()] += 1;
        }

        // count peptide sequence in each range bucket
        int size = ceil((end - start + 1.0) / range_);
        std::vector<std::unordered_map<std::string, int>> range_count;
        range_count.assign(size, std::unordered_map<std::string, int>());
        for(const auto& r : results)
        {
            double time = r.Retention();
            int index = Index(time, start, end);
            std::string peptide = r.Sequence();

            std::unordered_map<std::string, int>& counts = range_count[index];
            
            if (counts.find(peptide) == counts.end())
            {
                counts[peptide] = 0;
            }
            counts[peptide] += 1;
        }

        // udpate score
        for(auto& it : results)
        {
            
            double time = it.Retention();
            int index = Index(time, start, end);
            std::string s = it.Sequence();

            int count = range_count[index][s];
            // find neighbor counts
            if (index > 0)
            {
                std::unordered_map<std::string, int>& counts = range_count[index-1];
                if (counts.find(s) != counts.end())
                    count += counts[s];
            }
            if (index < size - 1)
            {
                std::unordered_map<std::string, int>& counts = range_count[index+1];
                if (counts.find(s) != counts.end())
                    count += counts[s];
            }

            double score = count * 1.0 / total[s];
            if (total[s] == 1)
                score = 0.3;
            it.set_score(it.Score() * score);
        }
    }


protected:
    double range_ = 1.0;
    int Index(double retention, int start, int end)
    {
        return floor(retention - start) * range_/ (end - start + 1);
    }
};


}   // namespace score 
}   // namespace engine



#endif