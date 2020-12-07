#ifndef ENGINE_ANALYSIS_SEARCH_RESULT_H_
#define ENGINE_ANALYSIS_SEARCH_RESULT_H_


#include "../../model/spectrum/peak.h"
#include "../../model/glycan/glycan.h"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath> 
#include <memory>
#include <algorithm>



namespace engine{
namespace analysis{

class SearchResult
{
public:
    SearchResult() = default;

    int Scan() const { return scan_; }
    double Retention() const { return retention_; }
    int ModifySite() const { return pos_; }
    std::string Sequence() const { return peptide_; }
    std::string Glycan() const { return glycan_; }
    double Score() const { return score_; }

    void set_scan(int scan) { scan_ = scan; }
    void set_retention(double retention) { retention_ = retention; }
    void set_site(int pos) { pos_ = pos; }
    void set_peptide(std::string seq) { peptide_ = seq; }
    void set_glycan(std::string glycan) { glycan_ = glycan; }
    void set_score(double score) { score_ = score; }

protected:
    int scan_;
    double retention_;
    std::string peptide_;
    std::string glycan_;
    int pos_;
    double score_;
};


} // namespace analysis
} // namespace search

#endif