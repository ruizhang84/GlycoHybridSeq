#ifndef ENGINE_SEARCH_GLYCAN_H_
#define ENGINE_SEARCH_GLYCAN_H_

#include <string>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "../../algorithm/search/search.h"
#include "../../model/glycan/glycan.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/ion.h"
#include "../../util/mass/spectrum.h"
#include "search_helper.h"


namespace engine{
namespace search{


class GlycanSearch
{
public:
    GlycanSearch(std::unique_ptr<algorithm::search::ISearch<std::string>> searcher):
        searcher_(std::move(searcher)){}

    // peptide seq, glycan*
    void Search(
        const std::vector<model::spectrum::Peak>& peaks, int max_charge, 
        const std::unordered_map<std::string, std::vector<model::glycan::Glycan*>>& candidate)
    {

    }
    

protected:
    std::unique_ptr<algorithm::search::ISearch<std::string>> searcher_;
};



} // namespace engine
} // namespace search

#endif