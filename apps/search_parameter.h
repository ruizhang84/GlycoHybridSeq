#ifndef APP_SEARCH_PARAMETER_H_
#define APP_SEARCH_PARAMETER_H_

#include <deque>

#include "../model/spectrum/spectrum.h"
#include "../engine/protein/protein_digest.h"

struct SearchParameter
{

    // upper bound of glycan seaerch
    int n_thread = 6;
    int hexNAc_upper_bound = 12;
    int hex_upper_bound = 12;
    int fuc_upper_bound = 5;
    int neuAc_upper_bound = 4;
    int neuGc_upper_bound = 0;
    // searching precision
    double ms1_tol = 10;
    model::spectrum::ToleranceBy ms1_by =
        model::spectrum::ToleranceBy::PPM;
    double ms2_tol = 0.01;
    model::spectrum::ToleranceBy ms2_by = 
        model::spectrum::ToleranceBy::Dalton;
    // fdr
    double fdr_rate = 0.01;
    // protease
    std::deque<engine::protein::Proteases> proteases
    {
        engine::protein::Proteases::Trypsin,
        engine::protein::Proteases::GluC
    };
    int miss_cleavage = 2;
    // glycan type
    bool complex = true;
    bool hybrid = false;
    bool highmannose = false;
    // dynamic modification
    bool oxidation = false;
    bool deamidation = false;

};



#endif