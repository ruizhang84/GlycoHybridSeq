#ifndef MODEL_GLYCAN_NGLYCAN_HYBRID_H_
#define MODEL_GLYCAN_NGLYCAN_HYBRID_H_
// The hybrid type of n-glycan

#include <sstream>
#include <iterator>
#include <typeinfo>
#include "glycan.h"

//GlcNAc(2) - Man(3) - Fuc(1) - GlcNAc(bisect,1) -0,1,2,3
//[Man(branch1) - Man(branch2)]  4, 5
//[GlcNAc(branch1) - GlcNAc(branch2)] 6, 7
//[Gal(branch1) - Gal(branch2)]   8, 9
//[Fuc, Fuc]         10, 11
//[NeuAc, NeuAc]     12, 13
//[NeuGc, NeuGc]     14, 15


namespace model {
namespace glycan {
class NGlycanHybrid : public Glycan 
{
public:
    NGlycanHybrid()
    { 
        table_.assign(16, 0);
    }
    ~NGlycanHybrid(){}

    std::string Type() override { return "Hybrid"; }
    
    std::vector<std::unique_ptr<Glycan>> Grow(Monosaccharide suger) override;

protected:
    void AddMonosaccharide(Monosaccharide suger)
    {
        auto it = composite_.find(suger);
        if (it != composite_.end())
        {
            composite_[suger] += 1;
        }
        else
        {
            composite_[suger] = 1;
        }
    }

    bool ValidAddGlcNAcCore();
    std::unique_ptr<NGlycanHybrid> CreateByAddGlcNAcCore();
    bool ValidAddGlcNAc();
    bool ValidAddGlcNAcBisect();
    std::unique_ptr<NGlycanHybrid> CreateByAddGlcNAcBisect();
    bool ValidAddGlcNAcBranch();
    std::vector<std::unique_ptr<NGlycanHybrid>> CreateByAddGlcNAcBranch();

    bool ValidAddManCore();
    std::unique_ptr<NGlycanHybrid> CreateByAddManCore();
    bool ValidAddManBranch();
    std::vector<std::unique_ptr<NGlycanHybrid>> CreateByAddManBranch();

    bool ValidAddGal();
    std::vector<std::unique_ptr<NGlycanHybrid>> CreateByAddGal();

    bool ValidAddFucCore();
    std::unique_ptr<NGlycanHybrid> CreateByAddFucCore();

    bool ValidAddFucTerminal();
    std::vector<std::unique_ptr<NGlycanHybrid>> CreateByAddFucTerminal();

    bool ValidAddNeuAc();
    std::vector<std::unique_ptr<NGlycanHybrid>> CreateByAddNeuAc();

    bool ValidAddNeuGc();
    std::vector<std::unique_ptr<NGlycanHybrid>> CreateByAddNeuGc();

}; 

}  //  namespace glycan
}  //  namespace model


#endif