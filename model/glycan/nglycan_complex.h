#ifndef MODEL_GLYCAN_NGLYCAN_COMPLEX_H_
#define MODEL_GLYCAN_NGLYCAN_COMPLEX_H_
// The complex type of n-glycan

#include <sstream>
#include <iterator>
#include <typeinfo>
#include "glycan.h"

//GlcNAc(2) - Man(3) - Fuc(1) - GlcNAc(bisect,1) -0,1,2,3
//[GlcNAc(branch1) - GlcNAc(branch2) - GlcNAc(branch3) - GlcNAc(branch4)] -4,5,6,7
//[Gal(branch1) - Gal(branch2) - Gal(branch3) - Gal(branch4)] -8,9,10,11
//[Fuc(branch1) - Fuc(branch2) - Fuc(branch3) - Fuc(branch4)] -12,13,14,15
//[NeuAc(branch1) - NeuAc(branch2) - NeuAc(branch3) - NeuAc(branch4)] -16,17,18,19
//[NeuGc(branch1) - NeuGc(branch2) - NeuGc(branch3) - NeuGc(branch4)] -20,21,22,23

namespace model {
namespace glycan {
class NGlycanComplex : public Glycan 
{
public:
    NGlycanComplex()
    { 
        table_.assign(24, 0);
    }
    ~NGlycanComplex(){}

    GlycanType Type() override { return GlycanType::NComplex; }
    
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
    std::unique_ptr<NGlycanComplex> CreateByAddGlcNAcCore();
    bool ValidAddGlcNAc();
    bool ValidAddGlcNAcBisect();
    std::unique_ptr<NGlycanComplex> CreateByAddGlcNAcBisect();
    bool ValidAddGlcNAcBranch();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddGlcNAcBranch();

    bool ValidAddMan();
    std::unique_ptr<NGlycanComplex> CreateByAddMan();

    bool ValidAddGal();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddGal();

    bool ValidAddFucCore();
    std::unique_ptr<NGlycanComplex> CreateByAddFucCore();

    bool ValidAddFucTerminal();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddFucTerminal();

    bool ValidAddNeuAc();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddNeuAc();

    bool ValidAddNeuGc();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddNeuGc();

}; 

}  //  namespace glycan
}  //  namespace model


#endif