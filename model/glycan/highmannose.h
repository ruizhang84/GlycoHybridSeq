#ifndef MODEL_GLYCAN_NGLYCAN_HIGH_MANNOSE_H_
#define MODEL_GLYCAN_NGLYCAN_HIGH_MANNOSE_H_
// The high mannose type of n-glycan

#include <sstream>
#include <iterator>
#include <typeinfo>
#include "glycan.h"

// GlcNAc(2) - Man(3) - Fuc 0 1 2
// [Man(branch1) - Man(branch2) - Man(branch3)]  3 4 5

namespace model {
namespace glycan {
class HighMannose : public Glycan 
{
public:
    HighMannose()
    { 
        table_.assign(6, 0);
    }
    ~HighMannose(){}

    std::string Type() { return "High Mannose"; }

    std::vector<std::unique_ptr<Glycan>> Grow(Monosaccharide suger) override
    {
        std::vector<std::unique_ptr<Glycan>>  glycans;
        switch (suger)
        {   
        case Monosaccharide::GlcNAc:
            if (ValidAddGlcNAc()){
                std::unique_ptr<HighMannose> ptr = CreateByAddGlcNAc();
                glycans.push_back(std::move(ptr));
            }
            break;

        case Monosaccharide::Man:
            if (ValidAddManCore()){ 
                std::unique_ptr<HighMannose> ptr = CreateByAddManCore();
                glycans.push_back(std::move(ptr));
            }else if (ValidAddManBranch())
            {
                std::vector<std::unique_ptr<HighMannose>> gs = CreateByAddManBranch();
                for (auto& ptr : gs){
                    glycans.push_back(std::move(ptr));
                }
            }
            break;

        case Monosaccharide::Fuc:
            if (ValidAddFucCore()){
                std::unique_ptr<HighMannose> ptr = CreateByAddFucCore();
                glycans.push_back(std::move(ptr));
            }
            break;

        default:
            break;
        }
        return glycans;
    }

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

    bool ValidAddGlcNAc()
    {
        return table_[0] < 2;
    }
    std::unique_ptr<HighMannose> CreateByAddGlcNAc()
    {
        auto g = std::make_unique<HighMannose>();
        g->set_table(table_);
        g->set_table(0, table_[0]+1);
        g->set_composition(composite_);
        g->AddMonosaccharide(Monosaccharide::GlcNAc);
        return g;
    }

    bool ValidAddManCore()
     {
        if (table_[0] == 2 && table_[1] < 3)
            return true;
        return false;
    }

    std::unique_ptr<HighMannose> CreateByAddManCore()
    {
        auto g = std::make_unique<HighMannose>();
        g->set_table(table_);
        g->set_table(1, table_[1]+1);
        g->set_composition(composite_);
        g->AddMonosaccharide(Monosaccharide::Man);
        return g;
    }

    bool ValidAddManBranch()
    {
        if (table_[0] == 2 && table_[1] == 3)
        {
            return true;
        }
        return false;
    }

    std::vector<std::unique_ptr<HighMannose>> CreateByAddManBranch()
    {    
        std::vector<std::unique_ptr<HighMannose>> glycans;
        for (int i = 0; i < 3; i++)
        {
            if (i == 0 || table_[i + 3] < table_[i + 2]) // make it order
            {
                
                auto g = std::make_unique<HighMannose>();
                g->set_table(table_);
                g->set_table(i + 3, table_[i + 3] + 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::Man);
                glycans.push_back(std::move(g));
            }
        }
        return glycans;
    }


    bool ValidAddFucCore()
    {
        if (table_[1] == 0 && table_[2] == 0)
        {
            return true;
        }
        return false;
    }
    std::unique_ptr<HighMannose> CreateByAddFucCore()
    {
        auto g = std::make_unique<HighMannose>();
        g->set_table(table_);
        g->set_table(2, 1);
        g->set_composition(composite_);
        g->AddMonosaccharide(Monosaccharide::Fuc);
        return g;
    }


}; 

}  //  namespace glycan
}  //  namespace model


#endif