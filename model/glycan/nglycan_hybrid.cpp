#include "nglycan_hybrid.h"

namespace model {
namespace glycan {


std::vector<std::unique_ptr<Glycan>> NGlycanHybrid::Grow(Monosaccharide suger){
   std::vector<std::unique_ptr<Glycan>>  glycans;
    switch (suger)
    {   
    case Monosaccharide::GlcNAc:
        if (ValidAddGlcNAcCore()){
            std::unique_ptr<NGlycanHybrid> ptr = CreateByAddGlcNAcCore();
            glycans.push_back(std::move(ptr));
        }else if (ValidAddGlcNAc()){
            if (ValidAddGlcNAcBisect()){
                std::unique_ptr<NGlycanHybrid> ptr = CreateByAddGlcNAcBisect();
                glycans.push_back(std::move(ptr));
            }
            if (ValidAddGlcNAcBranch()){
                std::vector<std::unique_ptr<NGlycanHybrid>> gs = CreateByAddGlcNAcBranch();
                for (auto& ptr : gs){
                    glycans.push_back(std::move(ptr));
                }
            }
        }
        break;

    case Monosaccharide::Man:
        if (ValidAddManCore()){ 
            std::unique_ptr<NGlycanHybrid> ptr = CreateByAddManCore();
            glycans.push_back(std::move(ptr));
        }else if (ValidAddManBranch())
        {
            std::vector<std::unique_ptr<NGlycanHybrid>> gs = CreateByAddManBranch();
            for (auto& ptr : gs){
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::Gal:
        if (ValidAddGal()){
            std::vector<std::unique_ptr<NGlycanHybrid>> gs = CreateByAddGal();
            for (auto& ptr : gs){
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::Fuc:
        if (ValidAddFucCore()){
            std::unique_ptr<NGlycanHybrid> ptr = CreateByAddFucCore();
            glycans.push_back(std::move(ptr));
        }
        if (ValidAddFucTerminal()){
            std::vector<std::unique_ptr<NGlycanHybrid>> gs = CreateByAddFucTerminal();
            for (auto& ptr : gs){
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::NeuAc:
        if (ValidAddNeuAc()){
            std::vector<std::unique_ptr<NGlycanHybrid>> gs = CreateByAddNeuAc();
            for (auto& ptr : gs){
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::NeuGc:
        if (ValidAddNeuGc()){
            std::vector<std::unique_ptr<NGlycanHybrid>> gs = CreateByAddNeuGc();
            for (auto& ptr : gs){
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    default:
        break;
    }
    return glycans;
}


bool NGlycanHybrid::ValidAddGlcNAcCore()
{
    return table_[0] < 2;
}

bool NGlycanHybrid::ValidAddGlcNAc()
{
    return (table_[0] == 2 && table_[1] == 3);
}

std::unique_ptr<NGlycanHybrid> NGlycanHybrid::CreateByAddGlcNAcCore()
{
    auto g = std::make_unique<NGlycanHybrid>();
    g->set_table(table_);
    g->set_table(0, table_[0]+1);
    g->set_composition(composite_);
    g->AddMonosaccharide(Monosaccharide::GlcNAc);
    return g;
}

bool NGlycanHybrid::ValidAddGlcNAcBisect()
{
    //bisect 0, not extanding on GlcNAc
    return (table_[1] == 3 && table_[3] == 0 && table_[4] == 0);
}

std::unique_ptr<NGlycanHybrid> NGlycanHybrid::CreateByAddGlcNAcBisect(){
    auto g = std::make_unique<NGlycanHybrid>();
    g->set_table(table_);
    g->set_table(3, 1);
    g->set_composition(composite_);
    g->AddMonosaccharide(Monosaccharide::GlcNAc);
    return g;
}

bool NGlycanHybrid::ValidAddGlcNAcBranch(){
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 6] < table_[i + 5]) // make it order
        {
            if (table_[i + 6] == table_[i + 8] && table_[i + 10] == 0 && table_[i + 12] == 0 && table_[i + 14] == 0)
            //equal GlcNAc Gal, no Fucose attached at terminal, no terminal NeuAc, NeuGc
            {
                return true;
            }
        }
    }
    return false;

}

std::vector<std::unique_ptr<NGlycanHybrid>> NGlycanHybrid::CreateByAddGlcNAcBranch(){
    std::vector<std::unique_ptr<NGlycanHybrid>> glycans;
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 6] < table_[i + 5]) // make it order
        {
            if (table_[i + 6] == table_[i + 8] && table_[i + 10] == 0 && table_[i + 12] == 0 && table_[i + 14] == 0)
            {
                auto g = std::make_unique<NGlycanHybrid>();
                g->set_table(table_);
                g->set_table(i + 6, table_[i + 6] + 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::GlcNAc);
                glycans.push_back(std::move(g));
            }
        }
    }
    return glycans;
}


bool NGlycanHybrid::ValidAddManCore()
{
    if (table_[0] == 2 && table_[1] < 3)
        return true;
    return false;
}

std::unique_ptr<NGlycanHybrid> NGlycanHybrid::CreateByAddManCore()
{
    auto g = std::make_unique<NGlycanHybrid>();
    g->set_table(table_);
    g->set_table(1, table_[1] + 1);
    g->set_composition(composite_);
    g->AddMonosaccharide(Monosaccharide::Man);
    return g;
}

bool NGlycanHybrid::ValidAddManBranch()
{
    if (table_[0] == 2 && table_[1] == 3)
    {
        return true;
    }
    return false;
}

std::vector<std::unique_ptr<NGlycanHybrid>> NGlycanHybrid::CreateByAddManBranch(){
    std::vector<std::unique_ptr<NGlycanHybrid>> glycans;
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 4] < table_[i + 3]) // make it order
        {
            
            auto g = std::make_unique<NGlycanHybrid>();
            g->set_table(table_);
            g->set_table(i + 4, table_[i + 4] + 1);
            g->set_composition(composite_);
            g->AddMonosaccharide(Monosaccharide::Man);
            glycans.push_back(std::move(g));
        }
    }
    return glycans;
}


bool NGlycanHybrid::ValidAddGal()
{
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
        {
            if (table_[i + 6] == table_[i + 8] + 1)
            {
                return true;
            }
        }
    }
    return false;
}
    
std::vector<std::unique_ptr<NGlycanHybrid>> NGlycanHybrid::CreateByAddGal(){
    std::vector<std::unique_ptr<NGlycanHybrid>> glycans;
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
        {
            if (table_[i + 6] == table_[i + 8] + 1)
            {
                auto g = std::make_unique<NGlycanHybrid>();
                g->set_table(table_);
                g->set_table(i + 8, table_[i + 8] + 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::Gal);
                glycans.push_back(std::move(g));
            }
        }
    }
    return glycans;
}

bool NGlycanHybrid::ValidAddFucCore()
{
    return (table_[2] == 0); 
    // return (table_[0] == 1 && table_[1] == 0 && table_[2] == 0);  //core
}

std::unique_ptr<NGlycanHybrid> NGlycanHybrid::CreateByAddFucCore()
{
    auto g = std::make_unique<NGlycanHybrid>();
    g->set_table(table_);
    g->set_table(2, 1);
    g->set_composition(composite_);
    g->AddMonosaccharide(Monosaccharide::Fuc);
    return g;
}

bool NGlycanHybrid::ValidAddFucTerminal()
{
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 10] < table_[i + 9]) // make it order
        {
            if (table_[i + 10] == 0 && table_[i + 6] > 0)
            {
                return true;
            }
        }
    }
    return false;
}
        
std::vector<std::unique_ptr<NGlycanHybrid>> NGlycanHybrid::CreateByAddFucTerminal()
{
    std::vector<std::unique_ptr<NGlycanHybrid>> glycans;
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 10] < table_[i + 9]) // make it order
        {
            if (table_[i + 10] == 0 && table_[i + 6] > 0)
            {
                auto g = std::make_unique<NGlycanHybrid>();
                g->set_table(table_);
                g->set_table(i + 10, 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::Fuc);
                glycans.push_back(std::move(g));
            }
        }
    }
    return glycans;
}

bool NGlycanHybrid::ValidAddNeuAc()
{
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 12] < table_[i + 11]) // make it order
        {
            if (table_[i + 6] > 0 && table_[i + 6] == table_[i + 8] && table_[i + 12] == 0 && table_[i + 14] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

std::vector<std::unique_ptr<NGlycanHybrid>> NGlycanHybrid::CreateByAddNeuAc()
{
    std::vector<std::unique_ptr<NGlycanHybrid>> glycans;
     for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 12] < table_[i + 11]) // make it order
        {
            if (table_[i + 6] > 0 && table_[i + 6] == table_[i + 8] && table_[i + 12] == 0 && table_[i + 14] == 0)
            {
                auto g = std::make_unique<NGlycanHybrid>();
                g->set_table(table_);
                g->set_table(i + 12, 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::NeuAc);
                glycans.push_back(std::move(g));
            }
        }
    }
    return glycans;
}

bool NGlycanHybrid::ValidAddNeuGc()
{
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 14] < table_[i + 13]) // make it order
        {
            if (table_[i + 6] > 0 && table_[i + 6] == table_[i + 8] && table_[i + 12] == 0 && table_[i + 14] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

std::vector<std::unique_ptr<NGlycanHybrid>> NGlycanHybrid::CreateByAddNeuGc()
{
    std::vector<std::unique_ptr<NGlycanHybrid>> glycans;
    for (int i = 0; i < 2; i++)
    {
        if (i == 0 || table_[i + 14] < table_[i + 13]) // make it order
        {
            if (table_[i + 6] > 0 && table_[i + 6] == table_[i + 8] && table_[i + 12] == 0 && table_[i + 14] == 0)
            {
                auto g = std::make_unique<NGlycanHybrid>();
                g->set_table(table_);
                g->set_table(i + 14, 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::NeuGc);
                glycans.push_back(std::move(g));
            }
        }
    }
    return glycans;
}


}  //  namespace glycan
}  //  namespace model
