#include "nglycan_complex.h"

namespace model {
namespace glycan {


std::vector<std::unique_ptr<Glycan>> NGlycanComplex::Grow(Monosaccharide suger){
   std::vector<std::unique_ptr<Glycan>>  glycans;
    switch (suger)
    {   
    case Monosaccharide::GlcNAc:
        if (ValidAddGlcNAcCore()){
            std::unique_ptr<NGlycanComplex> ptr = CreateByAddGlcNAcCore();
            glycans.push_back(std::move(ptr));
        }else if (ValidAddGlcNAc()){
            if (ValidAddGlcNAcBisect()){
                std::unique_ptr<NGlycanComplex> ptr = CreateByAddGlcNAcBisect();
                glycans.push_back(std::move(ptr));
            }
            if (ValidAddGlcNAcBranch()){
                std::vector<std::unique_ptr<NGlycanComplex>> gs = CreateByAddGlcNAcBranch();
                for (auto& ptr : gs){
                    glycans.push_back(std::move(ptr));
                }
            }
        }
        break;

    case Monosaccharide::Man:
        if (ValidAddMan()){ 
            std::unique_ptr<NGlycanComplex> ptr = CreateByAddMan();
            glycans.push_back(std::move(ptr));
        }
        break;

    case Monosaccharide::Gal:
        if (ValidAddGal()){
            std::vector<std::unique_ptr<NGlycanComplex>> gs = CreateByAddGal();
            for (auto& ptr : gs){
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::Fuc:
        if (ValidAddFucCore()){
            std::unique_ptr<NGlycanComplex> ptr = CreateByAddFucCore();
            glycans.push_back(std::move(ptr));
        }
        else if (ValidAddFucTerminal()){
            std::vector<std::unique_ptr<NGlycanComplex>> gs = CreateByAddFucTerminal();
            for (auto& ptr : gs){
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::NeuAc:
        if (ValidAddNeuAc()){
            std::vector<std::unique_ptr<NGlycanComplex>> gs = CreateByAddNeuAc();
            for (auto& ptr : gs){
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::NeuGc:
        if (ValidAddNeuGc()){
            std::vector<std::unique_ptr<NGlycanComplex>> gs = CreateByAddNeuGc();
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


bool NGlycanComplex::ValidAddGlcNAcCore()
{
    return table_[0] < 2;
}

bool NGlycanComplex::ValidAddGlcNAc()
{
    return (table_[0] == 2 && table_[1] == 3);
}

std::unique_ptr<NGlycanComplex> NGlycanComplex::CreateByAddGlcNAcCore()
{
    auto g = std::make_unique<NGlycanComplex>();
    g->set_table(table_);
    g->set_table(0, table_[0]+1);
    g->set_composition(composite_);
    g->AddMonosaccharide(Monosaccharide::GlcNAc);
    return g;
}

bool NGlycanComplex::ValidAddGlcNAcBisect()
{
    //bisect 0, not extanding on GlcNAc
    return (table_[1] == 3 && table_[3] == 0 && table_[4] == 0);
}

std::unique_ptr<NGlycanComplex> NGlycanComplex::CreateByAddGlcNAcBisect(){
    auto g = std::make_unique<NGlycanComplex>();
    g->set_table(table_);
    g->set_table(3, 1);
    g->set_composition(composite_);
    g->AddMonosaccharide(Monosaccharide::GlcNAc);
    return g;
}

bool NGlycanComplex::ValidAddGlcNAcBranch(){
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 4] < table_[i + 3]) // make it order
        {
            if (table_[i + 4] == table_[i + 8] && table_[i + 12] == 0 && table_[i + 16] == 0 && table_[i + 20] == 0)
            //equal GlcNAc Gal, no Fucose attached at terminal, no terminal NeuAc, NeuGc
            {
                return true;
            }
        }
    }
    return false;

}

std::vector<std::unique_ptr<NGlycanComplex>> NGlycanComplex::CreateByAddGlcNAcBranch(){
    std::vector<std::unique_ptr<NGlycanComplex>> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 4] < table_[i + 3]) // make it order
        {
            if (table_[i + 4] == table_[i + 8] && table_[i + 12] == 0 && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                auto g = std::make_unique<NGlycanComplex>();
                g->set_table(table_);
                g->set_table(i + 4, table_[i + 4] + 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::GlcNAc);
                glycans.push_back(std::move(g));
            }
        }
    }
    return glycans;
}

bool NGlycanComplex::ValidAddMan()
{
    return (table_[0] == 2 && table_[1] < 3);
}

std::unique_ptr<NGlycanComplex> NGlycanComplex::CreateByAddMan()
{
    auto g = std::make_unique<NGlycanComplex>();
    g->set_table(table_);
    g->set_table(1, table_[1] + 1);
    g->set_composition(composite_);
    g->AddMonosaccharide(Monosaccharide::Man);
    return g;
}

bool NGlycanComplex::ValidAddGal()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
        {
            if (table_[i + 4] == table_[i + 8] + 1)
            {
                return true;
            }
        }
    }
    return false;
}
    
std::vector<std::unique_ptr<NGlycanComplex>> NGlycanComplex::CreateByAddGal(){
    std::vector<std::unique_ptr<NGlycanComplex>> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
        {
            if (table_[i + 4] == table_[i + 8] + 1)
            {
                auto g = std::make_unique<NGlycanComplex>();
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

bool NGlycanComplex::ValidAddFucCore()
{
    return (table_[0] == 1 && table_[1] == 0 && table_[2] == 0);  //core
}

std::unique_ptr<NGlycanComplex> NGlycanComplex::CreateByAddFucCore()
{
    auto g = std::make_unique<NGlycanComplex>();
    g->set_table(table_);
    g->set_table(2, 1);
    g->set_composition(composite_);
    g->AddMonosaccharide(Monosaccharide::Fuc);
    return g;
}

bool NGlycanComplex::ValidAddFucTerminal()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 12] < table_[i + 11]) // make it order
        {
            if (table_[i + 12] == 0 && table_[i + 4] > 0)
            {
                return true;
            }
        }
    }
    return false;
}
        
std::vector<std::unique_ptr<NGlycanComplex>> NGlycanComplex::CreateByAddFucTerminal()
{
    std::vector<std::unique_ptr<NGlycanComplex>> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 12] < table_[i + 11]) // make it order
        {
            if (table_[i + 12] == 0 && table_[i + 4] > 0)
            {
                auto g = std::make_unique<NGlycanComplex>();
                g->set_table(table_);
                g->set_table(i + 12, 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::Fuc);
                glycans.push_back(std::move(g));
            }
        }
    }
    return glycans;
}

bool NGlycanComplex::ValidAddNeuAc()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 16] < table_[i + 15]) // make it order
        {
            if (table_[i + 4] > 0 && table_[i + 4] == table_[i + 8] && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

std::vector<std::unique_ptr<NGlycanComplex>> NGlycanComplex::CreateByAddNeuAc()
{
    std::vector<std::unique_ptr<NGlycanComplex>> glycans;
     for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 16] < table_[i + 15]) // make it order
        {
            if (table_[i + 4] > 0 && table_[i + 4] == table_[i + 8] && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                auto g = std::make_unique<NGlycanComplex>();
                g->set_table(table_);
                g->set_table(i + 16, 1);
                g->set_composition(composite_);
                g->AddMonosaccharide(Monosaccharide::NeuAc);
                glycans.push_back(std::move(g));
            }
        }
    }
    return glycans;
}

bool NGlycanComplex::ValidAddNeuGc()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 20] < table_[i + 19]) // make it order
        {
            if (table_[i + 4] > 0 && table_[i + 4] == table_[i + 8] && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

std::vector<std::unique_ptr<NGlycanComplex>> NGlycanComplex::CreateByAddNeuGc()
{
    std::vector<std::unique_ptr<NGlycanComplex>> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 20] < table_[i + 19]) // make it order
        {
            if (table_[i + 4] > 0 && table_[i + 4] == table_[i + 8] && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                auto g = std::make_unique<NGlycanComplex>();
                g->set_table(table_);
                g->set_table(i + 20, 1);
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
