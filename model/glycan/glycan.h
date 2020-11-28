#ifndef MODEL_GLYCAN_GLYCAN_H_
#define MODEL_GLYCAN_GLYCAN_H_

#include <string>
#include <vector>
#include <map> 
#include <set>
#include <memory>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <regex>
#include <iostream>

namespace model {
namespace glycan {

enum class Monosaccharide
{ GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc};

class Glycan
{
public:
    Glycan() = default;
    virtual ~Glycan(){}

    // Mass
    double Mass() { return mass_; }
    void set_mass(double mass) { mass_ = mass; }
    // std::vector<double> Fragments() { return std::vector<double>(fragments_.begin(), fragments_.end()); }
    // std::set<double> FragmentSet() { return fragments_; } 
    // void set_fragments(std::set<double> fragments) 
    //     { fragments_ = fragments; }
    // void AddFragments(double mass)
    //     { fragments_.insert(mass); }
    
    // for print 
    std::string Name() const 
    { 
        std::string name = "";
        for (const auto& it : composite_)
        {
            switch (it.first)
            {
            case Monosaccharide::GlcNAc:
                name += "GlcNAc-" + std::to_string(it.second) + " ";
                break;
            case Monosaccharide::Man:
                name += "Man-" + std::to_string(it.second) + " ";
                break;
            case Monosaccharide::Gal:
                name += "Gal-" + std::to_string(it.second) + " ";
                break;
            case Monosaccharide::Fuc:
                name += "Fuc-" + std::to_string(it.second) + " ";
                break;    
            case Monosaccharide::NeuAc:
                name += "NeuAc-" + std::to_string(it.second) + " ";
                break;
            case Monosaccharide::NeuGc:
                name += "NeuGc-" + std::to_string(it.second) + " ";
                break;        
            default:
                break;
            }
        }
        return name;
    } 
    void set_name(const std::string& name) 
        { name_ = name; }
    
    // use as key
    std::string ID() const 
    { 
        std::stringstream result;
        std::copy(table_.begin(), table_.end(), 
            std::ostream_iterator<int>(result, " "));
        return result.str();
    }  

    std::vector<int>& Table() { return table_; }
    void set_table(const std::vector<int>& table) 
        { table_ = table; }
    void set_table(int index, int num)
    {
        if (index >= 0 && index < (int) table_.size())
            table_[index] = num;
    }

    std::map<Monosaccharide, int>&  Composition()
        { return composite_; }
    void set_composition(const std::map<Monosaccharide, int>& composite)
        { composite_ = composite; }

    const std::map<Monosaccharide, int>&  CompositionConst() const
        { return composite_; }

    virtual std::vector<std::unique_ptr<Glycan>> Grow(Monosaccharide suger)
    {
        std::vector<std::unique_ptr<Glycan>> result;
        return result;
    }

protected:
    double mass_ = -1;
    // std::set<double> fragments_;
    std::string name_;
    std::vector<int> table_;
    std::map<Monosaccharide, int> composite_; 

};


}  //  namespace glycan
}  //  namespace model

#endif

