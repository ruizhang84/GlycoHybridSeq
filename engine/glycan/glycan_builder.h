#ifndef ENGINE_GLYCAN_GLYCAN_BUILDER_H_
#define ENGINE_GLYCAN_GLYCAN_BUILDER_H_

#include <deque>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include "../../model/glycan/nglycan_complex.h"
#include "../../model/glycan/nglycan_hybrid.h"
#include "../../model/glycan/highmannose.h"
#include "../../util/mass/glycan.h"

namespace engine{
namespace glycan {

using namespace model::glycan;

class GlycanBuilder
{

public:
    GlycanBuilder(int hexNAc, int hex, int fuc, int neuAc, int neuGc):
        hexNAc_(hexNAc), hex_(hex), fuc_(fuc), neuAc_(neuAc), neuGc_(neuGc),
            candidates_({Monosaccharide::GlcNAc, Monosaccharide::Man, Monosaccharide::Gal,
                Monosaccharide::Fuc, Monosaccharide::NeuAc}){}
    virtual ~GlycanBuilder(){};

    std::unordered_map<double, std::vector<std::string>> Glycans() 
        { return glycans_; }

    std::unordered_map<std::string, std::unique_ptr<Glycan>> GlycanMaps()
        { return std::move(glycans_map_); }

    std::unordered_map<std::string, std::unique_ptr<Glycan>>& GlycanMapsRef()
        { return glycans_map_; }

    std::vector<Monosaccharide> Candidates() { return candidates_; }
    int HexNAc() { return hexNAc_; }
    int Hex() { return hex_; }
    int Fuc() { return fuc_; }
    int NeuAc() { return neuAc_; }
    int NeuGc() { return neuGc_; }
    void set_candidates(std::vector<Monosaccharide> sugers) { candidates_ = sugers; }
    void set_HexNAc(int num) { hexNAc_ = num; }
    void set_Hex(int num) { hex_ = num; }
    void set_Fuc(int num) { fuc_ = num; }
    void set_NeuAc(int num) { neuAc_ = num; }
    void set_NeuGc(int num) { neuGc_ = num; }

    virtual void Build()
    {
        std::deque<Glycan*> queue;
        InitQueue(queue);

        while (!queue.empty())
        {
            Glycan* node = queue.front();
            queue.pop_front();

            // update table id
            double mass = util::mass::GlycanMass::Compute(node->Composition());
            std::string table_id = node->ID();
            node->set_mass(mass);
            if (glycans_.find(mass) == glycans_.end())
            {
                glycans_[mass] = std::vector<std::string>();
            }
            glycans_[mass].push_back(table_id);

            // next
            for(const auto& it : candidates_)
            {
                std::vector<std::unique_ptr<Glycan>> res = node->Grow(it);

                for(auto& g : res)
                {
                    if (SatisfyCriteria(g.get()))
                    {
                        std::string id = g->ID();
                        if (glycans_map_.find(id) == glycans_map_.end())
                        {
                            // g->set_fragments(node->FragmentSet());
                            // g->AddFragments(mass);
                            node->Add(g.get()); 
                            glycans_map_[id] = std::move(g);
                            queue.push_back(glycans_map_[id].get());                         
                        }
                        else
                        {
                            node->Add(glycans_map_[id].get());
                            // glycans_map_[id]->AddFragments(mass); 
                        }
                    }
                }
            }

        }
    }

protected:
    void InitQueue(std::deque<Glycan*>& queue)
    {
        std::unique_ptr<NGlycanComplex> root = 
            std::make_unique<NGlycanComplex>();

        std::string root_id = root->ID();
        glycans_map_[root_id] = std::move(root);
        queue.push_back(glycans_map_[root_id].get());

        std::unique_ptr<NGlycanHybrid> root2 = std::make_unique<NGlycanHybrid>();
        root_id = root2->ID();
        glycans_map_[root_id] = std::move(root2);
        queue.push_back(glycans_map_[root_id].get());

        std::unique_ptr<HighMannose> root3 = std::make_unique<HighMannose>();
        root_id = root3->ID();
        glycans_map_[root_id] = std::move(root3);
        queue.push_back(glycans_map_[root_id].get());
    }

    bool SatisfyCriteria(const Glycan* glycan) const
    {
        int hexNAc = 0, hex = 0, fuc = 0, neuAc = 0, neuGc = 0;
        for(auto& it : glycan->CompositionConst())
        {
            switch (it.first)
            {
            case Monosaccharide::GlcNAc:
                hexNAc += it.second;
                break;
            case Monosaccharide::Gal:
                hex += it.second;
                break;
            case Monosaccharide::Man:
                hex += it.second;
                break;    
            case Monosaccharide::Fuc:
                fuc += it.second;
                break;   
            case Monosaccharide::NeuAc:
                neuAc += it.second;
                break;   
            case Monosaccharide::NeuGc:
                neuGc += it.second;
                break;           
            default:
                break;
            }
        }
        return (hexNAc <= hexNAc_ && hex <= hex_ && fuc <= fuc_
                && neuAc <= neuAc_ && neuGc <= neuGc_);
    }

    int hexNAc_;
    int hex_;
    int fuc_;
    int neuAc_;
    int neuGc_;
    std::unordered_map<double, std::vector<std::string>> glycans_; // glycan mass, glycan id
    std::unordered_map<std::string, std::unique_ptr<Glycan>> glycans_map_; // glycan id -> glycan
    std::vector<Monosaccharide> candidates_;

};


} // namespace engine
} // namespace glycan




#endif