#ifndef MODEL_PROTEIN_PROTEIN_H_
#define MODEL_PROTEIN_PROTEIN_H_

#include <string>

namespace model {
namespace protein {

class Protein
{
public:
    Protein() = default;
    Protein(std::string seq): seq_(seq){}
    Protein(std::string seq, std::string id):
        seq_(seq), id_(id){}

    // Mass
    double Mass() { return mass_; }
    void set_mass(double mass) { mass_ = mass; }

    std::string ID() const { return id_; }
    std::string Sequence() const { return seq_; }
    void set_id(std::string id) { id_ = id; }
    void set_sequence(std::string seq) { seq_ = seq; }

protected:
    double mass_ = -1;
    std::string seq_;
    std::string id_;
};



} // namespace protein
} // namespace model


#endif