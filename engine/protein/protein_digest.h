#ifndef ENGINE_PROTEIN_PROTEIN_DIGEST_H_
#define ENGINE_PROTEIN_PROTEIN_DIGEST_H_

#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <unordered_set>
#include "protein_ptm.h"

namespace engine {
namespace protein {

enum class Proteases { Trypsin, Pepsin, Chymotrypsin, GluC };

std::string ReverseNGlycopeptide(const std::string& sequence)
{
    // reverse but last
    std::string seq = sequence;
    std::reverse(seq.begin(), seq.end()-1);

    // Asn–X–T
    int pos = ProteinPTM::FindNGlycanSite(sequence).front();
    std::swap(seq[seq.length()-2-pos], seq[seq.length()-4-pos]);
    
    return seq;
}

class Digestion
{
public:
    Digestion(): miss_cleavage_(2), min_length_(5), 
        enzyme_(Proteases::Trypsin){}

    int MissCleavage() { return miss_cleavage_; }
    int MinLength() { return min_length_; }
    Proteases Enzyme() { return enzyme_; }
    void set_min_length(int length) { min_length_ = length; }
    void set_miss_cleavage(int num) { miss_cleavage_ = num; }
    void SetProtease(Proteases enzyme) { enzyme_ = enzyme; }

    std::unordered_set<std::string> Sequences
        (const std::string seq, std::function<bool(const std::string&)> filter)
    {
        std::unordered_set<std::string> seq_list;
        std::vector<int> cutoffs = FindCutOffPosition(seq);

        //generate substring from sequences
        for (int i = 0; i <= miss_cleavage_; i++)
        {
            for (int j = 0; j <  (int) cutoffs.size() - i - 1; j++)
            {
                int start = cutoffs[j] + 1;
                int end = cutoffs[j + 1 + i];
                if (end - start + 1 >= min_length_)  // put minimum length in place
                {
                    std::string sub = seq.substr(start, end - start + 1);
                    if (filter(sub))
                        seq_list.insert(sub);
                }
            }
        }
        return seq_list;
    }

protected:
    std::vector<int> FindCutOffPosition(const std::string& sequence)
    {
        //get cleavable position, make all possible peptide cutoff  positoins
        std::vector<int> cutoffs;

        cutoffs.push_back(-1); //trivial to include starting place

       
        for (int i = 0; i < (int) sequence.length(); i++)
        {
            if (IsCleavablePosition(sequence, i))    //enzyme
            {
                cutoffs.push_back(i);
            }
        }
        if (!IsCleavablePosition(sequence, sequence.length() - 1))
        {
            cutoffs.push_back(sequence.length()- 1); //trivial to include ending place
        }

        return cutoffs;
    }
        
    bool IsCleavablePosition(const std::string& sequence, int index)
    {
        char s = sequence[index];
        switch (enzyme_)
        {
            //cleaves peptides on the C-terminal side of lysine and arginine
            case Proteases::Trypsin:
                //proline residue is on the carboxyl side of the cleavage site
                if (index < (int) sequence.length() - 1 && sequence[index + 1] == 'P')
                {
                    return false;
                }
                else if (s == 'K' || s == 'R')
                {
                    return true;
                }
                break;

            //cuts after aromatic amino acids such as phenylalanine, tryptophan, and tyrosine.
            case Proteases::Pepsin:
                if (s == 'W' || s == 'F' || s == 'Y')
                {
                    return true;
                }
                break;

            case Proteases::Chymotrypsin:
                if (index < (int) sequence.length() - 1 && sequence[index + 1] == 'P')
                {
                    return false;
                }
                else if (s == 'W' || s == 'F' || s == 'Y')
                {
                    return true;
                }
                break;

            case Proteases::GluC:
                if (index < (int) sequence.length() - 1 && sequence[index + 1] == 'P')
                {
                    return false;
                }
                else if (s == 'E' || s == 'D')
                {
                    return true;
                }
                break;
            
            default:
                break;
        }

        return false;
    }


    int miss_cleavage_;
    int min_length_;
    Proteases enzyme_;

};

} // namespace protein
} // namespace engine

#endif
