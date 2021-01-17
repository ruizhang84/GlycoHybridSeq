#ifndef ENGINE_PROTEIN_MODIFICATION_H_
#define ENGINE_PROTEIN_MODIFICATION_H_

#include <vector>
#include <string>
#include <regex>
#include <unordered_set>

namespace engine {
namespace protein {

class Modifier
{
public:
    static std::unordered_set<std::string> DynamicModification(
        std::unordered_set<std::string>& peptides, 
        std::function<bool(const std::string&)> filter,
        bool oxidation=true, bool deamidation=true)
    {
        std::unordered_set<std::string> peptides_modified;

        if (oxidation)
            peptides = Oxidation(peptides);
        if (deamidation)
            peptides = Deamidation(peptides);

        for(const auto& it: peptides)
        {
            if(filter(it))
            {
                peptides_modified.insert(it);
            }
        }
        return peptides_modified;
    }

    static std::unordered_set<std::string> Oxidation(const std::unordered_set<std::string>& peptides)
    {
        std::unordered_set<std::string> peptides_modified;
        for(const auto& it : peptides)
        {
            std::unordered_set<std::string> seq = Oxidation(it);
            peptides_modified.insert(seq.begin(), seq.end());
        }
        return peptides_modified;
    }

    static std::unordered_set<std::string> Oxidation(const std::string& sequence)
    {
        return Modification(sequence, 'M', '$');
    }


    static std::unordered_set<std::string> Deamidation(const std::unordered_set<std::string>& peptides)
    {
        std::unordered_set<std::string> peptides_modified;
        for(const auto& it : peptides)
        {
            std::unordered_set<std::string> seq = Deamidation(it);
            peptides_modified.insert(seq.begin(), seq.end());
        }
        return peptides_modified;
    }

    static std::unordered_set<std::string> Deamidation(const std::string& sequence)
    {
        std::unordered_set<std::string> N = Modification(sequence, 'N', '@');
        std::unordered_set<std::string> Q = Modification(sequence, 'Q', '#');
        N.insert(Q.begin(), Q.end());
        return N;
    }

    static std::string Interpret(const std::string& sequence) {

        std::string s = sequence;
        std::regex m ("\\$");
        std::regex n ("@");
        std::regex q ("#");
        s = std::regex_replace (s,m,"M*");
        s = std::regex_replace (s,n,"N^");
        s = std::regex_replace (s,q,"Q^");
        return s;
    }

    static std::unordered_set<std::string> Modification(const std::string& sequence, 
        const char& origin, const char& replace)
    {
        std::unordered_set<std::string> modified;
        std::vector<int> index = FindChar(sequence, origin);
        std::vector<std::vector<int>> permutes = Combinatation(index);
        for(const auto& vec : permutes)
        {
            modified.insert(ReplaceString(sequence, vec, replace));
        }
        return modified;
    }

    static std::vector<std::vector<int>> Combinatation(const std::vector<int>& sequence)
    {
        std::vector<std::vector<int>> results;
        std::vector<int> temp;
        RecurCombination(sequence, 0, temp, results);
        return results;
    }

protected:
    static std::string ReplaceString(const std::string& sequence, const std::vector<int>& index, const char& replace)
    {
        std::string str = sequence;
        for(const auto& i : index)
        {
            str[i] = replace;
        }
        return str;
    }

    static std::vector<int> FindChar(const std::string& seq, const char& c)
    {
        std::vector<int> results;
        for(int i = 0; i < (int)seq.size(); i++)
        {
            if (seq[i] == c)
                results.push_back(i);
        }
        return results;
    }

    static void RecurCombination(const std::vector<int>& sequence, int idx, std::vector<int>& temp, std::vector<std::vector<int>>& ans)
    {
        if (idx == (int) sequence.size())
        {
            ans.push_back(temp);
            return;
        }
        temp.push_back(sequence[idx]);
        RecurCombination(sequence, idx + 1, temp, ans);
        temp.pop_back();
        RecurCombination(sequence, idx + 1, temp, ans);

    }

};


} // namespace protein
} // namespace engine 

#endif
