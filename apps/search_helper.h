#ifndef APP_SEARCH_HELPER_H_
#define APP_SEARCH_HELPER_H_

#include <map>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <fstream>

#include "../util/io/fasta_reader.h"
#include "../engine/protein/protein_digest.h"
#include "../engine/protein/protein_ptm.h"
#include "../engine/analysis/search_result.h"
#include "../engine/analysis/search_result.h"


// generate peptides by digestion
std::vector<engine::analysis::SearchResult>ConvertComposition(
    const std::vector<engine::analysis::SearchResult>& results, 
    const std::unordered_map<std::string, std::unique_ptr<model::glycan::Glycan>>& glycans_map)
{
    std::vector<engine::analysis::SearchResult> res;
    std::unordered_set<std::string> seen;
    for(const auto& it : results)
    {
        std::string glycan = glycans_map.find(it.Glycan())->second->Name();
        std::string key = std::to_string(it.Scan()) + "|" +  std::to_string(it.ModifySite()) + "|" +  
            it.Sequence()+ "|" + glycan;
        if (seen.find(key) == seen.end())
        {
            seen.insert(key);
            engine::analysis::SearchResult r = it;
            r.set_glycan(glycan);
            res.push_back(r);
        }     
    }
    return res;
}


std::unordered_set<std::string> PeptidesDigestion
    (const std::string& fasta_path, const SearchParameter& parameter)
{
    util::io::FASTAReader fasta_reader(fasta_path);
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
   
    engine::protein::Digestion digest;
    digest.set_miss_cleavage(parameter.miss_cleavage);
    std::unordered_set<std::string> peptides;

    // digestion
    std::deque<engine::protein::Proteases> proteases(parameter.proteases);
    engine::protein::Proteases enzyme = proteases.front();
    digest.SetProtease(enzyme);
    proteases.pop_front();
    for(const auto& protein : proteins)
    {
        std::unordered_set<std::string> seqs = digest.Sequences(protein.Sequence(),
            engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(seqs.begin(), seqs.end());
    }
        
    // double digestion or more
    while (proteases.size() > 0)
    {
        std::unordered_set<std::string> double_seqs;
        enzyme = proteases.front();
        digest.SetProtease(enzyme);
        proteases.pop_front();
        for(const auto& seq : peptides)
        {
            std::unordered_set<std::string> seqs = digest.Sequences(seq,
                engine::protein::ProteinPTM::ContainsNGlycanSite);
            double_seqs.insert(seqs.begin(), seqs.end());
        }
        peptides.insert(double_seqs.begin(), double_seqs.end());
    }   

    return peptides;
}

// report glycopeptide identification of spectrum
void ReportResults(const std::string& out_path,
    const std::vector<engine::analysis::SearchResult>&  results)
{
    std::ofstream outfile;
    outfile.open (out_path);
    outfile << "scan#,peptide,glycan,score\n";

    for(auto it : results)
    {
        outfile << it.Scan() << ",";
        outfile << it.Sequence() << ",";
        outfile << it.Glycan() << ",";
        outfile << it.Score() << "\n";
    }
    outfile.close();
}


#endif