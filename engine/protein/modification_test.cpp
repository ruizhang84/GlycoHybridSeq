#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "modification.h"
#include <iterator>

#include <iostream>

#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../protein/protein_digest.h"
#include "../protein/protein_ptm.h"
#include "../glycan/glycan_builder.h"
#include "../../algorithm/search/bucket_search.h"
#include "../../algorithm/search/binary_search.h"

// BOOST_AUTO_TEST_CASE( Combinator ) 
// {
//     std::vector<int> seq {1, 2, 3};
//     engine::protein::Modifier modifier; 

//     std::vector<std::vector<int>> my_vectors = modifier.Combinatation(seq);

//     // print out
//     for(auto & it : my_vectors)
//     {
//         std::stringstream result;
//         std::copy(it.begin(), it.end(), std::ostream_iterator<int>(result, " "));
//         std::cout << result.str() << std::endl;
//     }
//     std::cout << my_vectors.size() << std::endl;

// }


// BOOST_AUTO_TEST_CASE( ReplaceSequence ) 
// {
   
//     engine::protein::Modifier modifier; 
//     std::string seq = "VVLHP$NYSQVD$#@NQ";
//     std::vector<std::string> result = modifier.Modification(seq, 'V', 'X');
//     for(auto & s : result)
//     {
//         std::cout << s << std::endl;
//     }
//     std::cout << modifier.Interpret(seq) << std::endl;

// }


BOOST_AUTO_TEST_CASE( digest_test ) 
{
    std::string seq = "VVLHPNYSQVD";
     
     // read fasta and build peptides
    util::io::FASTAReader fasta_reader("/home/ruiz/Documents/Glycoseq-Cpp/data/haptoglobin.fasta");
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
 
    engine::protein::Digestion digest;
    digest.SetProtease(engine::protein::Proteases::Trypsin);
    std::unordered_set<std::string> seqs = digest.Sequences(proteins.front().Sequence(),
         engine::protein::ProteinPTM::ContainsNGlycanSite);
    digest.SetProtease(engine::protein::Proteases::GluC);
    std::unordered_set<std::string> peptides;
    for(auto& it : seqs)
    {
        std::unordered_set<std::string> seq = digest.Sequences(it,
         engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(seq.begin(), seq.end());
    }

    std::unordered_set<std::string> peptides_modified = 
        engine::protein::Modifier::DynamicModification(peptides, 
        engine::protein::ProteinPTM::ContainsNGlycanSite);
    for(const auto& it : peptides_modified)
    {
        std::cout << it << std::endl;
    }

    
}