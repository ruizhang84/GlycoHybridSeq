#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <iostream>

#include "protein_digest.h"

#include <fstream>
#include <iostream>
#include "../../util/io/fasta_reader.h"


BOOST_AUTO_TEST_CASE( Monosaccharide_test ) 
{
    std::string seq = "VVLHPNYSQVD";
    std::cout << engine::protein::ReverseNGlycopeptide(seq) << std::endl;

    // std::string seq = "MSALGAVIALLLWGQLFAVDSGNDSVTDIADDGCP"
    //                 "KPPEIAHGYVEHSVRYQCKNYYKLRTEGDGVYTLND";
    // engine::protein::Digestion digest;
    // std::unordered_set<std::string> seqs = digest.Sequences(seq, engine::protein::ProteinPTM::ContainsNGlycanSite);
    
    // for(auto& s: seqs)
    // {
    //     std::cout << s << std::endl;
    // }
}










