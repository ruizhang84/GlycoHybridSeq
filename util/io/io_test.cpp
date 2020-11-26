#define BOOST_TEST_MODULE MGFParserTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <chrono> 
#include <vector>

#include "mgf_parser.h"
#include "fasta_reader.h"

namespace util {
namespace io {

BOOST_AUTO_TEST_CASE( mgf_read_test ) 
{
    auto start = std::chrono::high_resolution_clock::now(); 

    MGFParser parser("/home/ruiz/Documents/GlycoCrushSeq/data/ZC_20171218_C22_R1.mgf");
    parser.Init();
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    BOOST_CHECK( parser.GetFirstScan() == 3); 
    BOOST_CHECK( parser.ParentMZ(64) == 716.1069); 
    BOOST_CHECK( parser.ParentCharge(64) == 2); 
    BOOST_CHECK( parser.RTFromScanNum(64) == 25.5720644); 
    BOOST_CHECK( parser.GetScanInfo(64) == "C:\\Users\\Rui Zhang\\Downloads\\ZC_20171218_C22_R1.raw"); 

    Peak pk = parser.Peaks(64).front();
    BOOST_CHECK( pk.MZ() == 113.0671); 
    BOOST_CHECK( pk.Intensity() == 254.9); 

}

BOOST_AUTO_TEST_CASE( fasta_read_test ) 
{
    FASTAReader fasta_reader("/home/ruiz/Documents/GlycoCrushSeq/data/haptoglobin.fasta");
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
    BOOST_CHECK( proteins.front().Sequence() ==  
        "MSALGAVIALLLWGQLFAVDSGNDVTDIADDGCPKPPEIAHGYVEHSVRYQCKNYYKLRTEGDGVYTLND"
        "KKQWINKAVGDKLPECEADDGCPKPPEIAHGYVEHSVRYQCKNYYKLRTEGDGVYTLNNEKQWINKAVGD"
        "KLPECEAVCGKPKNPANPVQRILGGHLDAKGSFPWQAKMVSHHNLTTGATLINEQWLLTTAKNLFLNHSE"
        "NATAKDIAPTLTLYVGKKQLVEIEKVVLHPNYSQVDIGLIKLKQKVSVNERVMPICLPSKDYAEVGRVGY"
        "VSGWGRNANFKFTDHLKYVMLPVADQDQCIRHYEGSTVPEKKTPKSPVGVQPILNEHTFCAGMSKYQEDT"
        "CYGDAGSAFAVHDLEEDTWYATGILSFDKSCAVAEYGVYVKVTSIQDWVQKTIAEN"); 
}

} // namespace io
} // namespace util


