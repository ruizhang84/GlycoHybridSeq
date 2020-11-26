#ifndef UTIL_IO_FASTA_READER_H_
#define UTIL_IO_FASTA_READER_H_

#include <fstream>
#include "protein_reader.h"

namespace util {
namespace io {

class FASTAReader: public ProteinReader
{
public:
    FASTAReader(std::string path): ProteinReader(path) {}

    std::vector<model::protein::Protein> Read() override
    {
        std::vector<model::protein::Protein> result;

        std::ifstream file(path_);
        std::string line;
        model::protein::Protein protein;
        std::string seq;

        if (file.is_open()){
            while(std::getline(file, line))
            {
                // ignore comment lines
                if (line.compare(0, 1, ";") == 0)
                {
                    continue;
                }

                //e.g. >gi|186681228|ref|YP_001864424.1| phycoerythrobilin:ferredoxin oxidoreductase
                else if (line.compare(0, 1, ">") == 0)
                {
                    if (seq.length() > 0)
                    {
                        protein.set_sequence(seq);
                        result.push_back(protein);
                        seq.clear();
                    }
                    protein = model::protein::Protein();
                    protein.set_id(line);
                }
                else
                {
                    seq += trim(line);
                }
            }

            if (seq.length() > 0)
            {
                protein.set_sequence(seq);
                result.push_back(protein);
            }

        }
        return result;
    }

protected:
    std::string trim(const std::string& str)
    {
        size_t first = str.find_first_not_of(" \t\n\r\f\v");
        if (std::string::npos == first)
        {
            return str;
        }
        size_t last = str.find_last_not_of(" \t\n\r\f\v");
        return str.substr(first, (last - first + 1));
    }
};

} // namespace io
} // namespace util


#endif