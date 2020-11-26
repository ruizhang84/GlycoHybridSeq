#ifndef UTIL_IO_PROTEIN_READER_H_
#define UTIL_IO_PROTEIN_READER_H_

#include <string>
#include <vector>
#include "../../model/protein/protein.h"

namespace util {
namespace io {

class ProteinReader
{
public:
    ProteinReader(std::string path): path_(path){}
    virtual ~ProteinReader(){}
    
    virtual std::vector<model::protein::Protein> Read()
    {
         std::vector<model::protein::Protein> proteins;
         return proteins;
    }

    std::string Path() const { return path_; }
    void set_path(std::string path) { path_ = path; }

protected:
    std::string path_;

};

} // namespace io
} // namespace util


#endif