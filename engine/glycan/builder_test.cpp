#define BOOST_TEST_MODULE BuilderTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <chrono> 
#include "glycan_builder.h"


namespace engine{
namespace glycan {

BOOST_AUTO_TEST_CASE( glycan_builder_test ) 
{
    GlycanBuilder builder(2, 3, 1, 1, 0);

    auto start = std::chrono::high_resolution_clock::now(); 
    builder.Build();
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    std::unordered_map<double, std::vector<std::string>> glycans = builder.Glycans();
    std::unordered_map<std::string, std::unique_ptr<Glycan>> glycan_map = builder.GlycanMaps();

    for(const auto& it : glycan_map)
    {
        std::cout << it.first << std::endl;
        for(const auto& g : it.second->Children())
        {
            std::cout << g->Name() << std::endl;
        }
        std::cout << std::endl;
    }

    // std::cout << glycan_map["1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "]->Name() << std::endl;
    
    // for(const auto& it : glycans)
    // {
    //     std::cout << it.first << ": " << std::endl;
        
    //     for(const auto& id : it.second)
    //     {
    //         std::cout << id << std::endl;
    //         BOOST_CHECK_CLOSE(it.first, glycan_map[id]->Mass(), 0.001);
    //         for(double fg : glycan_map[id]->Fragments())
    //         {
    //             std::cout << fg << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    std::cout << glycan_map.size() << std::endl;
    // BOOST_CHECK(glycans.size() > 10);
}


}
}