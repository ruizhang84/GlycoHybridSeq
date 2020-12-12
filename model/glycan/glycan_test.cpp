#define BOOST_TEST_MODULE GlycanTest
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "nglycan_complex.h"
#include "nglycan_hybrid.h"
#include "highmannose.h"

using namespace std;
using namespace model::glycan;

// BOOST_AUTO_TEST_CASE( Glycan_test ) 
// {
//     std::unique_ptr<Glycan> nglycan = std::make_unique<NGlycanComplex>();
//     nglycan->set_table(0, 2); // core 
//     nglycan->set_table(1, 3); // 
//     nglycan->set_table(2, 1); // fuc
//     nglycan->set_table(3, 1); // bisect
//     nglycan->set_table(4, 1); // branch 
//     nglycan->set_table(5, 1);
//     nglycan->set_table(8, 1);
//     nglycan->set_table(9, 1);  
//     nglycan->Composition()[Monosaccharide::GlcNAc] = 3;
//     nglycan->Composition()[Monosaccharide::Gal] = 0;
//     nglycan->Composition()[Monosaccharide::Fuc] = 1;
//     nglycan->Composition()[Monosaccharide::Man] = 3;
//     nglycan->Composition()[Monosaccharide::NeuAc] = 0;

//     std::cout << nglycan->Name() << std::endl;
//     std::cout << nglycan->ID() << std::endl;

//     std::vector<std::unique_ptr<Glycan>> glycans = nglycan->Grow(Monosaccharide::GlcNAc);
//     BOOST_CHECK(glycans.size() == 2);
//     for (int i = 0; i < (int) glycans.size(); i++){
//         std::cout << glycans[i]->ID() << std::endl;
//         std::cout << glycans[i]->Name() << std::endl;
//     }

//     NGlycanComplex nglycan_1;
//     std::map<Monosaccharide, int> composite;
//     composite[Monosaccharide::GlcNAc] = 12;
//     composite[Monosaccharide::Gal] = 12;
//     composite[Monosaccharide::Fuc] = 1;
//     composite[Monosaccharide::Man] = 3;
//     composite[Monosaccharide::NeuAc] = 6;

//     nglycan_1.set_composition(composite);
//     std::string compos = nglycan_1.Name();
//     std::cout << compos << std::endl;
// }

// BOOST_AUTO_TEST_CASE( Glycan_test ) 
// {
//     std::unique_ptr<Glycan> nglycan = std::make_unique<NGlycanHybrid>();
//     nglycan->set_table(0, 2); // core 
//     nglycan->set_table(1, 3); // 
//     nglycan->set_table(2, 1); // fuc
//     nglycan->set_table(4, 1); // man branch 
//     nglycan->set_table(6, 1); // GlcNAc
//     nglycan->set_table(7, 1);  
//     nglycan->Composition()[Monosaccharide::GlcNAc] = 4;
//     nglycan->Composition()[Monosaccharide::Gal] = 0;
//     nglycan->Composition()[Monosaccharide::Fuc] = 1;
//     nglycan->Composition()[Monosaccharide::Man] = 4;
//     nglycan->Composition()[Monosaccharide::NeuAc] = 0;

//     // std::cout << nglycan->Name() << std::endl;
//     // std::cout << nglycan->ID() << std::endl;

//     std::vector<std::unique_ptr<Glycan>> glycans = nglycan->Grow(Monosaccharide::Gal);
//     // BOOST_CHECK(glycans.size() == 2);
//     for (int i = 0; i < (int) glycans.size(); i++){
//         std::cout << glycans[i]->ID() << std::endl;
//         std::cout << glycans[i]->Name() << std::endl;
//     }
// }


BOOST_AUTO_TEST_CASE( Glycan_test ) 
{
    std::unique_ptr<Glycan> nglycan = std::make_unique<HighMannose>();
    nglycan->set_table(0, 2); // core 
    nglycan->set_table(1, 3); // 
    nglycan->set_table(2, 1); // fuc
    nglycan->set_table(3, 1); // man branch 
    nglycan->set_table(4, 1); // 
    nglycan->set_table(5, 1);  
    nglycan->Composition()[Monosaccharide::GlcNAc] = 2;
    nglycan->Composition()[Monosaccharide::Gal] = 0;
    nglycan->Composition()[Monosaccharide::Fuc] = 1;
    nglycan->Composition()[Monosaccharide::Man] = 6;
    nglycan->Composition()[Monosaccharide::NeuAc] = 0;

    // std::cout << nglycan->Name() << std::endl;
    // std::cout << nglycan->ID() << std::endl;

    std::vector<std::unique_ptr<Glycan>> glycans = nglycan->Grow(Monosaccharide::Man);
    // BOOST_CHECK(glycans.size() == 2);
    for (int i = 0; i < (int) glycans.size(); i++){
        std::cout << glycans[i]->ID() << std::endl;
        std::cout << glycans[i]->Name() << std::endl;
    }
}


// BOOST_AUTO_TEST_CASE( glycan_row_test ) 
// {
//     NGlycanComplex nglycan;
//     std::vector<std::unique_ptr<Glycan>> glycans = nglycan.Grow(Monosaccharide::GlcNAc);
//     std::vector<std::unique_ptr<Glycan>> glycans_2 = glycans.front()->Grow(Monosaccharide::GlcNAc);
//     std::vector<std::unique_ptr<Glycan>> glycans_3 = glycans_2.front()->Grow(Monosaccharide::Man);
//     std::string name = glycans_3.front()->Name();
//     std::cout << name << std::endl;
//     BOOST_CHECK(glycans_3.front()->Composition()[Monosaccharide::GlcNAc] == 2);
// }








