// #ifndef ENGINE_DECOY_DECOY_BUILDER_H_
// #define ENGINE_DECOY_DECOY_BUILDER_H_

// #include <string>
// #include <vector>
// #include <queue>
// #include "../../util/mass/peptide.h"

// namespace engine {
// namespace decoy {

// class DecoyBuilder
// {
// public:
//     DecoyBuilder() = default;

//     void Build() 
//     {
//         std::queue<string> queue;
//         std::vector<std::string> results;
//         // bfs search
//         for(const auto& it : kAmino)
//         {
//             PeptideNode node;
//             node.sequence = it;
//             node.miss_cleavage = -1;
//             queue.push(node);
//         }

//         // grow beyond min_length
//         while (queue.size() > 0)
//         {
//             PeptideNode node = queue.front();
//             queue.pop();
            
//             // miss cleavage
//             if (node.miss_cleavage == miss_cleavage)
//                 continue;

//             // max allowed length
//             if (node.sequence.length() == kMaxLength) 
//                 continue;

//             //grow
//             for(const auto& it : kAmino)
//             {
//                 PeptideNode child;
//                 child.sequence = node.sequence + it;
//                 child.miss_cleavage = node.miss_cleavage;

//                 // if cleavable, check if results
//                 if (cleavable(child.sequence))
//                 {
//                     // minLength
//                     if (child.sequence.length() > min_length)
//                     {
//                         if (filter(child.sequence))
//                         {
//                             results.push_back(child.sequence);
//                         }
//                     }
                    
//                     child.miss_cleavage += 1;
//                 }
//                 queue.push(child);
//             }
//         }


//     }

// protected:
//     std::vector<std::string> peptide_; 
//     const int kMaxLength =  30;
//     const std::vector<std::string> kAmino 
//     {
//         "A", "C", "D", "E", "F",
//         "G", "H", "I", "K", "L",
//         "M", "N", "P", "Q", "R",
//         "S", "T", "V", "W", "Y"
//     };

// };


// }
// }


// #endif