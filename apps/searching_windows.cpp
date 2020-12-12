#include <algorithm>
#include <fstream>
#include <iostream>
#include <thread>  
#include <mutex> 
#include <chrono> 
#include <map>

#include "search_parameter.h"
#include "search_dispatcher.h"
#include "search_helper.h"

#include "../util/io/mgf_parser.h"
#include "../util/io/fasta_reader.h"
#include "../engine/protein/protein_digest.h"
#include "../engine/protein/protein_ptm.h"
#include "../engine/glycan/glycan_builder.h"
#include "../engine/search/precursor_match.h"
#include "../engine/analysis/multi_comparison.h"
#include "../engine/analysis/fdr_filter.h"
#include "../engine/analysis/score_elution.h"
#include "../engine/spectrum/lsh_clustering.h"

int main(int argc, char *argv[])
{
    // parse arguments
    extern char *optarg;
    int opt;
    std::string spectra_path = "/home/ruiz/Documents/Glycoseq-Cpp/data/ZC_20171218_C16_R1.mgf";
    std::string fasta_path = "/home/ruiz/Documents/GlycoCrushSeq/data/haptoglobin.fasta";
    std::string out_path = "result.csv";
    std::string decoy_path = "/home/ruiz/Documents/GlycoCrushSeq/data/titin.fasta";
    bool decoy_set = false;
    SearchParameter parameter;
    std::string protease = "TG";
    std::string glycan_type = "CHM";

    // pharser parameter
    while ((opt = getopt(argc, argv, ":i:f:d:o:g:k:l:m:n:e:r:s:u:w:x:y:p:h")) != EOF)
        switch(opt)
        {
            case 'i': 
                spectra_path = optarg;
                std::cout <<"The spectrum file located at " << spectra_path << std::endl; 
                break;
            case 'f':
                fasta_path = optarg;
                std::cout <<"The fasta file located at " << fasta_path << std::endl;
                break;
            case 'd': 
                decoy_set = true;
                decoy_path = optarg;
                std::cout <<"The fasta decoy file the file at " << decoy_path << std::endl; 
                break;
            case 'o': 
                out_path = optarg;
                std::cout <<"Output the file at " << out_path << std::endl; 
                break;

            case 'g':
                glycan_type = optarg;
                break;
            
            case 'k':
                parameter.ms1_by = atoi(optarg) == 0 ?
                    model::spectrum::ToleranceBy::PPM :
                    model::spectrum::ToleranceBy::Dalton;;
                break;
            case 'l':
                parameter.ms2_by = atoi(optarg)== 0 ?
                    model::spectrum::ToleranceBy::PPM :
                    model::spectrum::ToleranceBy::Dalton;
                break;
            case 'm':
                parameter.ms1_tol = atof(optarg);
                break;
            case 'n':
                parameter.ms2_tol = atof(optarg);
                break;

            case 'e':
                protease = optarg;
                break;

            case 'r':
                parameter.fdr_rate = atof(optarg);
                break;
            
            case 's':
                parameter.miss_cleavage = atoi(optarg);
                break;

            case 'u':
                parameter.neuAc_upper_bound = atoi(optarg);
                break;

            case 'w':
                parameter.neuGc_upper_bound = atoi(optarg);
                break;

            case 'x':
                parameter.hexNAc_upper_bound = atoi(optarg);
                break;

            case 'y':
                parameter.hex_upper_bound = atoi(optarg);
                break;

            case 'p': 
                parameter.n_thread =  atoi(optarg);
                std::cout <<"default thread number " << parameter.n_thread << std::endl; 
                break;
            case 'h': 
                std::cout <<"clustering mgf spectrum, parameter:\n" 
                    << "-f [file] -o [output] \n"
                    << "-t [tolerance] -l [lower_bound] -u [upper_bound] \n"
                    << "-K [hash_function_num] -L [iterations] -s [cosine] \n" 
                    << "-p [thread_num] " << std::endl; 
                return 1;
            case ':':  
                std::cout << "option needs a value" << std::endl;  
                return 1;
            case '?':  
                std::cout << "unknown option " << std::endl; 
                return 1;
            default: 
                std::cout <<"-h for help" << std::endl;
                return 1;
        }

    for(const char& c : protease)
    {
        switch (c)
        {
        case 'T': case 't':
            parameter.proteases.push_back(engine::protein::Proteases::Trypsin);
            break;

        case 'G': case 'g':
            parameter.proteases.push_back(engine::protein::Proteases::GluC);
            break;

        case 'P': case 'p':
            parameter.proteases.push_back(engine::protein::Proteases::Pepsin);
            break;
        case 'C': case 'c':
            parameter.proteases.push_back(engine::protein::Proteases::Chymotrypsin);
            break;

        default:
            break;
        }
    }

    for(const char& c : glycan_type)
    {
        switch (c)
        {
        case 'C': case 'c':
            parameter.complex = true;
            break;

        case 'H': case 'h':
            parameter.hybrid = true;
            break;

        case 'M': case 'm':
            parameter.highmannose = true;
            break;

        default:
            break;
        }
    }

    // read spectrum
    std::unique_ptr<util::io::SpectrumParser> parser = std::make_unique<util::io::MGFParser>();
    util::io::SpectrumReader spectrum_reader(spectra_path, std::move(parser));

    // read fasta and build peptides
    std::vector<std::string> peptides, decoy_peptides;
    std::vector<model::protein::Protein> proteins = ReadProteins(fasta_path);
    std::unordered_set<std::string> seqs = PeptidesDigestion(proteins, parameter);
    peptides.insert(peptides.end(), seqs.begin(), seqs.end());
    if (decoy_set)
    {
        std::vector<model::protein::Protein> decoy_proteins = ReadProteins(decoy_path);
        std::unordered_set<std::string> decoy_seqs = PeptidesDigestion(decoy_proteins, parameter);
        decoy_peptides.insert(decoy_peptides.end(), decoy_seqs.begin(), decoy_seqs.end());
    }
    else
    {
        for(auto& p: proteins)
        {
            std::string seq = p.Sequence();
            std::reverse(seq.begin(), seq.end());
            p.set_sequence(seq);
        }
        std::unordered_set<std::string> decoy_seqs = PeptidesDigestion(proteins, parameter);
        decoy_peptides.insert(decoy_peptides.end(), decoy_seqs.begin(), decoy_seqs.end());
    }
   
    // // build glycans
    std::unique_ptr<engine::glycan::GlycanBuilder> builder =
        std::make_unique<engine::glycan::GlycanBuilder>(parameter.hexNAc_upper_bound, 
            parameter.hex_upper_bound, parameter.fuc_upper_bound, 
                parameter.neuAc_upper_bound, parameter.neuGc_upper_bound,
                    parameter.complex, parameter.hybrid, parameter.highmannose);
    builder->Build();

    // search
    std::cout << "Start to scan\n"; 
    auto start = std::chrono::high_resolution_clock::now();

    // engine::spectrum::LSHClustering cluster(model::spectrum::ToleranceBy::Dalton, 0.01);
    // std::vector<engine::spectrum::EmbededSpectrum> embeds = cluster.Embed(spectrum_reader.GetSpectrum());
    // std::unordered_map<int, std::vector<engine::spectrum::EmbededSpectrum>> hash_table =  cluster.Hashing(embeds);
    // std::vector<model::spectrum::Spectrum> spectra = cluster.filter(hash_table);

    // std::cout << spectra.size() << std::endl;

    // seraching targets 
    SearchDispatcher target_searcher(spectrum_reader.GetSpectrum(), builder.get(), peptides, parameter);
    std::vector<engine::analysis::SearchResult> targets = target_searcher.Dispatch();

    // seraching decoys
    SearchDispatcher decoy_searcher(spectrum_reader.GetSpectrum(), builder.get(), decoy_peptides, parameter);
    std::vector<engine::analysis::SearchResult> decoys = decoy_searcher.Dispatch();

    std::cout << "Total target:" << targets.size() <<" decoy:" << decoys.size() << std::endl;

    // compute p value
    engine::analysis::FDRFilter tester(parameter.fdr_rate);
    tester.set_data(targets, decoys);
    tester.Init();
    std::vector<engine::analysis::SearchResult> results = tester.Filter();
    // engine::analysis::MultiComparison tester(parameter.fdr_rate);
    // std::vector<engine::analysis::SearchResult> results = tester.Tests(targets, decoys);


    // output analysis results
    ReportResults(out_path, ConvertComposition(results, builder->GlycanMapsRef()));
    // ReportResults(out_path, results);

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << "Total Time: " << duration.count() << std::endl; 

}