#include <algorithm>
#include <fstream>
#include <iostream>
#include <thread>  
#include <mutex> 
#include <chrono> 
#include <map>

#include <argp.h>

#include "search_parameter.h"
#include "search_dispatcher.h"
#include "search_helper.h"

#include "../util/io/mgf_parser.h"
#include "../util/io/fasta_reader.h"
#include "../engine/protein/protein_digest.h"
#include "../engine/protein/protein_ptm.h"
#include "../engine/glycan/glycan_builder.h"
#include "../engine/search/precursor_match.h"
#include "../engine/analysis/fdr_filter.h"
#include "../engine/analysis/score_elution.h"
#include "../engine/spectrum/lsh_clustering.h"

const char *argp_program_version =
  "GlycoCrushseq v1.0";
const char *argp_program_bug_address =
  "<rz20@iu.edu>";

static char doc[] =
  "GlycoCrushseq -- a program to search glycopeptide from high thoughput LS-MS/MS";

static struct argp_option options[] = {
    {"spath", 'i',    "spectrum.mgf",  0,  "mgf, Spectrum MS/MS Input Path" },
    {"fpath", 'f',    "protein.fasta",  0,  "fasta, Protein Sequence Input Path" },
    {"dpath", 'd',    "reversed",  0,  "fasta, Protein Sequence for Decoy" },
    {"output",    'o',    "result.csv",   0,  "csv, Results Output Path" },
    {"pthread",   'p',  "6",  0,  "Number of Searching Threads" },
    {"enzyme",   'e',  "TG",  0,  "The Digestion Enzyme, Trypsin (T), Pepsin (P), Chymotrypsin (C), GluC (G)" }, 
    {"miss_cleavage",   's',  "2",  0,  "The Missing Cleavage Upto" },    
    {"HexNAc",   'x',  "12",  0,  "Search Up to Number of HexNAc" },
    {"HexNA",   'y',  "12",  0,  "Search Up to Number of Hex" },
    {"Fuc",   'z',  "5",  0,  "Search Up to Number of Fuc" },
    {"NeuAc",   'u',  "4",  0,  "Search Up to Number of NeuAc" },
    {"NeuGc",   'w',  "0",  0,  "Search Up to Number of NeuGc" },
    {"ms1_tol",   'm',  "10",  0,  "MS Tolereance" },
    {"ms2_tol",   'n',  "0.01",  0,  "MS2 Tolereance" },
    {"ms1_by",   'k',  "0",  0, "MS Tolereance By Int: PPM (0) or Dalton (1)" },
    {"ms2_by",   'l',  "1",  0, "MS2 Tolereance By Int: PPM (0) or Dalton (1)" },
    {"fdr_rate",   'r',  "0.01",  0, "FDR rate" },
    {"glycan_type",   'g',  "C",  0, "The Searching Glycan Type, Complex (C), Hybrid (H), High Mannose(M)"},
    { 0 }
};

static std::string default_spectra_path = "ZC_20171218_C22_R1.mgf";
static std::string default_fasta_path = "haptoglobin.fasta";
static std::string default_decoy_path = "titin.fasta";
static std::string default_out_path = "result.csv";
static std::string default_digestion = "TG";
static std::string default_glycan_type = "CHM";

struct arguments
{
    char * spectra_path = const_cast<char*> (default_spectra_path.c_str());
    char * fasta_path = const_cast<char*> (default_fasta_path.c_str());
    char * out_path = const_cast<char*> (default_out_path.c_str());
    // decoy
    bool decoy_set = false;
    char * decoy_path = const_cast<char*> (default_decoy_path.c_str());
    //digestion
    int miss_cleavage = 2;
    char * digestion = const_cast<char*> (default_digestion.c_str());
    // upper bound of glycan seaerch
    int n_thread = 6;
    int hexNAc_upper_bound = 12;
    int hex_upper_bound = 12;
    int fuc_upper_bound = 5;
    int neuAc_upper_bound = 4;
    int neuGc_upper_bound = 0;
    // searching precision
    double ms1_tol = 10;
    double ms2_tol = 0.01;
    int ms1_by = 0;
    int ms2_by = 1;
    // fdr
    double fdr_rate = 0.01;
    // glycan type
    char * glycan_type = const_cast<char*> (default_glycan_type.c_str());
};


static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
    error_t err = 0;
    struct arguments *arguments =  static_cast<struct arguments*>(state->input);

    switch (key)
    {
    case 'd':
        arguments->decoy_set = true;
        arguments->decoy_path = arg;
        break;

    case 'e':
        arguments->digestion = arg;
        break;

    case 'f':
        arguments->fasta_path = arg;
        break;

    case 'g':
        arguments->glycan_type = arg;
        break;

    case 'i':
        arguments->spectra_path = arg;
        break;

    case 'k':
        arguments->ms1_by = atoi(arg);
        break;

    case 'l':
        arguments->ms2_by = atoi(arg);
        break;

    case 'm':
        arguments->ms1_tol = atof(arg);
        break;

    case 'n':
        arguments->ms2_tol = atof(arg);
        break;

    case 'o':
        arguments->out_path = arg;
        break;
    
    case 'p':
        arguments->n_thread = atoi(arg);
        break;
    
    case 'r':
        arguments->fdr_rate = atof(arg);
        break;
    
    case 's':
        arguments->miss_cleavage = atoi(arg);
        break;

    case 'u':
        arguments->neuAc_upper_bound = atoi(arg);
        break;

    case 'w':
        arguments->neuGc_upper_bound = atoi(arg);
        break;

    case 'x':
        arguments->hexNAc_upper_bound = atoi(arg);
        break;

    case 'y':
        arguments->hex_upper_bound = atoi(arg);
        break;

    case 'z':
        arguments->fuc_upper_bound = atoi(arg);
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return err;
}

static struct argp argp = { options, parse_opt, 0, doc };


SearchParameter GetParameter(const struct arguments& arguments)
{
    SearchParameter parameter;
    parameter.n_thread = arguments.n_thread;
    parameter.miss_cleavage = arguments.miss_cleavage;
    parameter.hexNAc_upper_bound = arguments.hexNAc_upper_bound;
    parameter.hex_upper_bound = arguments.hex_upper_bound;
    parameter.neuAc_upper_bound = arguments.neuAc_upper_bound;
    parameter.neuGc_upper_bound = arguments.neuGc_upper_bound;
    parameter.ms1_tol = arguments.ms1_tol;
    parameter.ms1_by = arguments.ms1_by == 0 ?
        model::spectrum::ToleranceBy::PPM :
        model::spectrum::ToleranceBy::Dalton;
    parameter.ms2_tol = arguments.ms2_tol;
    parameter.ms2_by = arguments.ms2_by == 0 ?
        model::spectrum::ToleranceBy::PPM :
        model::spectrum::ToleranceBy::Dalton;
    parameter.fdr_rate = arguments.fdr_rate;
    std::string protease(arguments.digestion);
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
    std::string glycan_type(arguments.glycan_type);
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

    return parameter;
}

int main(int argc, char *argv[])
{
    // parse arguments
    struct arguments arguments;
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    std::string spectra_path(arguments.spectra_path) ;
    std::string fasta_path(arguments.fasta_path);
    std::string decoy_path(arguments.decoy_path); 
    std::string out_path(arguments.out_path);
    SearchParameter parameter = GetParameter(arguments);

    // read spectrum
    std::unique_ptr<util::io::SpectrumParser> parser = std::make_unique<util::io::MGFParser>();
    util::io::SpectrumReader spectrum_reader(spectra_path, std::move(parser));

    // read fasta and build peptides
    std::vector<std::string> peptides, decoy_peptides;
    std::vector<model::protein::Protein> proteins = ReadProteins(fasta_path);
    std::unordered_set<std::string> seqs = PeptidesDigestion(proteins, parameter);
    peptides.insert(peptides.end(), seqs.begin(), seqs.end());
    if (arguments.decoy_set)
    {
        std::vector<model::protein::Protein> decoy_proteins = ReadProteins(decoy_path);
        std::unordered_set<std::string> decoy_seqs = PeptidesDigestion(decoy_proteins, parameter);
        decoy_peptides.insert(decoy_peptides.end(), decoy_seqs.begin(), decoy_seqs.end());
    }
    else
    {
        // for(const auto& s : peptides)
        // {
        //     decoy_peptides.push_back(engine::protein::ReverseNGlycopeptide(s));
        // }
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
