#ifndef APP_SEARCH_WORK_DISTRIBUTOR_H_
#define APP_SEARCH_WORK_DISTRIBUTOR_H_

#include <deque>
#include <thread>  
#include <mutex> 

#include "search_parameter.h"
#include "../algorithm/search/bucket_search.h"
#include "../engine/glycan/glycan_builder.h"
#include "../engine/analysis/search_result.h"
#include "../engine/search/precursor_match.h"
#include "../engine/search/search_glycan.h"
#include "../engine/search/search_sequence.h"

class SearchQueue
{
public:
    SearchQueue(const std::vector<model::spectrum::Spectrum>& spectra)
        { GenerateQueue(spectra); }

    SearchQueue(const SearchQueue& other)
    {
        queue_ = other.queue_;
    }

    void GenerateQueue(
        std::vector<model::spectrum::Spectrum> spectra)
    {
        for(const auto& it : spectra)
        {
            queue_.push_back(it);
        }
    }

    model::spectrum::Spectrum TryGetSpectrum()
    {
        model::spectrum::Spectrum spec;
        mutex_.lock();
            if (! queue_.empty())
            {
                spec = queue_.front();
                queue_.pop_front();
            }
            else
            {
                spec.set_scan(-1);
            }
        mutex_.unlock();
        return spec;
    }
    


protected:
    std::deque<model::spectrum::Spectrum> queue_;
    std::mutex mutex_; 
};



class SearchDispatcher
{
public:
    SearchDispatcher(
        const std::vector<model::spectrum::Spectrum>& spectra, 
        engine::glycan::GlycanBuilder* builder, const std::vector<std::string>& peptides, 
            SearchParameter parameter): queue_(SearchQueue(spectra)), builder_(builder), 
                peptides_(peptides), parameter_(parameter){}


    std::vector<engine::analysis::SearchResult> Dispatch()
    {
        std::vector<engine::analysis::SearchResult> results;
        std::vector< std::thread> thread_pool;
        for (int i = 0; i < parameter_.n_thread; i ++)
        {
            std::thread worker(&SearchDispatcher::SearchingWorker, this, std::ref(results));
            thread_pool.push_back(std::move(worker));
        }
        for (auto& worker : thread_pool)
        {
            worker.join();
        }
        return results;
    }

protected:
    void SearchingWorker(
        std::vector<engine::analysis::SearchResult>& results)
    {
        std::unique_ptr<algorithm::search::ISearch<std::string>> searcher =
            std::make_unique<algorithm::search::BucketSearch<std::string>>(parameter_.ms1_by, parameter_.ms1_tol);
        
        std::unique_ptr<algorithm::search::ISearch<std::string>> more_searcher =
            std::make_unique<algorithm::search::BucketSearch<std::string>>(parameter_.ms2_by, parameter_.ms2_tol);    
        std::unique_ptr<algorithm::search::ISearch<int>> extra_searcher =
            std::make_unique<algorithm::search::BucketSearch<int>>(parameter_.ms2_by, parameter_.ms2_tol);
    

        engine::search::PrecursorMatcher precursor_runner(std::move(searcher));
        precursor_runner.Init(peptides_, builder_->GlycanMapsRef());

        engine::search::SequenceSearch spectrum_sequencer(std::move(more_searcher));
        engine::search::GlycanSearch spectrum_searcher(std::move(extra_searcher), builder_->GlycanMapsRef());

        std::vector<engine::analysis::SearchResult> temp_result;
        engine::analysis::SearchAnalyzer analyzer;
        
        while (true)
        {
            model::spectrum::Spectrum spectrum = queue_.TryGetSpectrum();
            if (spectrum.Scan() < 0) break;

            // precusor
            auto results = precursor_runner.Match(spectrum.PrecursorMZ(), spectrum.PrecursorCharge());
            if (results.empty()) continue;

            // msms
            auto peptide_results = spectrum_sequencer.Search(spectrum.Peaks(), spectrum.PrecursorCharge(), results);
            if (peptide_results.empty()) continue;

            auto glycan_results = spectrum_searcher.Search(spectrum.Peaks(), spectrum.PrecursorCharge(), results);
            if (glycan_results.empty()) continue;

            auto searched = analyzer.Analyze(spectrum.Scan(), spectrum.Peaks(), peptide_results, glycan_results);
            temp_result.insert(temp_result.end(), searched.begin(), searched.end());
        }
        
        mutex_.lock();
            results.insert(results.end(), temp_result.begin(), temp_result.end());
        mutex_.unlock();
    }

    std::mutex mutex_; 
    SearchQueue queue_;
    engine::glycan::GlycanBuilder* builder_;
    std::vector<std::string> peptides_;
    SearchParameter parameter_;

};

#endif