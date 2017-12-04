#ifndef XOLIK_SEARCH_H
#define XOLIK_SEARCH_H

#include <future>
#include <algorithm>
#include <functional>
#include <MzLoader.h>
#include <PPData.h>
#include "Params.h"
#include "XCorr.h"
#include "Match.h"
#include "ScoreArray.h"
#include "Evalue.h"
#include "PPArray.h"

double DetermineThreshold(ScoreArray& scores, int rank, double default_threshold) {
    if (scores.size() <= rank) {
        return default_threshold;
    }
    std::vector<double> local_scores(scores.size(), 0);
    for (int i = 0; i < scores.size(); ++i) {
        local_scores[i] = scores[i];
    }
    std::nth_element(local_scores.begin(), local_scores.begin() + rank - 1, local_scores.end(), std::greater<double>());
    return local_scores[rank - 1];
}

// subroutine for searching one spectrum
bool SearchSpectrum(const MzLoader::Spectrum& spectrum, const std::vector<double>& mass_array,
                    const PPArray& pparray, const Params& params, Record& record_out) {

    // preprocess experimental spectrum
    double precursor_mass = (spectrum.precursor_mz - PROTON_MASS) * spectrum.precursor_charge;
    std::vector<double> processed_peaks = Preprocess(spectrum.peaks, params.ms2_tolerance);
    double left_tol = precursor_mass - precursor_mass / (1 + params.ms1_tolerance / 1000000);
    double right_tol = precursor_mass / (1 - params.ms1_tolerance / 1000000) - precursor_mass;

    // calculate end_idx
    double max_allowed_peptide_mass = precursor_mass - params.xlmass - params.min_allowed_mass + right_tol;
    auto last_iter = std::upper_bound(mass_array.begin(), mass_array.end(), max_allowed_peptide_mass);
    size_t end_idx = std::distance(mass_array.begin(), last_iter);

    // build score array
    int maximum_charge = spectrum.precursor_charge > 1 ? spectrum.precursor_charge - 1 : 1;
    maximum_charge = maximum_charge > 6 ? 6 : maximum_charge;
    ScoreArray score_array(processed_peaks, precursor_mass, params.ms2_tolerance, end_idx, maximum_charge, pparray);

    // determine threshold
    double threshold = params.threshold;
    if (params.enable_rank) {
        // if scores.size() <= rank, fallback to default_threshold params.threshold
        threshold = DetermineThreshold(score_array, params.rank, params.threshold);
        if (threshold < params.threshold) {
            threshold = params.threshold;
        }
    }

    // match algorithm
    std::tuple<CandIdx, CandIdx, Score> max_match;
    std::vector<double> collected_scores; // for evalue estimation
    int match_count = 0;  // for evalue estimation
    if (params.use_LimXL_match) {
        max_match = XolikMatch(mass_array, score_array, precursor_mass, params.xlmass, threshold, left_tol, right_tol, match_count);
    } else {
        max_match = NaiveMatch(mass_array, score_array, precursor_mass, params.xlmass, threshold, left_tol, right_tol,
                               collected_scores, false, 0, match_count);
    }
    if (std::get<0>(max_match) == score_array.size() || std::get<1>(max_match) == score_array.size()) {
        return false;
    }
    double report_score = std::get<2>(max_match); // raw Xcorr score

    // estimate evalue
    if (params.use_E_value) {
        double additional_tol = 0.0;
        int collect_count = 0;
        NaiveMatch(mass_array, score_array, precursor_mass, params.xlmass, threshold,
                   left_tol, right_tol, collected_scores, true, params.histogram_size, collect_count);
        while (collected_scores.size() < params.histogram_size && additional_tol + std::max(left_tol, right_tol) <= 20.0) {
            additional_tol += 1.0;
            NaiveMatch(mass_array, score_array, precursor_mass, params.xlmass, threshold,
                       additional_tol + left_tol, -(additional_tol - 1.0 + left_tol), collected_scores, true, params.histogram_size, collect_count);
            if (collected_scores.size() >= params.histogram_size) {
                break;
            }
            NaiveMatch(mass_array, score_array, precursor_mass, params.xlmass, threshold,
                       -(right_tol + additional_tol - 1.0), right_tol + additional_tol, collected_scores, true, params.histogram_size, collect_count);
        }
        collect_count = collected_scores.size();
        double evalue = CalculateEValue(std::get<2>(max_match), collected_scores);
        report_score = -evalue -log10(double(match_count) / double(collect_count));
    }

    record_out = {
        spectrum.scan_num,  // SpecIdx
        report_score,  // Score  
        pparray[std::get<0>(max_match)].raw_index,  // PeptIdx
        pparray[std::get<0>(max_match)].link_site,  // link site
        pparray[std::get<1>(max_match)].raw_index,  // PeptIdx
        pparray[std::get<1>(max_match)].link_site,  // link site
        score_array[std::get<0>(max_match)],  // Score
        score_array[std::get<1>(max_match)]  // Score
    };
    return true;
}

// for multithread computing
std::vector<Record> SearchBatch(const std::vector<MzLoader::Spectrum>& spectrum_list, 
                                const std::vector<double>& mass_array,
                                const PPArray& pparray, 
                                const Params& params) {
    std::vector<Record> records;
    Record record_buffer;
    for (int i = 0; i < spectrum_list.size(); ++i) {
        const MzLoader::Spectrum& spectrum_buffer = spectrum_list[i];
        bool success = SearchSpectrum(spectrum_buffer, mass_array, pparray, params, record_buffer);
        if (success) {
            records.push_back(record_buffer);
        }
    }
    return records;
}

// The tasks of the caller of Search():
//     1. prepare resource (MzLoader & PPData)
//     2. prepare runtime parameters (params)
//     3. handle exceptions thrown by the function
//     4. receive the searching results returned by the function
std::vector<Record> Search(MzLoader& loader, const PPData& ppdata, const Params& params) {
    std::vector<Record> records;

    PPArray pparray(ppdata, params);
    std::vector<double> peptide_masses = pparray.GetMassArray();
    
    if (!params.enable_parallel) {
        MzLoader::Spectrum spectrum_buffer;
        Record record_buffer;
        while (loader.LoadNext(spectrum_buffer)) { // loop each spectrum
            bool success = SearchSpectrum(spectrum_buffer, peptide_masses, pparray, params, record_buffer);
            if (success) {
                records.push_back(record_buffer);
            }
        }
    } else {  // multithread version
        std::vector<std::vector<MzLoader::Spectrum>> task_list(params.thread, std::vector<MzLoader::Spectrum>());
        int idx = 0;
        MzLoader::Spectrum spectrum_buffer;
        while (loader.LoadNext(spectrum_buffer)) {  // loop each spectrum
            task_list[idx % params.thread].push_back(spectrum_buffer);
            ++idx;
        }
        std::vector<std::future<std::vector<Record>>> futures;
        for (int i = 0; i < params.thread; ++i) {
            std::future<std::vector<Record>> f = std::async(std::launch::async, 
                    SearchBatch, task_list[i], peptide_masses, pparray, params);
            futures.push_back(std::move(f));
        }
        std::vector<std::vector<Record>> results;  // TODO: better multithread code, this version is quite ugly
        for (int i = 0; i < params.thread; ++i) {
            std::vector<Record> r = futures[i].get();
            results.push_back(r);
//            records.insert(records.end(), r.begin(), r.end());
        }
        std::vector<int> indexes(params.thread, 0);
        bool flag = true;
        while (flag) {
            flag = false;
            for (int i = 0; i < params.thread; ++i) {
                if (indexes[i] < results[i].size()) {
                    flag = true;
                    records.push_back(results[i][indexes[i]++]);
                }
            }
        }
    }

    return records;
}

#endif // XOLIK_SEARCH_H
