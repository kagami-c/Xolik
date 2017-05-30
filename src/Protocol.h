#pragma once

#include <future>
#include <algorithm>
#include <functional>
#include <MzLoader.h>
#include <PPData.h>
#include "Params.h"
#include "XCorr.h"
#include "LimXL.h"
#include "Scores.h"
#include "Evalue.h"

inline std::vector<size_t> FindLinkSites(char xlsite, const char* sequence, size_t sequence_length, 
                                         size_t offset, size_t protein_length) {  
                                         // for xlsite == K, being cross-linked should not be digested
    std::vector<size_t> sites;
    for (auto i = 0; i < sequence_length; ++i) {
        if (sequence[i] == xlsite) {
            // protect K & R being digested if it is cross-linked
            if ((xlsite == 'K' || xlsite == 'R') && i + 1 == sequence_length // last char of the sequence
                    && offset + sequence_length < protein_length) {  // and is not the end of the whole protein
                break;
            }
            sites.push_back(i);
        }
    }
    return sites;
}

double DetermineThreshold(Scores scores, int rank, double default_threshold) {
    if (scores.size() <= rank) { return default_threshold; }
    std::vector<double> local_scores(scores.size(), 0);
    for (int i = 0; i < scores.size(); ++i) {
        local_scores[i] = scores[i];
    }
    std::nth_element(local_scores.begin(), local_scores.begin() + rank - 1, local_scores.end(), std::greater<double>());
    return local_scores[rank - 1];
}

// The tasks of the caller of Search():
//     1. prepare resource (MzLoader & PPData)
//     2. prepare runtime parameters (params)
//     3. handle exceptions thrown by the function
//     4. receive the searching results returned by the function
std::vector<Record> Search(MzLoader& loader, const PPData& ppdata, const Params& params) {
    std::vector<Record> results;

    std::vector<std::pair<size_t /*peptide idx*/, size_t /*site idx*/>> generalized_peptides;  // peptide site combinations
    std::vector<double> peptide_masses;  // a separate peptide mass array, already sorted
    for (auto i = 0; i < ppdata.size(); ++i) {
        auto sites = FindLinkSites(params.xlsite, ppdata[i].sequence, ppdata[i].sequence_length, ppdata[i].offset, ppdata[i].protein->sequence_length);
        for (auto site : sites) {
            generalized_peptides.push_back(std::make_pair(i, site));
            peptide_masses.push_back(ppdata[i].mass);
        }
    }

    // multithread version
    if (params.enable_parallel) {
        // function declaration
        std::vector<Record> SearchBatch(std::vector<MzLoader::Spectrum> spectrum_list, const PPData& ppdata, const Params& params,
                                        std::vector<double> peptide_masses, std::vector<std::pair<size_t, size_t>> generalized_peptides);
        std::vector< std::vector<MzLoader::Spectrum> > task_list(params.thread, std::vector<MzLoader::Spectrum>());
        int idx = 0;
        MzLoader::Spectrum spectrum_buffer;
        while (loader.LoadNext(spectrum_buffer)) {  // loop each spectrum
            task_list[idx % params.thread].push_back(spectrum_buffer);
            ++idx;
        }
        std::vector< std::future< std::vector<Record> > > futures;
        for (int i = 0; i < params.thread; ++i) {
            std::future< std::vector<Record> > f = std::async(std::launch::async, SearchBatch,
                                                              task_list[i], ppdata, params, peptide_masses, generalized_peptides);
            futures.push_back(std::move(f));
        }
        std::vector<Record> records;
        for (int i = 0; i < params.thread; ++i) {
            std::vector<Record> r = futures[i].get();
            records.insert(records.end(), r.begin(), r.end());
        }
        return records;
    }

    MzLoader::Spectrum spectrum_buffer;
    while (loader.LoadNext(spectrum_buffer)) {  // loop each spectrum

        // preprocess experimental spectrum
        auto precursor_mass = (spectrum_buffer.precursor_mz - PROTON_MASS) * spectrum_buffer.precursor_charge;
        auto processed_peaks = Preprocess(spectrum_buffer.peaks, params.ms2_tolerance);
        //double tolerance_in_da = precursor_mass * params.ms1_tolerance / 1000000;
        double left_tol = precursor_mass - precursor_mass / (1 + params.ms1_tolerance / 1000000);
        double right_tol = precursor_mass / (1 - params.ms1_tolerance / 1000000) - precursor_mass;

        // calculate end_idx
        auto max_allowed_peptide_mass = precursor_mass - params.xlmass - params.min_allowed_mass + right_tol;
        auto last_iter = std::upper_bound(peptide_masses.begin(), peptide_masses.end(), max_allowed_peptide_mass);
        auto end_idx = std::distance(peptide_masses.begin(), last_iter);

        // build scores
        int maximum_charge = spectrum_buffer.precursor_charge > 1 ? spectrum_buffer.precursor_charge - 1 : 1;
        maximum_charge = maximum_charge > 6 ? 6 : maximum_charge;
        Scores scores(processed_peaks, precursor_mass, ppdata, generalized_peptides, params.ms2_tolerance, end_idx,
                      maximum_charge);

        // determine threshold
        double threshold = params.threshold;
        if (params.enable_rank) {
            // if scores.size() <= rank, fallback to default_threshold params.threshold
            threshold = DetermineThreshold(scores, params.rank, params.threshold);
            if (threshold < params.threshold) {
                threshold = params.threshold;
            }
        }

        // match algorithm
        std::tuple<CandIdx, CandIdx, Score> max_match;
        std::vector<double> collected_scores;  // for evalue estimation
        if (params.use_LimXL_match) {
            max_match = LimXLMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold, left_tol, right_tol);
        }
        else {
            max_match = NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold, left_tol, right_tol, 
                                   collected_scores, false, 0);
        }
        if (std::get<0>(max_match) == scores.size() || std::get<1>(max_match) == scores.size()) {
            continue;
        }
        double report_score = std::get<2>(max_match);  // raw Xcorr score

        // estimate evalue
        if (params.use_E_value) {
            // GARBAGE, I really don't understand why PEOPLE will be proud of features using such an STUPID solution
            // Both solution and the value towards this solution are GARBAGE.
            double additional_tol = 0.0;
            NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold,
                       left_tol, right_tol, collected_scores, true, params.histogram_size);
            while (collected_scores.size() < params.histogram_size && additional_tol + std::max(left_tol, right_tol) <= 20.0) {
                additional_tol += 1.0;
                NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold,
                           additional_tol + left_tol, -(additional_tol - 1.0 + left_tol), collected_scores, true, params.histogram_size);
                if (collected_scores.size() >= params.histogram_size) {
                    break;
                }
                NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold,
                           -(right_tol + additional_tol - 1.0), right_tol + additional_tol, collected_scores, true, params.histogram_size);
            }

            double evalue = CalculateEValue(std::get<2>(max_match), collected_scores);
            report_score = -log10(evalue);  // - log 10 evalue to make it compatible with FDR control
        }

        results.push_back({ spectrum_buffer.scan_num,  // SpecIdx
                            report_score,  // Score  
                            generalized_peptides[std::get<0>(max_match)].first,  // PeptIdx
                            generalized_peptides[std::get<0>(max_match)].second,  // link site
                            generalized_peptides[std::get<1>(max_match)].first,  // PeptIdx
                            generalized_peptides[std::get<1>(max_match)].second,  // link site
                            scores[std::get<0>(max_match)],  // Score
                            scores[std::get<1>(max_match)] });  // Score
    }

    return results;
}

// for multithread computing, need a refactor
std::vector<Record> SearchBatch(std::vector<MzLoader::Spectrum> spectrum_list, const PPData& ppdata, const Params& params,
                                std::vector<double> peptide_masses, std::vector<std::pair<size_t, size_t>> generalized_peptides) {
    std::vector<Record> records;
    for (int i = 0; i < spectrum_list.size(); ++i) {
        MzLoader::Spectrum spectrum_buffer = spectrum_list[i];

        // preprocess experimental spectrum
        auto precursor_mass = (spectrum_buffer.precursor_mz - PROTON_MASS) * spectrum_buffer.precursor_charge;
        auto processed_peaks = Preprocess(spectrum_buffer.peaks, params.ms2_tolerance);
        //double tolerance_in_da = precursor_mass * params.ms1_tolerance / 1000000;
        double left_tol = precursor_mass - precursor_mass / (1 + params.ms1_tolerance / 1000000);
        double right_tol = precursor_mass / (1 - params.ms1_tolerance / 1000000) - precursor_mass;

        // calculate end_idx
        auto max_allowed_peptide_mass = precursor_mass - params.xlmass - params.min_allowed_mass + right_tol;
        auto last_iter = std::upper_bound(peptide_masses.begin(), peptide_masses.end(), max_allowed_peptide_mass);
        auto end_idx = std::distance(peptide_masses.begin(), last_iter);

        // build scores
        int maximum_charge = spectrum_buffer.precursor_charge > 1 ? spectrum_buffer.precursor_charge - 1 : 1;
        maximum_charge = maximum_charge > 6 ? 6 : maximum_charge;
        Scores scores(processed_peaks, precursor_mass, ppdata, generalized_peptides, params.ms2_tolerance, end_idx,
                      maximum_charge);

        // determine threshold
        double threshold = params.threshold;
        if (params.enable_rank) {
            // if scores.size() <= rank, fallback to default_threshold params.threshold
            threshold = DetermineThreshold(scores, params.rank, params.threshold);
            if (threshold < params.threshold) {
                threshold = params.threshold;
            }
        }

        // match algorithm
        std::tuple<CandIdx, CandIdx, Score> max_match;
        std::vector<double> collected_scores;  // for evalue estimation
        if (params.use_LimXL_match) {
            max_match = LimXLMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold, left_tol, right_tol);
        }
        else {
            max_match = NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold, left_tol, right_tol,
                                   collected_scores, false, 0);
        }
        if (std::get<0>(max_match) == scores.size() || std::get<1>(max_match) == scores.size()) {
            continue;
        }
        double report_score = std::get<2>(max_match);  // raw Xcorr score

                                                       // estimate evalue
        if (params.use_E_value) {
            // GARBAGE, I really don't understand why PEOPLE will be proud of features using such an STUPID solution
            // Both solution and the value towards this solution are GARBAGE.
            double additional_tol = 0.0;
            NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold,
                       left_tol, right_tol, collected_scores, true, params.histogram_size);
            while (collected_scores.size() < params.histogram_size && additional_tol + std::max(left_tol, right_tol) <= 20.0) {
                additional_tol += 1.0;
                NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold,
                           additional_tol + left_tol, -(additional_tol - 1.0 + left_tol), collected_scores, true, params.histogram_size);
                if (collected_scores.size() >= params.histogram_size) {
                    break;
                }
                NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, threshold,
                           -(right_tol + additional_tol - 1.0), right_tol + additional_tol, collected_scores, true, params.histogram_size);
            }

            double evalue = CalculateEValue(std::get<2>(max_match), collected_scores);
            report_score = -log10(evalue);  // - log 10 evalue to make it compatible with FDR control
        }

        records.push_back({ spectrum_buffer.scan_num,  // SpecIdx
                          report_score,  // Score  
                          generalized_peptides[std::get<0>(max_match)].first,  // PeptIdx
                          generalized_peptides[std::get<0>(max_match)].second,  // link site
                          generalized_peptides[std::get<1>(max_match)].first,  // PeptIdx
                          generalized_peptides[std::get<1>(max_match)].second,  // link site
                          scores[std::get<0>(max_match)],  // Score
                          scores[std::get<1>(max_match)] });  // Score
    }
    return records;
}