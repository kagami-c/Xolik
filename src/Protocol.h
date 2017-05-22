#pragma once

#include <algorithm>
#include <functional>
#include <MzLoader.h>
#include <PPData.h>
#include "Params.h"
#include "XCorr.h"
#include "LimXL.h"
#include "Scores.h"

inline std::vector<size_t> FindLinkSites(char xlsite, const char* sequence, size_t sequence_length) {
    std::vector<size_t> sites;
    for (auto i = 0; i < sequence_length; ++i) {
        if (sequence[i] == xlsite) {
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
        auto sites = FindLinkSites(params.xlsite, ppdata[i].sequence, ppdata[i].sequence_length);
        for (auto site : sites) {
            generalized_peptides.push_back(std::make_pair(i, site));
            peptide_masses.push_back(ppdata[i].mass);
        }
    }

    MzLoader::Spectrum spectrum_buffer;
    while (loader.LoadNext(spectrum_buffer)) {  // loop each spectrum

        // preprocess experimental spectrum
        auto precursor_mass = (spectrum_buffer.precursor_mz - PROTON_MASS) * spectrum_buffer.precursor_charge;
        auto processed_peaks = Preprocess(spectrum_buffer.peaks, params.ms2_tolerance);
        double tolerance_in_da = precursor_mass * params.ms1_tolerance / 1000000;

        // calculate end_idx
        auto max_allowed_peptide_mass = precursor_mass - params.xlmass - params.min_allowed_mass + tolerance_in_da;
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
        if (params.use_LimXL_match) {
            max_match = LimXLMatch(peptide_masses, scores, precursor_mass, params.xlmass, tolerance_in_da, threshold);
        }
        else {
            max_match = NaiveMatch(peptide_masses, scores, precursor_mass, params.xlmass, tolerance_in_da, threshold);
        }
        if (std::get<0>(max_match) == scores.size() || std::get<1>(max_match) == scores.size()) {
            continue;
        }

        results.push_back({ spectrum_buffer.scan_num,  // SpecIdx
                            std::get<2>(max_match),  // Score
                            generalized_peptides[std::get<0>(max_match)].first,  // PeptIdx
                            generalized_peptides[std::get<0>(max_match)].second,  // link site
                            generalized_peptides[std::get<1>(max_match)].first,  // PeptIdx
                            generalized_peptides[std::get<1>(max_match)].second,  // link site
                            scores[std::get<0>(max_match)],  // Score
                            scores[std::get<1>(max_match)] });  // Score
    }

    return results;
}
