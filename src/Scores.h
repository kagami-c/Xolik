#pragma once

#include <PPData.h>
#include <vector>
#include "XCorr.h"

class Scores {
public:
    Scores(const std::vector<double>& processed_peaks, double precursor_mass,
           const PPData& ppdata, const std::vector<std::pair<size_t, size_t>>& generalized_peptides,
           double ms2_tolerance, size_t total_size, int maximum_charge)
            : precursor_mass_(precursor_mass), processed_peaks_(processed_peaks), ppdata_(ppdata),
              generalized_peptides_(generalized_peptides), ms2_tolerance_(ms2_tolerance), maximum_charge_(maximum_charge),
              size_(total_size), cache_flags_(total_size, 0), cache_scores_(total_size, 0.0) {}

    double operator[](int idx) {
        if (cache_flags_[idx] == 1) {
            return cache_scores_[idx];
        }
        else {  // cache_flags == 0
            const auto& peptide = ppdata_[generalized_peptides_[idx].first];
            auto site = generalized_peptides_[idx].second;
            auto mass_shift = precursor_mass_ - peptide.mass;
            auto xcorr = XCorr(processed_peaks_, ms2_tolerance_, peptide.sequence,
                               peptide.sequence_length, site, mass_shift, maximum_charge_);
            cache_scores_[idx] = xcorr;
            cache_flags_[idx] = 1;
            return xcorr;
        }
    }

    size_t size() const {
        return size_;
    }

private:
    const double precursor_mass_;
    const std::vector<double>& processed_peaks_;
    const PPData& ppdata_;
    const std::vector<std::pair<size_t, size_t>>& generalized_peptides_;
    const double ms2_tolerance_;
    const int maximum_charge_;

    size_t size_;
    std::vector<int> cache_flags_;
    std::vector<double> cache_scores_;
};
