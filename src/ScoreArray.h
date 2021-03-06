#ifndef XOLIK_SCOREARRAY_H
#define XOLIK_SCOREARRAY_H

#include <vector>
#include "XCorr.h"
#include "PPArray.h"

class ScoreArray {
public:
    ScoreArray(const std::vector<double>& processed_peaks, double precursor_mass, 
               const PPArray& pparray, size_t total_size,
               double ms2_tolerance, int maximum_charge)
             : precursor_mass_(precursor_mass), processed_peaks_(processed_peaks), 
               ms2_tolerance_(ms2_tolerance), maximum_charge_(maximum_charge), 
               pparray_(pparray), size_(total_size), 
               cache_flags_(total_size, 0), cache_scores_(total_size, 0.0)  {}

    double operator[](size_t idx) {
        if (cache_flags_[idx] == 1) {
            return cache_scores_[idx];
        } else {  // cache_flags == 0
            const auto& peptide = *(pparray_[idx].raw_peptide);
            auto site = pparray_[idx].link_site;
            auto mass_shift = precursor_mass_ - peptide.mass;
//            auto xcorr = XCorr(processed_peaks_, ms2_tolerance_, peptide.sequence,
//                               peptide.sequence_length, site, mass_shift, maximum_charge_);
            const auto& mods = pparray_[idx].mods;
            auto xcorr = ModXCorr(processed_peaks_, ms2_tolerance_, 
                                  peptide.sequence, peptide.sequence_length, 
                                  site, mass_shift, maximum_charge_, mods);
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
    const double ms2_tolerance_;
    const int maximum_charge_;

    const PPArray& pparray_;
    size_t size_;
    std::vector<int> cache_flags_;
    std::vector<double> cache_scores_;
};

#endif // XOLIK_SCOREARRAY_H
