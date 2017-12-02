#ifndef XOLIK_PPARRAY_H
#define XOLIK_PPARRAY_H

#include <vector>
#include <unordered_map>
#include <PPData.h>
#include "Params.h"

struct Peptide {
    size_t raw_index;
    const PPData::Peptide* raw_peptide;
    size_t link_site;
    double mass;
    std::vector<std::pair<size_t /* site */, double /* shift */>> mods;
};

class PPArray {
public:
    PPArray(const PPData& ppdata, const Params& params) : peptide_array_(), mass_array_() {
        using std::vector;
        using std::pair;

        for (auto i = 0; i < ppdata.size(); ++i) {
            auto sites = FindLinkSites(params.xlsite, ppdata[i].sequence, ppdata[i].sequence_length,
                                       ppdata[i].offset, ppdata[i].protein->sequence_length);

            vector<vector<pair<size_t, double>>> mods_buffer(1, vector<pair<size_t, double>>());
            EnumerateMods(mods_buffer, 0, ppdata[i].sequence, ppdata[i].sequence_length, 
                          params.fix_mods, params.var_mods);

            for (auto site : sites) {
                for (auto& mods : mods_buffer) {
                    double total_shifts = 0;
                    for (auto& p : mods) {
                        total_shifts += p.second;
                    }
                    // PPData internally has +57.021464 at C, so remove it here
                    double total_mass = ppdata[i].mass + total_shifts - 57.021464;  
                    peptide_array_.push_back({ size_t(i), &ppdata[i], site, total_mass, mods });
                }
            }
        }

        std::stable_sort(peptide_array_.begin(), peptide_array_.end(),
                         [](const auto& a, const auto& b) { return a.mass < b.mass; });
        for (auto& p : peptide_array_) {
            mass_array_.push_back(p.mass);
        }
    }

    const Peptide& operator[](const size_t index) const {
        return peptide_array_[index];
    }

    const std::vector<double>& GetMassArray() const {
        return mass_array_;
    }

private:
    std::vector<Peptide> peptide_array_;
    std::vector<double> mass_array_;

    // helpers
    static std::vector<size_t> FindLinkSites(char xlsite, const char* sequence, size_t sequence_length,
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

    static void EnumerateMods(std::vector<std::vector<std::pair<size_t, double>>>& mods_out, size_t index,
                       const char* sequence, size_t sequence_length,
                       const std::unordered_map<char, double>& fix_mods,
                       const std::unordered_map<char, double>& var_mods) {
        if (index < sequence_length) {
            char current_aa = sequence[index];
            if (fix_mods.find(current_aa) != fix_mods.end()) {
                for (auto& mods : mods_out) {
                    mods.push_back(std::make_pair(index, fix_mods.at(current_aa)));
                }
            }
            else if (var_mods.find(current_aa) != var_mods.end()) {
                size_t current_mods_size = mods_out.size();
                for (int i = 0; i < current_mods_size; ++i) {
                    mods_out.push_back(mods_out[i]);
                    mods_out.back().push_back(std::make_pair(index, var_mods.at(current_aa)));
                }
            }
            return EnumerateMods(mods_out, index + 1, sequence, sequence_length, fix_mods, var_mods);
        }
    }
};

#endif // XOLIK_PPARRAY_H
