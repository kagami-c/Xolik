// Copyright (C) 2016

#pragma once

#include <vector>
#include <algorithm>
#include <cstring>
#include <PPData.h>
#include "Params.h"

inline bool IsDecoy(const char* protein_name) {
    return 0 == strncmp("DECOY_", protein_name, 6);
}

std::vector<double> CalculateFDR(std::vector<Record>& records, const PPData& ppdata) {
    std::stable_sort(records.begin(), records.end(), 
                     [](const auto& one, const auto& another) {
                         return one.score > another.score;
                     });
    auto target_count = 0;
    auto semi_decoy_count = 0;
    auto decoy_count = 0;
    std::vector<double> fdr(records.size());
    for (auto i = 0; i < records.size(); ++i) {
        const auto& record = records[i];
        auto alpha_decoy = IsDecoy(ppdata[record.alpha_idx].protein->name);
        auto beta_decoy = IsDecoy(ppdata[record.beta_idx].protein->name);
        if (alpha_decoy && beta_decoy) {
            decoy_count += 1;
        }
        else if (alpha_decoy || beta_decoy) {
            semi_decoy_count += 1;
        }
        else {
            target_count += 1;
        }
        if (target_count == 0) {
            if (semi_decoy_count == 0 && decoy_count == 0) {
                throw std::runtime_error("Impossible path");
            }
            fdr[i] = 1;
        }
        else {
            if (semi_decoy_count <= decoy_count) {
                fdr[i] = double(decoy_count) / double(target_count);
            }
            else {
                fdr[i] = double(semi_decoy_count - decoy_count) / double(target_count);
            }
        }
    }

    // calculate q value, backtracing
    std::vector<double> q_values(records.size());
    double current_min = 1;
    for (int i = records.size() - 1; i >= 0; --i) {
        auto current_fdr = fdr[i];
        if (current_fdr < current_min) {
            current_min = current_fdr;
        }
        q_values[i] = current_min;
    }
    return q_values;
}
