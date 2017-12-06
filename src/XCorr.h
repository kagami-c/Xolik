#ifndef XOLIK_XCORR_H
#define XOLIK_XCORR_H

// novel preprocessing and scoring
std::vector<double> Preprocess(const std::vector<std::pair<double, double>>& peaks, 
                               double resolution) {
    if (peaks.empty()) {
        return std::vector<double>();
    }
    int vector_size = static_cast<int>(peaks[peaks.size() - 1].first / resolution) + 1;
    std::vector<double> processed_peaks(vector_size);

    // normalize in local range
    double min_mass = peaks[0].first;
    double max_mass = peaks[peaks.size() - 1].first;
    const int num_step = 10;
    int step_unit = vector_size / 10 + 1;

    int start = 0;
    int end = 0;  // exclude
    int step = 0;
    do {
        double lower_bound = step * step_unit * resolution;
        double upper_bound = lower_bound + step_unit * resolution;
        start = end;
        while (end < peaks.size() && peaks[end].first <= upper_bound) {
            ++end;
        }
        // from start to end (excluded) is the range for normalization
        if (start == end) {
            continue;
        }
        double local_max_intensity = 0;
        for (int i = start; i < end; ++i) {
            if (peaks[i].second > local_max_intensity) {
                local_max_intensity = peaks[i].second;
            }
        }

        // fill processed peaks
        double normalized_factor = 50 / sqrt(local_max_intensity);
        for (auto i = start; i < end; ++i) {
            double normalized_intensity = sqrt(peaks[i].second) * normalized_factor;
            int index = int(peaks[i].first / resolution);
            if (normalized_intensity > processed_peaks[index]) {
                processed_peaks[index] = normalized_intensity;
            }
        }
    } while ((++step) < num_step);

    // mean filter preprocess, always int +- 75 bins range
    const int offset = 75;
    const int window_size = 2 * offset + 1;
    if (processed_peaks.size() <= window_size) {  // ignore small spectrum
        return processed_peaks;
    }
    std::vector<double> filtered(processed_peaks.size());

    double sum = 0;
    for (size_t i = 0; i < offset; ++i) {
        sum += processed_peaks[i];
    }
    for (size_t i = offset; i < processed_peaks.size() + offset; ++i) {
        if (i < processed_peaks.size()) {
            sum += processed_peaks[i];
        }
        if (i >= window_size) {
            sum -= processed_peaks[i - window_size];
        }
        double old_intensity = processed_peaks[i - offset];
        double new_intensity = old_intensity - (sum - old_intensity) / (window_size - 1);
        filtered[i - offset] = new_intensity;
    }

    // add flanking peaks by default
    std::vector<double> leftshifted(filtered);
    for (size_t i = 0; i + 1 < leftshifted.size(); ++i) {
        leftshifted[i] = leftshifted[i + 1] * 0.5;
    }
    leftshifted[leftshifted.size() - 1] = 0;
    std::vector<double> rightshifted(filtered);
    for (int i = int(rightshifted.size()) - 1; i > 0; --i) {
        rightshifted[i] = rightshifted[i - 1] * 0.5;
    }
    rightshifted[0] = 0;
    for (int i = 0; i < filtered.size(); ++i) {
        filtered[i] += (leftshifted[i] + rightshifted[i]);
    }

    return filtered;
}

double XCorr(const std::vector<double>& spec, double resolution,
             const char* sequence, size_t sequence_length,
             size_t modified_site, double mass_shift, int maximum_charge) {
    double xcorr = 0;
    // merge two ion series separately
    double current_b_ion = 0;
    for (int b_ion_idx = 1; b_ion_idx <= sequence_length; ++b_ion_idx) {
        current_b_ion += MassTable.at(sequence[b_ion_idx - 1]);
        if (b_ion_idx - 1 == modified_site) { current_b_ion += mass_shift; }
        for (int charge_state = 1; charge_state <= maximum_charge; ++charge_state) {
            double mz_position = current_b_ion / float(charge_state) + PROTON_MASS;
            int index = int(mz_position / resolution);  // truncate
            if (index >= spec.size()) { continue; }
            xcorr += spec[index];
        }
    }
    double current_y_ion = WATER_MASS;
    for (int y_ion_idx = 1; y_ion_idx <= sequence_length; ++y_ion_idx) {
        current_y_ion += MassTable.at(sequence[sequence_length - y_ion_idx]);
        if (sequence_length - y_ion_idx == modified_site) { current_y_ion += mass_shift; }
        for (int charge_state = 1; charge_state <= maximum_charge; ++charge_state) {
            double mz_position = current_y_ion / float(charge_state) + PROTON_MASS;
            int index = int(mz_position / resolution);  // truncate
            if (index >= spec.size()) { continue; }
            xcorr += spec[index];
        }
    }
    return xcorr * 0.005;
}

double ModXCorr(const std::vector<double>& spec, double resolution,
                const char* sequence, size_t sequence_length,
                size_t modified_site, double mass_shift, int maximum_charge,
                const std::vector<std::pair<size_t, double>>& mods) {
    
    constexpr double MassTable[26] = {
        71.03712 /* A */, NAN /* B */, 103.00919 /* C */, 115.02695 /* D */,
        129.04260 /* E */, 147.06842 /* F */, 57.02147 /* G */, 137.05891 /* H */,
        113.08407 /* I */, NAN /* J */, 128.09497 /* K */, 113.08407 /* L */,
        131.04049 /* M */, 114.04293 /* N */, NAN /* O */, 97.05277 /* P */,
        128.05858 /* Q */, 156.10112 /* R */, 87.03203 /* S */, 101.04768 /* T */,
        NAN /* U */, 99.06842 /* V */, 186.07932 /* W */, NAN /* X */,
        163.06333 /* Y */, NAN /* Z */
    };
    constexpr double inv_charge[] = { 
        0.0, 1.0 / 1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0 
    };

    const int spec_size = static_cast<int>(spec.size());
    const int mods_size = static_cast<int>(mods.size());
    const double inv_resolution = 1 / resolution;

    double xcorr = 0.0;

    int b_mod_idx = 0;
    double current_b_ion = 0.0;
    for (int b_ion_idx = 1; b_ion_idx <= sequence_length; ++b_ion_idx) {
        current_b_ion += MassTable[sequence[b_ion_idx - 1] - 'A'];
        if (b_ion_idx - 1 == modified_site) {
            current_b_ion += mass_shift;
        }
        if (b_mod_idx < mods_size && b_ion_idx - 1 == mods[b_mod_idx].first) {
            current_b_ion += mods[b_mod_idx++].second;
        }
        for (int charge_state = 1; charge_state <= maximum_charge; ++charge_state) {
            double mz_position = current_b_ion * inv_charge[charge_state] + PROTON_MASS;
            int index = static_cast<int>(mz_position * inv_resolution);  // truncate
            if (index >= spec_size) {
                continue;
            }
            xcorr += spec[index];
        }
    }

    int y_mod_idx = mods_size - 1;  // already consider whether mods exist
    double current_y_ion = WATER_MASS;
    for (int y_ion_idx = 1; y_ion_idx <= sequence_length; ++y_ion_idx) {
        current_y_ion += MassTable[sequence[sequence_length - y_ion_idx] - 'A'];
        if (sequence_length - y_ion_idx == modified_site) {
            current_y_ion += mass_shift;
        }
        if (y_mod_idx >= 0 && sequence_length - y_ion_idx == mods[y_mod_idx].first) {
            current_y_ion += mods[y_mod_idx--].second;
        }
        // TODO: Optimize this line has to reorder the compute, lead to unstable numerical compute.
        for (int charge_state = 1; charge_state <= maximum_charge; ++charge_state) {  
            double mz_position = current_y_ion * inv_charge[charge_state] + PROTON_MASS;
            int index = static_cast<int>(mz_position * inv_resolution);  // truncate
            if (index >= spec_size) {
                continue;
            }
            xcorr += spec[index];
        }
    }

    return xcorr * 0.005;
}

#endif // XOLIK_XCORR_H
