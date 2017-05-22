#pragma once

// novel preprocessing and scoring
std::vector<double> Preprocess(const std::vector<std::pair<double, double>>& peaks,
                               double resolution) {
    if (peaks.empty()) { return std::vector<double>(); }
    size_t vector_size = size_t(peaks[peaks.size() - 1].first / resolution) + 1;
    std::vector<double> processed_peaks(vector_size);

    // normalize in local range
    auto min_mass = peaks[0].first;
    auto max_mass = peaks[peaks.size() - 1].first;
    const auto num_step = 10;
    auto step_unit = (max_mass - min_mass) / 10;

    int start = 0;
    int end = 0;  // exclude
    int step = 0;
    do {
        auto lower_bound = min_mass + step * step_unit;
        auto upper_bound = lower_bound + step_unit;
        start = end;
        while (end < peaks.size() && peaks[end].first <= upper_bound) {
            ++end;
        }
        // from start to end (excluded) is the range for normalization
        if (start == end) { continue; }
        double local_max_intensity = 0;
        for (auto i = start; i < end; ++i) {
            if (peaks[i].second > local_max_intensity) {
                local_max_intensity = peaks[i].second;
            }
        }

        // fill processed peaks
        auto normalized_factor = 50 / sqrt(local_max_intensity);
        for (auto i = start; i < end; ++i) {
            auto normalized_intensity = sqrt(peaks[i].second) * normalized_factor;
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
