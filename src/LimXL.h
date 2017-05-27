#pragma once

#include <tuple>
#include <vector>
#include <deque>
#include <algorithm>
#include "Scores.h"

typedef size_t CandIdx;
typedef double Score;
typedef double Mass;

std::tuple<CandIdx, CandIdx, Score> LimXLMatch(const std::vector<double>& masses, Scores& scores,
                                               double precursor_mass, double xlinker_mass, double threshold,
                                               double left_tol, double right_tol) {
    auto global_max_info = std::make_tuple(scores.size(), scores.size(), 0.0);  // out of range means no
    auto max_allowed = (precursor_mass - xlinker_mass + right_tol) / 2;
    auto last_iter = std::upper_bound(masses.begin(), masses.begin() + scores.size(), max_allowed);
    auto last_idx = std::distance(masses.begin(), last_iter);

    int forward_idx = 0;
    int backward_front = scores.size();
    int backward_end = scores.size() - 1;
    std::deque<size_t> deque;  // store indexes
    while (forward_idx < last_idx) {

        auto alpha_mass = masses[forward_idx];
        auto lower_bound = precursor_mass - alpha_mass - xlinker_mass - left_tol;
        auto upper_bound = precursor_mass - alpha_mass - xlinker_mass + right_tol;

        // Step 1: move backward_end and make sure elements in the front of deque are
        // within the allowed mass range, backward_end included
        while (0 <= backward_end && masses[backward_end] > upper_bound) {
            --backward_end;
        }
        while (!deque.empty() && deque.front() > backward_end) {  // compare indexes
            deque.pop_front();
        }

        if (backward_front > backward_end + 1) {
            backward_front = backward_end + 1;
        }

        // Step 2: move backward_front to append element to deque
        // backward_front is included in the deque
        while (0 < backward_front && masses[backward_front - 1] >= lower_bound) {
            // append backward_front
            while (!deque.empty() && scores[backward_front - 1] > scores[deque.back()]) {
                deque.pop_back();
            }
            deque.push_back(backward_front - 1);
            --backward_front;
        }

        // Step 3: get maximum, front of the deque
        // because of this specific case, it is possible that no element in the deque
        if (!deque.empty()) {  // if deque.empty() means
            if (scores[deque.front()] < threshold || scores[forward_idx] < threshold) { ++forward_idx; continue; }
            auto local_max = scores[forward_idx] + scores[deque.front()];
            if (std::get<0>(global_max_info) == scores.size()
                    || local_max > std::get<2>(global_max_info)) {
                global_max_info = std::make_tuple(forward_idx, deque.front(), local_max);
            }
        }
        ++forward_idx;
    }
    return global_max_info;
}

std::tuple<CandIdx, CandIdx, Score> NaiveMatch(const std::vector<double>& masses, Scores& scores,
                                               double precursor_mass, double xlinker_mass, double threshold,
                                               double left_tol, double right_tol,
                                               std::vector<double>& collected_scores) {
    auto global_max_info = std::make_tuple(scores.size(), scores.size(), 0.0);  // out of range means no
    auto max_allowed = (precursor_mass - xlinker_mass + right_tol) / 2;
    auto last_iter = std::upper_bound(masses.begin(), masses.begin() + scores.size(), max_allowed);
    auto last_idx = std::distance(masses.begin(), last_iter);

    for (auto i = 0; i < last_idx; ++i) {

        auto alpha_mass = masses[i];
        auto lower_bound = precursor_mass - alpha_mass - xlinker_mass - left_tol;
        auto upper_bound = precursor_mass - alpha_mass - xlinker_mass + right_tol;
        auto start = std::lower_bound(masses.begin(), masses.begin() + scores.size(), lower_bound);
        auto end = std::upper_bound(masses.begin(), masses.begin() + scores.size(), upper_bound);
        auto start_idx = std::distance(masses.begin(), start);
        auto end_idx = std::distance(masses.begin(), end);

        if (start_idx == end_idx) { continue; }
        auto local_max_idx = end_idx - 1;
        for (int j = end_idx - 1; j >= start_idx; --j) {
            collected_scores.push_back(scores[i] + scores[j]);
            if (scores[j] > scores[local_max_idx]) {
                local_max_idx = j;
            }
        }

        if (scores[local_max_idx] < threshold || scores[i] < threshold) { continue; }

        auto local_max = scores[i] + scores[local_max_idx];
        if (std::get<0>(global_max_info) == scores.size()
                || local_max > std::get<2>(global_max_info)) {
            global_max_info = std::make_tuple(i, local_max_idx, local_max);
        }
    }
    return global_max_info;
}
