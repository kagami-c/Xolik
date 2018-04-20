#ifndef XOLIK_MATCH_H
#define XOLIK_MATCH_H

#include <tuple>
#include <vector>
#include <deque>
#include <algorithm>
#include "ScoreArray.h"

typedef size_t CandIdx;
typedef double Score;

std::tuple<CandIdx, CandIdx, Score> XolikMatch(const std::vector<double>& mass_array, 
                                               ScoreArray& score_array, double precursor_mass, 
                                               double xlinker_mass, double threshold, 
                                               double left_tol, double right_tol, int& count_out) {

    const size_t scores_size = score_array.size();
    auto global_max_info = std::make_tuple(scores_size, scores_size, 0.0);  // out of range means no
    double max_allowed = (precursor_mass - xlinker_mass + right_tol) / 2;
    auto last_iter = std::upper_bound(mass_array.begin(), mass_array.begin() + scores_size, max_allowed);
    int last_idx = std::distance(mass_array.begin(), last_iter);

    int forward_idx = 0;
    int backward_front = static_cast<int>(scores_size);
    int backward_end = static_cast<int>(scores_size) - 1;
    std::deque<int> deque;  // store indexes
    while (forward_idx < last_idx) {

        auto alpha_mass = mass_array[forward_idx];
        auto lower_bound = precursor_mass - alpha_mass - xlinker_mass - left_tol;
        auto upper_bound = precursor_mass - alpha_mass - xlinker_mass + right_tol;

        // Step 1: move backward_end and make sure elements in the front of deque are
        // within the allowed mass range, backward_end included
        while (0 <= backward_end && mass_array[backward_end] > upper_bound) {
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
        while (0 < backward_front && mass_array[backward_front - 1] >= lower_bound) {
            // append backward_front
            while (!deque.empty() && score_array[backward_front - 1] > score_array[deque.back()]) {
                deque.pop_back();
            }
            deque.push_back(backward_front - 1);
            --backward_front;
        }

        // optional: to support evalue estimation
        count_out += (backward_end + 1 - backward_front);

        // Step 3: get maximum, front of the deque
        // because of this specific case, it is possible that no element in the deque
        if (!deque.empty()) {  // if deque.empty() means no match
            if (score_array[deque.front()] < threshold || score_array[forward_idx] < threshold) {
                ++forward_idx; 
                continue;
            }
            auto local_max = score_array[forward_idx] + score_array[deque.front()];
            if (std::get<0>(global_max_info) == scores_size
                    || local_max >= std::get<2>(global_max_info)) {
                global_max_info = std::make_tuple(forward_idx, deque.front(), local_max);
            }
        }
        ++forward_idx;
    }
    return global_max_info;
}

std::tuple<CandIdx, CandIdx, Score> NaiveMatch(const std::vector<double>& mass_array, 
                                               ScoreArray& score_array, double precursor_mass, 
                                               double xlinker_mass, double threshold, 
                                               double left_tol, double right_tol, int& count_out, 
                                               std::vector<double>& collected_scores, bool collect,
                                               int collect_size) {

    const size_t scores_size = score_array.size();
    auto global_max_info = std::make_tuple(scores_size, scores_size, 0.0);  // out of range means no
    double max_allowed = (precursor_mass - xlinker_mass + right_tol) / 2;
    auto last_iter = std::upper_bound(mass_array.begin(), mass_array.begin() + scores_size, max_allowed);
    int last_idx = std::distance(mass_array.begin(), last_iter);

    for (auto i = 0; i < last_idx; ++i) {

        auto alpha_mass = mass_array[i];
        auto lower_bound = precursor_mass - alpha_mass - xlinker_mass - left_tol;
        auto upper_bound = precursor_mass - alpha_mass - xlinker_mass + right_tol;
        auto start = std::lower_bound(mass_array.begin(), mass_array.begin() + scores_size, lower_bound);
        auto end = std::upper_bound(mass_array.begin(), mass_array.begin() + scores_size, upper_bound);
        int start_idx = std::distance(mass_array.begin(), start);
        int end_idx = std::distance(mass_array.begin(), end);

        if (start_idx == end_idx) { continue; }
        auto local_max_idx = end_idx - 1;
        for (int j = end_idx - 1; j >= start_idx; --j) {

            // collect mode
            if (collect) {
                collected_scores.push_back(score_array[i] + score_array[j]);
                if (collected_scores.size() >= collect_size) {
                    return global_max_info;
                }
            }

            // optional: to support evalue estimation
            ++count_out;

            if (score_array[j] > score_array[local_max_idx]) {
                local_max_idx = j;
            }
        }

        if (score_array[local_max_idx] < threshold || score_array[i] < threshold) {
            continue;
        }

        auto local_max = score_array[i] + score_array[local_max_idx];
        if (std::get<0>(global_max_info) == scores_size
                || local_max >= std::get<2>(global_max_info)) {
            global_max_info = std::make_tuple(i, local_max_idx, local_max);
        }
    }
    return global_max_info;
}

#endif // XOLIK_MATCH_H
