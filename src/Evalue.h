#ifndef XOLIK_EVALUE_H
#define XOLIK_EVALUE_H

#include <vector>
#include <cmath>

// return y = a + b * x as well as a flag to indicate whether the fitting succeeds
bool LeastSquares(const std::vector<double>& y, const std::vector<double>& x,
                  double& a_out, double& b_out) {
    if (y.size() != x.size() || y.empty()) {
        return false;
    }
    int size = y.size();
    double sum_x = 0.0;
    double sum_y = 0.0;
    for (int i = 0; i < size; ++i) {
        sum_x += x[i];
        sum_y += y[i];
    }
    double mean_x = sum_x / double(size);
    double mean_y = sum_y / double(size);
    double x2 = 0.0;
    double xy = 0.0;
    for (int i = 0; i < size; ++i) {
        double dx = x[i] - mean_x;
        double dy = y[i] - mean_y;
        x2 += dx * dx;
        xy += dx * dy;
    }
    b_out = xy / x2;
    a_out = mean_y - b_out * mean_x;
    return true;
}

// Comet's estimation procedure
double CalculateEValue(double target, const std::vector<double>& scores) {
    using std::vector;

    // generate histogram, TODO: do you need a different way to generate histogram?
    const size_t histogram_size = 152;
    vector<int> histogram(histogram_size, 0);
    int max_idx = 0;
    for (double xcorr : scores) {
        int idx = int(xcorr * 10.0 + 0.05);
        if (idx >= histogram_size) {
            idx = histogram_size - 1;
        }
        histogram[idx] += 1;
        max_idx = idx > max_idx ? idx : max_idx;
    }

    // calculate ending position
    bool from_zero = true;
    int end_idx = 0;
    int i = 0;
    for (; i < max_idx; ++i) {
        if (from_zero && histogram[i] > 0) {
            from_zero = false;
        }
        if (!from_zero && histogram[i] == 0) {
            if (i + 1 < histogram_size && histogram[i+1] == 0) {
                end_idx = i > 0 ? i - 1 : 0;
                break;
            }
        }
    }
    if (i == max_idx) {
        end_idx = max_idx > 12 ? max_idx - 2 : max_idx;
    }

    // generate survival
    vector<double> survival(histogram_size, 0.0);
    survival[end_idx] = histogram[end_idx];
    for (int j = end_idx - 1; j >= 0; --j) {
        survival[j] = survival[j+1] + histogram[j];
        if (histogram[j+1] == 0) {  // TODO: strange
            survival[j+1] = 0.0;
        }
    }

    // log10 transformation
    vector<double> log10survival(histogram_size, -999.0);
    for (int j = 0; j <= end_idx; ++j) {
        log10survival[j] = log10(survival[j]);
    }
    int start_idx = 0;
    if (end_idx >= 30) {
        start_idx = int(end_idx * 0.75);
    }
    else if (end_idx >= 15) {
        start_idx = int(end_idx * 0.5);
    }

    // prepare data points for fitting
    vector<double> y;
    vector<double> x;
    double a = 0.0;
    double b = 0.0;
    for (int j = start_idx; j <= end_idx; ++j) {
        if (survival[j] > 0) {
            x.push_back(j);
            y.push_back(log10survival[j]);
        }
    }
    while (!LeastSquares(y, x, a, b) || (start_idx > 0 && b >= 0)) {
        --start_idx;
        x.push_back(start_idx);
        y.push_back(log10survival[start_idx]);
    }

    // calculate the exponent of the e value
    double exponent = a + b * int(target * 10.0 + 0.05);
//    double evalue = pow(10.0, a + b * int(target * 10.0 + 0.05));
//    return evalue > 999.0 ? 999.0 : evalue;
    return exponent;
}

#endif // XOLIK_EVALUE_H
