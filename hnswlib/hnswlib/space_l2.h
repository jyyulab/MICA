#pragma once
#include "hnswlib.h"
#include <algorithm>
#include <cmath>

namespace hnswlib {

float entropy(const std::vector<float>& probabilities) {
    float entropy = 0.0;
    for (float prob : probabilities) {
        if (prob > 0) {
            entropy -= prob * log(prob);
        }
    }
    return entropy;
}

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& input) {
    std::vector<T> output;
    for (const auto& row : input) {
        output.insert(output.end(), row.begin(), row.end());
    }
    return output;
}

static float
L2Sqr(const void *pVect1v, const void *pVect2v, const void *qty_ptr) {
    float *x = (float *) pVect1v;
    float *y = (float *) pVect2v;
    size_t qty = *((size_t *) qty_ptr);
    size_t bins;
    if (bin_size_mi != 0){
        bins = bin_size_mi;
    }
    else if (bin_power_mi != 0){
        bins = static_cast<size_t>(pow(qty, 1.0 / bin_power_mi));
    }
    else{
        bins = static_cast<size_t>(pow(qty, 1.0 / 3.0));
    }
    // size_t bins = static_cast<size_t>(pow(qty, 1.0 / 3.0));
    

    // float res = 0;
    // for (size_t i = 0; i < qty; i++) {
    //     float t = *pVect1 - *pVect2;
    //     pVect1++;
    //     pVect2++;
    //     res += t * t;
    // }
    // return (res);

    // Create 2D histogram
    std::vector<std::vector<int>> hist(bins, std::vector<int>(bins, 0));

    // Populate histogram
    //for (size_t i = 0; i < qty; ++i) {
        //int bin_x = static_cast<int>(floor(x[i] * bins));
        //int bin_y = static_cast<int>(floor(y[i] * bins));

        // Ensure bins are within range
        //bin_x = std::max(0, std::min(static_cast<int>(bins - 1), bin_x));
        //bin_y = std::max(0, std::min(static_cast<int>(bins - 1), bin_y));

        //hist[bin_x][bin_y]++;
    //}
    float min_x = *std::min_element(x, x + qty);
    float max_x = *std::max_element(x, x + qty);
    float min_y = *std::min_element(y, y + qty);
    float max_y = *std::max_element(y, y + qty);

    // Calculate bin widths
    float bin_width_x = (max_x - min_x) / bins;
    float bin_width_y = (max_y - min_y) / bins;

    // Populate histogram
    for (size_t i = 0; i < qty; ++i) {
        int bin_x = static_cast<int>((x[i] - min_x) / bin_width_x);
        int bin_y = static_cast<int>((y[i] - min_y) / bin_width_y);

        // Ensure bins are within range
        bin_x = std::max(0, std::min(static_cast<int>(bins - 1), bin_x));
        bin_y = std::max(0, std::min(static_cast<int>(bins - 1), bin_y));

        hist[bin_x][bin_y]++;
    }

    // Calculate probabilities
    float total_samples = static_cast<float>(qty);
    std::vector<float> px(bins, 0.0);
    std::vector<float> py(bins, 0.0);
    std::vector<std::vector<float>> pxy(bins, std::vector<float>(bins, 0.0));

    for (size_t i = 0; i < bins; ++i) {
        for (size_t j = 0; j < bins; ++j) {
            px[i] += hist[i][j] / total_samples;
            py[j] += hist[i][j] / total_samples;
            pxy[i][j] = hist[i][j] / total_samples;
        }
    }

    // Calculate entropies
    float hxy = entropy(flatten(pxy));
    float hx = entropy(px);
    float hy = entropy(py);
    //if(hxy==0){
        //return (1.0);
    //}
    // Calculate mutual information
    float mutual_info = 2 - (hx + hy)/hxy;
    return (mutual_info);
}




class L2Space : public SpaceInterface<float> {
    DISTFUNC<float> fstdistfunc_;
    size_t data_size_;
    size_t dim_;

 public:
    L2Space(size_t dim) {
        fstdistfunc_ = L2Sqr;
        dim_ = dim;
        data_size_ = dim * sizeof(float);
    }

    size_t get_data_size() {
        return data_size_;
    }

    DISTFUNC<float> get_dist_func() {
        return fstdistfunc_;
    }

    void *get_dist_func_param() {
        return &dim_;
    }

    ~L2Space() {}
};

} // namespace hnswlib
