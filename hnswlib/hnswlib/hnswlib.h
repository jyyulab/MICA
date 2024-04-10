#pragma once

// https://github.com/nmslib/hnswlib/pull/508
// This allows others to provide their own error stream (e.g. RcppHNSW)
#ifndef HNSWLIB_ERR_OVERRIDE
  #define HNSWERR std::cerr
#else
  #define HNSWERR HNSWLIB_ERR_OVERRIDE
#endif


#include <queue>
#include <vector>
#include <iostream>
#include <string.h>

namespace hnswlib {
typedef size_t labeltype;

// This can be extended to store state for filtering (e.g. from a std::set)
class BaseFilterFunctor {
 public:
    virtual bool operator()(hnswlib::labeltype id) { return true; }
    virtual ~BaseFilterFunctor() {};
};

template<typename dist_t>
class BaseSearchStopCondition {
 public:
    virtual void add_point_to_result(labeltype label, const void *datapoint, dist_t dist) = 0;

    virtual void remove_point_from_result(labeltype label, const void *datapoint, dist_t dist) = 0;

    virtual bool should_stop_search(dist_t candidate_dist, dist_t lowerBound) = 0;

    virtual bool should_consider_candidate(dist_t candidate_dist, dist_t lowerBound) = 0;

    virtual bool should_remove_extra() = 0;

    virtual void filter_results(std::vector<std::pair<dist_t, labeltype >> &candidates) = 0;

    virtual ~BaseSearchStopCondition() {}
};

template <typename T>
class pairGreater {
 public:
    bool operator()(const T& p1, const T& p2) {
        return p1.first > p2.first;
    }
};

template<typename T>
static void writeBinaryPOD(std::ostream &out, const T &podRef) {
    out.write((char *) &podRef, sizeof(T));
}

template<typename T>
static void readBinaryPOD(std::istream &in, T &podRef) {
    in.read((char *) &podRef, sizeof(T));
}

template<typename MTYPE>
using DISTFUNC = MTYPE(*)(const void *, const void *, const void *);

template<typename MTYPE>
class SpaceInterface {
 public:
    // virtual void search(void *);
    virtual size_t get_data_size() = 0;

    virtual DISTFUNC<MTYPE> get_dist_func() = 0;

    virtual void *get_dist_func_param() = 0;

    virtual ~SpaceInterface() {}
};

template<typename dist_t>
class AlgorithmInterface {
 public:
    virtual void addPoint(const void *datapoint, labeltype label, bool replace_deleted = false) = 0;

    virtual std::priority_queue<std::pair<dist_t, labeltype>>
        searchKnn(const void*, size_t, BaseFilterFunctor* isIdAllowed = nullptr) const = 0;

    // Return k nearest neighbor in the order of closer fist
    virtual std::vector<std::pair<dist_t, labeltype>>
        searchKnnCloserFirst(const void* query_data, size_t k, BaseFilterFunctor* isIdAllowed = nullptr) const;

    virtual void saveIndex(const std::string &location) = 0;
    virtual ~AlgorithmInterface(){
    }
};

template<typename dist_t>
std::vector<std::pair<dist_t, labeltype>>
AlgorithmInterface<dist_t>::searchKnnCloserFirst(const void* query_data, size_t k,
                                                 BaseFilterFunctor* isIdAllowed) const {
    std::vector<std::pair<dist_t, labeltype>> result;

    // here searchKnn returns the result in the order of further first
    auto ret = searchKnn(query_data, k, isIdAllowed);
    {
        size_t sz = ret.size();
        result.resize(sz);
        while (!ret.empty()) {
            result[--sz] = ret.top();
            ret.pop();
        }
    }

    return result;
}
}  // namespace hnswlib

size_t bin_size_mi = 0;
size_t bin_power_mi = 0;

void set_bin_size(size_t bins){
    bin_size_mi = bins;
    if (bin_power_mi != 0){
        bin_power_mi = 0;
    }
}

void set_bin_power(size_t bins){
    bin_power_mi = bins;
    if (bin_size_mi != 0){
        bin_size_mi = 0;
    }
}

#include "space_l2.h"
#include "space_ip.h"
#include "stop_condition.h"
#include "bruteforce.h"
#include "hnswalg.h"
