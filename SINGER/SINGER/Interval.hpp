//
//  Interval.hpp
//  SINGER
//
//  Created by Yun Deng on 4/9/22.
//

#ifndef Interval_hpp
#define Interval_hpp

#include <stdio.h>
#include <map>
#include <numeric>
#include "Recombination.hpp"

using namespace std;

class Interval {
    
public:
    
    Branch branch = Branch();
    float lb = 0;
    float ub = 0;
    float weight = 0.0;
    float time = 0.0;
    int start_pos = 0;
    Node_ptr node = nullptr;
    float reduction = 1.0;
    
    vector<float> source_weights = {};
    vector<Interval *> source_intervals = {};
    vector<shared_ptr<Interval>> intervals = {};
    
    Interval();
    
    Interval(Branch b, float tl, float tu, int init_pos);
    
    void assign_weight(float w);
    
    void assign_time(float t);
    
    void fill_time();
    
    bool full(float t);

    bool operator<(const Interval &other) const;
    
    bool operator==(const Interval &other) const;
    
    bool operator!=(const Interval &other) const;
};

shared_ptr<Interval> create_interval(Branch b, float tl, float tu, int init_pos);

struct compare_interval {
    
    bool operator()(const Interval *i1, const Interval *i2) const {
            if (i1->branch != i2->branch) {
                return i1->branch < i2->branch;
            }
            if (i1->ub != i2->ub) {
                return i1->ub < i2->ub;
            }
            return i1->lb < i2->lb;
    }
};

struct equal_interval {
    
    bool operator()(const Interval *i1, const Interval *i2) const {
            if (i1->branch != i2->branch) {
                return false;
            }
            if (i1->ub != i2->ub) {
                return false;
            }
            return i1->lb == i2->lb;
    }
};


class Interval_info {
    
public:
    
    Branch branch;
    float lb = 0;
    float ub = 0;
    float time = 0;
    float seed_pos = 0;
    
    Interval_info();
    
    Interval_info(Branch b, float tl, float tu);
    
    bool operator==(const Interval_info& other) const;

    bool operator!=(const Interval_info& other) const;

    bool operator<(const Interval_info& other) const;
    
};

#endif /* Interval_hpp */
