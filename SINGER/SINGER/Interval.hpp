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
    double lb = 0;
    double ub = 0;
    double weight = 0.0;
    double time = 0.0;
    int start_pos = 0;
    int source_pos = 0;
    Node_ptr node = nullptr;
    double reduction = 1.0;
    
    vector<double> source_weights = {};
    vector<Interval *> source_intervals = {};
    vector<shared_ptr<Interval>> intervals = {};
    
    Interval();
    
    Interval(Branch b, double tl, double tu, int init_pos);
    
    void assign_weight(double w);
    
    void assign_time(double t);
    
    void fill_time();
    
    bool full(double t);

    bool operator<(const Interval &other) const;
    
    bool operator==(const Interval &other) const;
    
    bool operator!=(const Interval &other) const;
};

shared_ptr<Interval> create_interval(Branch b, double tl, double tu, int init_pos);

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
    double lb = 0;
    double ub = 0;
    double time = 0;
    double seed_pos = 0;
    
    Interval_info();
    
    Interval_info(Branch b, double tl, double tu);
    
    bool operator==(const Interval_info& other) const;

    bool operator!=(const Interval_info& other) const;

    bool operator<(const Interval_info& other) const;
    
};

#endif /* Interval_hpp */
