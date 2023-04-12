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
    
    Branch branch;
    float lb = 0;
    float ub = 0;
    float weight = 0.0;
    float time = 0.0;
    float start_pos = 0;
    vector<float> source_weights = {};
    vector<Interval *> source_intervals = {};
    Node *node = nullptr;
    float reduction = 1.0;
    
    Interval(Branch b, float tl, float tu, float init_pos);
    
    void assign_weight(float w);
    
    void assign_time(float t);
    
    void fill_time();
    
    bool full(float t);

    void set_source(vector<Interval *> intervals, vector<float> weights);
    
    void set_node(Node *n);
    
    Interval *sample_source();
};

struct compare_interval {
    
    bool operator()(const Interval *i1, const Interval *i2) const {
        if (i1->branch != i2->branch) {
            return i1->branch < i2->branch;
        }
        return i1->time < i2->time;
    }
};

class Interval_info {
    
public:
    
    Branch branch;
    float lb;
    float ub;
    
    Interval_info();
    
    Interval_info(Branch b, float tl, float tu);
    
};

bool operator==(const Interval_info& i1, const Interval_info& i2);

bool operator!=(const Interval_info& i1, const Interval_info& i2);

bool operator<(const Interval_info& i1, const Interval_info& i2);

#endif /* Interval_hpp */
