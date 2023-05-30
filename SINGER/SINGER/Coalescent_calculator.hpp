//
//  Coalescent_calculator.hpp
//  SINGER
//
//  Created by Yun Deng on 4/17/23.
//

#ifndef Coalescent_calculator_hpp
#define Coalescent_calculator_hpp

#include <stdio.h>
#include <map>
#include <math.h>
#include "Branch.hpp"

struct compare_time {
    
    bool operator()(const pair<float, float> p1, const pair<float, float> p2) const {
        return p1.first < p2.first;
    }
    
};

struct compare_prob {
    
    bool operator()(const pair<float, float> p1, const pair<float, float> p2) const {
            if (p1.second != p2.second) {
                return p1.second < p2.second;
            }
        return p1.first < p2.first;
    }
    
};

class Coalescent_calculator {
    
public:
    
    float cut_time;
    float min_time, max_time;
    map<float, int> rate_changes = {};
    map<float, int> rates = {};
    set<pair<float, float>, compare_time> probs = {};
    set<pair<float, float>, compare_prob> quantiles = {};
    
    Coalescent_calculator(float t);
    
    ~Coalescent_calculator();
    
    void compute(set<Branch> &branches);
    
    float weight(float lb, float ub);
    
    float time(float lb, float ub);
    
// private:
    
    void compute_rate_changes(set<Branch> &branches);
    
    void compute_rates();
    
    void compute_probs_quantiles();
    
    float prob(float x);
    
    float quantile(float p);
    
};

#endif /* Coalescent_calculator_hpp */
