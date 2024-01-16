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
    
    bool operator()(const pair<double, double> p1, const pair<double, double> p2) const {
        return p1.first < p2.first;
    }
    
};

struct compare_prob {
    
    bool operator()(const pair<double, double> p1, const pair<double, double> p2) const {
            if (p1.second != p2.second) {
                return p1.second < p2.second;
            }
        return p1.first < p2.first;
    }
    
};

class Coalescent_calculator {
    
public:
    
    double cut_time;
    double min_time, max_time;
    map<double, int> rate_changes = {};
    map<double, int> rates = {};
    set<pair<double, double>, compare_time> probs = {};
    set<pair<double, double>, compare_prob> quantiles = {};
    
    Coalescent_calculator(double t);
    
    ~Coalescent_calculator();
    
    void compute(set<Branch> &branches);
    
    double weight(double lb, double ub);
    
    double time(double lb, double ub);
    
// private:
    
    void compute_rate_changes(set<Branch> &branches);
    
    void compute_rates();
    
    void compute_probs_quantiles();
    
    double prob(double x);
    
    double quantile(double p);
    
};

#endif /* Coalescent_calculator_hpp */
