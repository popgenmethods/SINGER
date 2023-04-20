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

class Coalescent_calculator {
    
public:
    
    float cut_time;
    // set<Branch> branches = {};
    map<float, int> rate_changes = {};
    map<float, int> rates = {};
    map<float, float> probs = {};
    map<float, float> quantiles = {};
    
    Coalescent_calculator(float t);
    
    ~Coalescent_calculator();
    
    // void start(set<Branch> &inserted_branches);
    
    // void update(set<Branch> &deleted_branches, set<Branch> &inserted_branches);
    
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
