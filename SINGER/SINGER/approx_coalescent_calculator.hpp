//
//  approx_coalescent_calculator.hpp
//  SINGER
//
//  Created by Yun Deng on 6/25/23.
//

#ifndef approx_coalescent_calculator_hpp
#define approx_coalescent_calculator_hpp

#include <stdio.h>
#include <map>
#include <math.h>
#include "Branch.hpp"
#include "Tree.hpp"
#include "Recombination.hpp"

class approx_coalescent_calculator {
    
public:
    
    float cut_time;
    float rho = 0;
    float first_moment = 0;
    int n0 = 0;
    
    approx_coalescent_calculator(float t);
    
    ~approx_coalescent_calculator();
    
    void start(set<Branch> &branches);
    
    void start(Tree &tree);
    
    void update(Recombination &r);
    
    pair<float, float> compute_time_weights(float x, float y);
    
    void compute_first_moment();
    
    float prob(float x, float y);
    
    float prob_integral(float x);
    
    float find_median(float x, float y);
    
};

#endif /* approx_coalescent_calculator_hpp */
