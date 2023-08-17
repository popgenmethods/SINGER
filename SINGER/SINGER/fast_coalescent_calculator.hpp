//
//  fast_coalescent_calculator.hpp
//  SINGER
//
//  Created by Yun Deng on 6/15/23.
//

#ifndef fast_coalescent_calculator_hpp
#define fast_coalescent_calculator_hpp

#include <stdio.h>
#include <map>
#include <math.h>
#include "Branch.hpp"
#include "Recombination.hpp"

class fast_coalescent_calculator {
    
public:
    
    float cut_time;
    float first_moment = 0;
    float rho = 0;
    multiset<float> coalescence_times = {};
    
    fast_coalescent_calculator(float t);
    
    ~fast_coalescent_calculator();
    
    void start(set<Branch> &branches);
    
    void update(Recombination &r);
    
    void compute_first_moment();
    
    pair<float, float> compute_time_weights(float x, float y);
    
    float prob(float x, float y);
    
    float get_num_lineages(float x);
    
    float get_integral(float x);

};

#endif /* fast_coalescent_calculator_hpp */
