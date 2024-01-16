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
    
    double cut_time;
    double first_moment = 0;
    double rho = 0;
    multiset<double> coalescence_times = {};
    
    fast_coalescent_calculator(double t);
    
    ~fast_coalescent_calculator();
    
    void start(set<Branch> &branches);
    
    void update(Recombination &r);
    
    void compute_first_moment();
    
    pair<double, double> compute_time_weights(double x, double y);
    
    double prob(double x, double y);
    
    double get_num_lineages(double x);
    
    double get_integral(double x);

};

#endif /* fast_coalescent_calculator_hpp */
