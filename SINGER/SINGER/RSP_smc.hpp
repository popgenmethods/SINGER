//
//  RSP_smc.hpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#ifndef RSP_smc_hpp
#define RSP_smc_hpp

#include <stdio.h>
#include <math.h>
#include <map>
#include "random_utils.hpp"
#include "Branch.hpp"
#include "Tree.hpp"
#include "Recombination.hpp"

class RSP_smc {
    
public:
    
    RSP_smc();
    
    void sample_recombination(Recombination &r, double cut_time, Tree &tree);
    
    void approx_sample_recombination(Recombination &r, double cut_time);
    
    void adjust(Recombination &r, double cut_time);
    
    void approx_sample_recombination(Recombination &r, double cut_time, double n);
    
    void adjust(Recombination &r, double cut_time, double n);
    
private:
    
    map<double, int> coalescence_rates = {};
    
    double sample_start_time(Branch b, int density, double join_time, double cut_time);
    
    pair<Branch, double> sample_start_time(Branch b1, Branch b2, int density, double join_time, double cut_time);
    
    void get_coalescence_rate(Tree &tree, Recombination &r, double cut_time);
    
    double recomb_pdf(double s, double t);
    
    double random_time(double lb, double ub);
    
    double random_time(double lb, double ub, double q);
    
    double choose_time(double lb, double ub);
    
    double choose_time(double lb, double ub, double n);
    
};

#endif /* RSP_smc_hpp */
