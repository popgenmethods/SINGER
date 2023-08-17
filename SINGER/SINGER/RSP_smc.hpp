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
    
    void sample_recombination(Recombination &r, float cut_time, Tree &tree);
    
    void approx_sample_recombination(Recombination &r, float cut_time);
    
private:
    
    map<float, int> coalescence_rates = {};
    
    float sample_start_time(Branch b, int density, float join_time, float cut_time);
    
    pair<Branch, float> sample_start_time(Branch b1, Branch b2, int density, float join_time, float cut_time);
    
    void get_coalescence_rate(Tree &tree, Recombination &r, float cut_time);
    
    float recomb_pdf(float s, float t);
    
    float random();
    
    float random_time(float lb, float ub);
    
    float random_time(float lb, float ub, float q);
    
};

#endif /* RSP_smc_hpp */
