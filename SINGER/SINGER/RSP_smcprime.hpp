//
//  RSP_smcprime.hpp
//  SINGER
//
//  Created by Yun Deng on 3/19/23.
//

#ifndef RSP_smcprime_hpp
#define RSP_smcprime_hpp

#include <stdio.h>
#include <math.h>
#include <map>
#include "Branch.hpp"
#include "Tree.hpp"
#include "Recombination.hpp"

class RSP_smcprime {
    
public:
    
    RSP_smcprime();
    
    void sample_recombination(Recombination &r, float cut_time, Tree tree);
    
private:
    
    map<float, int> coalescence_rates = {};
    
    float sample_start_time(Branch b, int density, float join_time, float cut_time);
    
    pair<Branch, float> sample_start_time(Branch b1, Branch b2, int density, float join_time, float cut_time);
    
    void get_coalescence_rate(Tree tree, Recombination r, float cut_time);
    
    float recomb_pdf(float s, float t);
    
    float random();
    
};

#endif /* RSP_smcprime_hpp */
