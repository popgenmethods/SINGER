//
//  RSP.hpp
//  SINGER
//
//  Created by Yun Deng on 10/13/22.
//

#ifndef RSP_hpp
#define RSP_hpp

#include <stdio.h>
#include <math.h>
#include <map>
#include "Branch.hpp"
#include "Tree.hpp"

class RSP {
    
public:
    
    virtual void exact_sample_recombination(Recombination &r, float cut_time, Tree tree) = 0;
    
private:
    
    map<float, int> coalescence_rates = {};
    
    float exact_sample_start_time(Branch b, int density, float join_time, float cut_time);
    
    pair<Branch, float> exact_sample_start_time(Branch b1, Branch b2, int density, float join_time, float cut_time);
    
    void get_coalescence_rate(Tree tree, Recombination r, float cut_time);
    
    float recomb_pdf(float s, float t);
    
    float exact_recomb_pdf(float s, float t);
    
    float random();
    
};

#endif /* RSP_hpp */
