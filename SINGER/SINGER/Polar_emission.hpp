//
//  Polar_emission.hpp
//  SINGER
//
//  Created by Yun Deng on 6/14/23.
//

#ifndef Polar_emission_hpp
#define Polar_emission_hpp

#include <stdio.h>
#include <math.h>
#include "Emission.hpp"

using namespace std;

class Polar_emission : public Emission {
    
public:
    
    double penalty = 0.01;
    double ancestral_prob = 0.5;
    double root_reward = 1;
    
    vector<double> diff = vector<double>(4);
    
    Polar_emission();
    
    ~Polar_emission();
    
    double null_emit(Branch &branch, double time, double theta, Node_ptr node) override;
    
    double mut_emit(Branch &branch, double time, double theta, double bin_size, set<double> &mut_set, Node_ptr node) override;
    
    double emit(Branch &branch, double time, double theta, double bin_size, vector<double> &emissions, Node_ptr node) override;
    
    double mut_prob(double theta, double bin_size, double ll, double lu, double l0, int sl, int su, int s0);
    
    double null_prob(double theta, double ll, double lu, double l0);
    
    double mut_prob(double theta, double bin_size, int s);
    
    double null_prob(double theta);
    
    void get_diff(double m, Branch branch, Node_ptr node);
};


#endif /* Polar_emission_hpp */
