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
    
    float penalty = 0.001;
    float ancestral_prob = 0.99;
    float reverse_penalty = 1e-4;
    float root_reward = 1;
    
    vector<float> diff = vector<float>(4);
    
    map<float, float> num_unmapped = {};
    float num_muts = 0;
    
    Polar_emission();
    
    ~Polar_emission();
    
    float null_emit(Branch &branch, float time, float theta, Node_ptr node) override;
    
    float mut_emit(Branch &branch, float time, float theta, float bin_size, set<float> &mut_set, Node_ptr node) override;
    
    float emit(Branch &branch, float time, float theta, float bin_size, vector<float> &emissions, Node_ptr node) override;
    
    float mut_prob(float theta, float bin_size, float ll, float lu, float l0, int sl, int su, int s0);
    
    float null_prob(float theta, float ll, float lu, float l0);
    
    float mut_prob(float theta, float bin_size, int s);
    
    float null_prob(float theta);
    
    void get_diff(float m, Branch branch, Node_ptr node);
};


#endif /* Polar_emission_hpp */
