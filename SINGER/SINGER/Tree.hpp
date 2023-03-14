//
//  Tree.hpp
//  SINGER
//
//  Created by Yun Deng on 4/12/22.
//

#ifndef Tree_hpp
#define Tree_hpp

#include <stdio.h>
#include <math.h>
#include "Branch.hpp"
#include "Recombination.hpp"

class Tree {

public:
    
    set<Branch> branches = {};
    
    Tree();
    
    float length();
    
    void insert_branch(Branch b);
    
    void delete_branch(Branch b);
    
    void forward_update(Recombination &r);
    
    void backward_update(Recombination &r);
    
    Branch find_split_branch(Branch removed_branch);
    
    pair<Branch, float> sample_cut_point();
    
    float prior_likelihood();
    
    float data_likelihood(float theta, float pos);
    
    float null_likelihood(float theta);
    
    float data_likelihood(float theta, float bin_size, set<float> mutations);
    
    float transition_likelihood(Recombination& r);
    
private:
    
    float tree_length = 0.0f;
    
    float log_exp(float lambda, float x);
    
    float random();
    
};

#endif /* Tree_hpp */