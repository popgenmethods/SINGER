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
#include <map>
#include <unordered_map>
#include "random_utils.hpp"
#include "Branch.hpp"
#include "Recombination.hpp"

class Tree {

public:
    
    set<Branch> branches = {};
    unordered_map<Node_ptr , Node_ptr > parents = {};
    
    Tree();
    
    float length();
    
    void insert_branch(Branch b);
    
    void delete_branch(Branch b);
    
    void forward_update(Recombination &r);
    
    void backward_update(Recombination &r);
    
    void remove(Branch b, Node_ptr n);
    
    void add(Branch added_branch, Branch joining_branch, Node_ptr n);
    
    Node_ptr find_sibling(Node_ptr n);
    
    Branch find_joining_branch(Branch removed_branch);
    
    pair<Branch, float> sample_cut_point();
    
    float prior_likelihood();
    
    float data_likelihood(float theta, float pos);
    
    float null_likelihood(float theta);
    
    float data_likelihood(float theta, float bin_size, set<float> mutations);
    
    float transition_likelihood(Recombination& r);
    
// private:
    
    float tree_length = 0.0f;
    
    float log_exp(float lambda, float x);
    
    int depth(Node_ptr n);
    
    Node_ptr LCA(Node_ptr n1, Node_ptr n2);
    
    int distance(Node_ptr n1, Node_ptr n2);
    
    float random();
    
};

#endif /* Tree_hpp */
