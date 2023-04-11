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
#include "Branch.hpp"
#include "Recombination.hpp"

class Tree {

public:
    
    set<Branch> branches = {};
    map<Node *, Node *> parents = {};
    
    Tree();
    
    float length();
    
    void insert_branch(Branch b);
    
    void delete_branch(Branch b);
    
    void forward_update(Recombination &r);
    
    void backward_update(Recombination &r);
    
    Node *find_sibling(Node *n);
    
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
    
    int depth(Node *n);
    
    Node *LCA(Node *n1, Node *n2);
    
    int distance(Node *n1, Node *n2);
    
    float random();
    
};

#endif /* Tree_hpp */
