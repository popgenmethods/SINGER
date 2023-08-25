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
    
    // set<Branch> branches = {};
    // map<Node_ptr, Node_ptr, compare_node> parents = {};
    map<Node_ptr, Node_ptr, compare_node> parents = {};
    unordered_map<Node_ptr, unordered_set<Node_ptr>> children = {};
    
    float root_time = 0;
    
    Tree();
    
    float length();
    
    void insert_branch(const Branch &b);
    
    void delete_branch(const Branch &b);
    
    void internal_insert_branch(const Branch &b, float cut_time);
    
    void internal_delete_branch(const Branch &b, float cut_time);
    
    void forward_update(Recombination &r);
    
    void backward_update(Recombination &r);
    
    void remove(Branch b, Node_ptr n);
    
    void add(Branch added_branch, Branch joining_branch, Node_ptr n);
    
    Node_ptr find_sibling(Node_ptr n);
    
    Branch find_joining_branch(Branch removed_branch);
    
    pair<Branch, float> sample_cut_point();
    
    void internal_cut(float cut_time);
    
    void internal_forward_update(Recombination &r, float cut_time);
    
    void internal_backward_update(Recombination &r, float cut_time);
    
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
    
    void impute_states(float m, set<Branch> &mutation_branches);
    
    void impute_states_helper(Node_ptr n, map<Node_ptr, float> &states);
    
    float random();
    
};

#endif /* Tree_hpp */
