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
    
    map<Node_ptr, Node_ptr, compare_node> parents = {};
    unordered_map<Node_ptr, unordered_set<Node_ptr>> children = {};
    
    Tree();
    
    double length();
    
    void insert_branch(const Branch &b);
    
    void delete_branch(const Branch &b);
    
    void internal_insert_branch(const Branch &b, double cut_time);
    
    void internal_delete_branch(const Branch &b, double cut_time);
    
    void forward_update(Recombination &r);
    
    void backward_update(Recombination &r);
    
    void remove(Branch b, Node_ptr n);
    
    void add(Branch added_branch, Branch joining_branch, Node_ptr n);
    
    Node_ptr find_sibling(Node_ptr n);
    
    Branch find_joining_branch(Branch removed_branch);
    
    pair<Branch, double> sample_cut_point();
    
    void internal_cut(double cut_time);
    
    void internal_forward_update(Recombination &r, double cut_time);
    
    void internal_backward_update(Recombination &r, double cut_time);
    
    double prior_likelihood();
    
    double data_likelihood(double theta, double pos);
    
    double null_likelihood(double theta);
    
    double data_likelihood(double theta, double bin_size, set<double> mutations);
    
    double transition_likelihood(Recombination& r);
    
// private:
    
    double tree_length = 0.0f;
    
    double log_exp(double lambda, double x);
    
    int depth(Node_ptr n);
    
    Node_ptr LCA(Node_ptr n1, Node_ptr n2);
    
    int distance(Node_ptr n1, Node_ptr n2);
    
    void impute_states(double m, set<Branch> &mutation_branches);
    
    void impute_states_helper(Node_ptr n, map<Node_ptr, double> &states);
    
    double random();
    
};

#endif /* Tree_hpp */
