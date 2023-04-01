//
//  Parsimony_pruner.hpp
//  SINGER
//
//  Created by Yun Deng on 3/14/23.
//

#ifndef Parsimony_pruner_hpp
#define Parsimony_pruner_hpp

#include <stdio.h>
#include <algorithm>
#include "Pruner.hpp"
#include "ARG.hpp"

class Parsimony_pruner {
    
public:
    
    float max_mismatch = 0.999;
    float max_time = 100;
    float start = 0;
    float end = 0;
    map<float, Node *> nodes = {};
    set<float> private_mutations = {};
    map<Branch, float> curr_mismatch = {};
    map<float, float> match_map = {};
    map<float, float> potential_seeds = {};
    set<float> used_seeds = {};
    
    map<float, set<Branch>> reductions = {};
    map<float, set<Branch>> deleted_branches = {};
    map<float, set<Branch>> inserted_branches = {};
    
    map<Branch, set<Branch>> transitions = {};

    Parsimony_pruner();
    
    void prune_arg(ARG &a, map<float, Node *> base_nodes);
    
    void start_search(ARG &a, float m);
    
    void extend(ARG &a, float x);

    void mutation_forward(Node *n, float m);
    
    void mutation_backward(Node *n, float m);

    void recombination_forward(Recombination &r);
    
    void recombination_backward(Recombination &r);
    
    void check_reduction(map<float, pair<Branch, Node *>> joining_points);
    
    void write_reduction_size(string filename);

    // private:
    
    Node *get_node_at(float x);
    
    void build_match_map(ARG &a, map<float, Node *> base_nodes);
    
    float find_minimum_match();
    
    float count_mismatch(Branch branch, Node *n, float m);
    
    void transition_helper(Branch sb, Branch tb);
    
    void update_mismatch();
    
    void write_init_set();
    
    void write_reduction_change(float x, set<Branch> db, set<Branch> ib);
    
    void write_reductions(ARG &a);
    
    void extend_forward(ARG &a, float x);
    
    void extend_backward(ARG &a, float x);
    
};

#endif /* Parsimony_pruner_hpp */
