//
//  Trace_pruner.hpp
//  SINGER
//
//  Created by Yun Deng on 4/25/23.
//

#ifndef Trace_pruner_hpp
#define Trace_pruner_hpp

#include <stdio.h>
#include <algorithm>
#include "Pruner.hpp"
#include "Interval.hpp"
#include "ARG.hpp"

class Trace_pruner : public Pruner {
    
public:
    
    float cutoff = 0.001;
    float mut_prob = 4e-4;
    float max_time = 100;
    float start = 0;
    float end = 0;
    float cut_time = 0;
    
    map<float, Branch> queries = {};
    set<float> private_mutations = {};
    
    map<float, Tree> seed_trees = {};
    
    map<float, float> match_map = {};
    map<float, float> potential_seeds = {};
    set<float> used_seeds = {};
    
    map<Branch, float> seed_match = {};
    map<Interval_info, float> seed_scores = {};
    map<Interval_info, float> curr_scores = {};
    
    map<float, set<Branch>> reductions = {};
    map<float, set<Interval_info>> deletions = {};
    map<float, set<Interval_info>> insertions = {};
    
    map<Interval_info, float> transition_scores = {};
    
    Trace_pruner();
    
    void prune_arg(ARG &a);
    
    void start_search(ARG &a, float m);
    
    void extend(ARG &a, float x);
    
    void extend_forward(ARG &a, float x);
    
    void extend_backward(ARG &a, float x);
    
    void mutation_update(Node *n, float m);

    void recombination_forward(Recombination &r);
    
    void recombination_backward(Recombination &r);
    
    void write_reduction_distance(ARG &a, string filename);
    
    void write_reduction_size(string filename);
    
    void write_reductions(ARG &a);
    
    Node *get_node_at(float x);
    
    void build_match_map(ARG &a);
    
    float find_closest_reference(float x);
    
    float find_minimum_match();
    
    float count_mismatch(Branch branch, Node *n, float m);
    
    void forward_prune_states(float x);
    
    void backward_prune_states(float x);
    
    void delete_all(float x);
    
    void insert_all(float x);
    
    void forward_transition(Recombination &r, const Interval_info &interval);
    
    void backward_transition(Recombination &r, const Interval_info &interval);
    
    void forward_transition_helper(Interval_info prev_interval, Interval_info next_interval, float x, float p);
    
    void backward_transition_helper(Interval_info next_interval, Interval_info prev_interval, float x, float p);
    
    float exp_prob(float l, float u);
    
    float exp_prop(float l, float u, float x, float y);
    
    float exp_median(float l, float u);
    
    float forward_overwrite_prob(Recombination &r, float lb, float ub);
    
    float backward_overwrite_prob(Recombination &r, float lb, float ub);
};

#endif /* Trace_pruner_hpp */
