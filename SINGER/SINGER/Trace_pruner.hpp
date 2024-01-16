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
    
    double cutoff = 1e-2;
    double mut_prob = 4e-4;
    double max_time = 100;
    double start = 0;
    double end = 0;
    double cut_time = 0;
    int band_width = 10;
    
    double length = 0;
    
    map<double, Branch> queries = {};
    set<double> private_mutations = {};
    
    map<double, Tree> seed_trees = {};
    
    map<double, double> match_map = {};
    map<double, double> potential_seeds = {};
    set<double> used_seeds = {};
    
    map<Branch, double> seed_match = {};
    map<Interval_info, double> seed_scores = {};
    map<Interval_info, double> curr_scores = {};
    
    set<double> check_points;
    
    map<double, set<Branch>> reductions = {};
    map<double, set<Interval_info>> deletions = {};
    map<double, set<Interval_info>> insertions = {};
    
    map<Interval_info, double> transition_scores = {};
    
    set<pair<double, double>> segments = {};
    
    Trace_pruner();
    
    void prune_arg(ARG &a);
    
    void set_check_points(set<double> &p);
    
    void start_search(ARG &a, double m);
    
    void extend(ARG &a, double x);
    
    void extend_forward(ARG &a, double x);
    
    void extend_backward(ARG &a, double x);
    
    void mutation_update(Node_ptr n, double m);

    void recombination_forward(Recombination &r);
    
    void recombination_backward(Recombination &r);
    
    void write_reduction_distance(ARG &a, string filename);
    
    void write_reduction_size(string filename);
    
    double min_reduction_error();
    
    void write_reductions(ARG &a);
    
    Node_ptr get_node_at(double x);
    
    double get_match_time(set<Branch> &branches, double m, Node_ptr n);
    
    void build_match_map(ARG &a);
    
    double find_closest_reference(double x);
    
    double find_minimum_match();
    
    double count_mismatch(Branch branch, Node_ptr n, double m);
    
    void forward_prune_states(double x);
    
    void backward_prune_states(double x);
    
    void delete_all(double x);
    
    void insert_all(double x);
    
    void forward_transition(Recombination &r, const Interval_info &interval);
    
    void backward_transition(Recombination &r, const Interval_info &interval);
    
    void forward_transition_helper(Interval_info prev_interval, Interval_info next_interval, double x, double p);
    
    void backward_transition_helper(Interval_info next_interval, Interval_info prev_interval, double x, double p);
    
    double exp_prob(double l, double u);
    
    double exp_prop(double l, double u, double x, double y);
    
    double exp_median(double l, double u);
    
    double forward_overwrite_prob(Recombination &r, double lb, double ub);
    
    double backward_overwrite_prob(Recombination &r, double lb, double ub);
    
    void remove_segment(double x, double y);
    
    void restrict_search();
};

#endif /* Trace_pruner_hpp */
