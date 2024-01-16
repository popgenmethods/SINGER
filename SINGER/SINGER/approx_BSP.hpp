//
//  approx_BSP.hpp
//  SINGER
//
//  Created by Yun Deng on 8/15/23.
//

#ifndef approx_BSP_hpp
#define approx_BSP_hpp

#include <stdio.h>
#include <fstream>
#include "Tree.hpp"
#include "Coalescent_calculator.hpp"
#include "approx_coalescent_calculator.hpp"
#include "Interval.hpp"
#include "Emission.hpp"
#include "Binary_emission.hpp"

using Interval_ptr = shared_ptr<Interval>;

class approx_BSP {
    
public:
    
    // basic setup
    double cut_time = 0.0;
    double cutoff = 0;
    double epsilon = 1e-30;
    shared_ptr<Emission> eh;
    set<double> check_points = {};
    
    // pruning parameters
    double rho_unit = 0;
    int grace_period = 0;
    double penalty = 1;
    
    // hmm running results
    vector<double> rhos = {};
    vector<double> recomb_sums = {}; // length: number of blocks - 1
    vector<double> weight_sums = {}; // length: number of blocks
    
    // hmm states
    int curr_index = 0;
    map<int, vector<Interval_ptr>>  state_spaces = {{INT_MAX, {}}};
    vector<Interval_ptr> curr_intervals = {};
    vector<Interval_ptr> temp_intervals = {};
    map<int, vector<double>> times = {{INT_MAX, {}}};
    map<int, vector<double>> weights = {{INT_MAX, {}}};
    
    // coalescent computation
    shared_ptr<approx_coalescent_calculator> cc;
    
    // transfer at recombinations
    map<Interval_info, vector<Interval_ptr>> transfer_intervals = {};
    map<Interval_info, vector<double>> transfer_weights = {};
    
    // cache:
    double prev_rho = -1;
    double prev_theta = -1;
    Node_ptr prev_node = nullptr;
    
    // vector computation:
    int dim = 0;
    double recomb_sum = 0;
    double weight_sum = 0;
    vector<double> temp = {};
    vector<double> time_points = {};
    vector<double> raw_weights = {};
    vector<double> recomb_probs = {};
    vector<double> recomb_weights = {};
    vector<double> null_emit_probs = {};
    vector<double> mut_emit_probs = {};
    int sample_index = -1;
    vector<double> trace_back_probs = {};
    vector<vector<double>> forward_probs = {};
    
    // states after pruning:
    bool states_change = false;
    set<Branch> valid_branches = {};
    
    approx_BSP();
    
    ~approx_BSP();
    
    void reserve_memory(int length);
    
    void start(set<Branch> &branches, double t);
    
    void start(Tree &tree, double t);
    
    void set_cutoff(double x);
    
    void set_emission(shared_ptr<Emission> e);
    
    void set_check_points(set<double> &p);
    
    void forward(double rho); // forward pass when there is no recombination (without emission). Also update recomb_sums and weight_sums.
    
    void transfer(Recombination &r); // forward pass when there is a recombination (without emission), and add a transition object. Also update active intervals, recomb_sums and weight_sums.

    double get_recomb_prob(double rho, double t);
    
    void null_emit(double theta, Node_ptr query_node);
    
    void mut_emit(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node);
    
    map<double, Branch> sample_joining_branches(int start_index, vector<double> &coordinates);
    
    void set_dimensions();
    
    void compute_recomb_probs(double rho);
    
    void compute_recomb_weights(double rho);
    
    void compute_null_emit_prob(double theta, Node_ptr query_node);
    
    void compute_mut_emit_probs(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node);
    
    void transfer_helper(Interval_info &next_interval, Interval_ptr &prev_interval, double w);
    
    void transfer_helper(Interval_info &next_interval);
    
    void add_new_branches(Recombination &r);
    
    void compute_interval_info();
    
    void sanity_check(Recombination &r);
    
    void generate_intervals(Recombination &r);
    
    double get_overwrite_prob(Recombination &r, double lb, double ub);
    
    void process_interval(Recombination &r, int i);
    
    void process_source_interval(Recombination &r, int i);
    
    void process_target_interval(Recombination &r, int i);
    
    void process_other_interval(Recombination &r, int i);
    
    double random();
    
    int get_prev_breakpoint(int x);
    
    vector<Interval_ptr> &get_state_space(int x);
    
    vector<double> &get_time_points(int x);
    
    vector<double> &get_raw_weights(int x);
    
    int get_interval_index(Interval_ptr interval, vector<Interval_ptr> &intervals);
    
    void simplify(map<double, Branch> &joining_branches);
    
    Interval_ptr sample_curr_interval(int x);
    
    Interval_ptr sample_prev_interval(int x);
    
    Interval_ptr sample_source_interval(Interval_ptr interval, int x);
    
    int trace_back_helper(Interval_ptr interval, int x);
    
    double avg_num_states();
    
};

#endif /* approx_BSP_hpp */
