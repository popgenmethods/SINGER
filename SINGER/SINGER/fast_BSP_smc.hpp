//
//  fast_BSP_smc.hpp
//  SINGER
//
//  Created by Yun Deng on 4/21/23.
//

#ifndef fast_BSP_smc_hpp
#define fast_BSP_smc_hpp

#include <stdio.h>
#include "Tree.hpp"
#include "Emission.hpp"
#include "Interval.hpp"
#include "Coalescent_calculator.hpp"

class fast_BSP_smc {
    
public:
    
    // basic setup
    double cut_time = 0.0;
    double cutoff = 0.0;
    double epsilon = 1e-6;
    shared_ptr<Emission> eh;
    set<double> check_points = {};
    equal_interval ei = equal_interval();
    
    // hmm running results
    vector<double> rhos = {};
    vector<double> recomb_sums = {}; // length: number of blocks - 1
    vector<double> weight_sums = {}; // length: number of blocks
    
    // hmm states
    int curr_index = 0;
    map<int, vector<Interval *>>  state_spaces = {{INT_MAX, {}}};
    vector<Interval *> curr_intervals = {};
    vector<double> temp_probs = {};
    vector<Interval *> temp_intervals = {};
    
    // coalescent computation
    shared_ptr<Coalescent_calculator> cc;
    
    // transfer at recombinations
    map<Interval_info, vector<Interval *>> transfer_intervals = {};
    map<Interval_info, vector<double>> transfer_weights = {};
    
    // cache:
    double prev_rho = -1;
    double prev_theta = -1;
    Node_ptr prev_node = nullptr;
    
    // vector computation:
    int dim = 0;
    double recomb_sum = 0;
    double weight_sum = 0;
    vector<double> recomb_probs = {};
    vector<double> recomb_weights = {};
    vector<double> null_emit_probs = {};
    vector<double> mut_emit_probs = {};
    int sample_index = -1;
    vector<double> trace_back_probs = {};
    vector<vector<double>> forward_probs = {};
    
    // states after pruning:
    bool states_change = false;
    set<Branch> covered_branches = {};
    set<Branch> valid_branches = {};
    
    fast_BSP_smc();
    
    ~fast_BSP_smc();
    
    void reserve_memory(int length);
    
    void start(set<Branch> &branches, double t);
    
    void set_cutoff(double x);
    
    void set_emission(shared_ptr<Emission> e);
    
    void set_check_points(set<double> &p);
    
    void forward(double rho); // forward pass when there is no recombination (without emission). Don't forget to update recomb_sums and weight_sums.
    
    void transfer(Recombination &r); // forward pass when there is a recombination (without emission), and add a transition object. Don't forget to update active intervals, recomb_sums and weight_sums.
    
    void regular_forward(double rho);
    
    void update(double rho);
    
    double get_recomb_prob(double rho, double t);
    
    void null_emit(double theta, Node_ptr query_node);
    
    void mut_emit(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node);
    
    map<double, Branch> sample_joining_branches(int start_index, vector<double> &coordinates);
    
    void write_forward_probs(string filename);
    
    // private methods:
    
    void update_states(set<Branch> &deletions, set<Branch> &insertions);
    
    void set_states(set<Branch> &branches);
    
    void set_dimensions();
    
    void compute_recomb_probs(double rho);
    
    void compute_recomb_weights(double rho);
    
    void compute_null_emit_prob(double theta, Node_ptr query_node);
    
    void compute_mut_emit_probs(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node);
    
    void transfer_helper(Interval_info next_interval, Interval *prev_interval, double w);
    
    void transfer_helper(Interval_info next_interval);
    
    double compute_transfer_prob();
    
    Interval *duplicate_interval(Interval *interval);
    
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
    
    vector<Interval *> &get_state_space(int x);
    
    int get_interval_index(Interval *interval, vector<Interval *> &intervals);
    
    void simplify(map<double, Branch> &joining_branches);
    
    Interval *sample_curr_interval(int x);
    
    Interval *sample_prev_interval(int x);
    
    Interval *sample_source_interval(Interval *interval, int x);
    
    Interval *sample_connection_interval(Interval *interval, int x);
    
    int trace_back_helper(Interval *interval, int x);
    
    void check_intervals();
    
    void check_recomb_sums();
    
};

#endif /* fast_BSP_smc_hpp */
