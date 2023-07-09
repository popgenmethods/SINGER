//
//  TSP_smc.hpp
//  SINGER
//
//  Created by Yun Deng on 4/4/23.
//

#ifndef TSP_smc_hpp
#define TSP_smc_hpp

#include <stdio.h>
#include "random_utils.hpp"
#include "Interval.hpp"
#include "Emission.hpp"

class TSP_smc {
    
public:
    
    float cut_time = 0;
    float lower_bound = 0;
    float gap = 0;
    int min_num = 1;
    float epsilon = 1e-7;
    set<float> check_points = {};
    shared_ptr<Emission> eh;
    
    TSP_smc();
    
    ~TSP_smc();
    
    void set_gap(float q);
    
    void set_emission(shared_ptr<Emission> e);
    
    void set_check_points(set<float> &p);
    
    void reserve_memory(int length);
    
    void start(Branch &branch, float t);
    
    void set_interval_constraint(Recombination &r);
    
    void set_point_constraint(Recombination &r);
    
    Interval *search_point_interval(Recombination &r);
    
    void transfer(Recombination &r, Branch &prev_branch, Branch &next_branch);
    
    void recombine(Branch &prev_branch, Branch &next_branch);
    
    void forward(float rho);
    
    void null_emit(float theta, Node_ptr query_node);
    
    void mut_emit(float theta, float bin_size, set<float> &mut_set, Node_ptr query_node);
    
    map<float, Node_ptr > sample_joining_nodes(int start_index, vector<float> &coordinates);
    
// private:

    int curr_index = 0;
    Branch curr_branch = Branch();
    vector<Interval *> curr_intervals = {};
    map<int, vector<Interval *>>  state_spaces = {{INT_MAX, {}}};
    map<Interval *, Interval *> source_interval = {};
    
    vector<float> rhos = {}; // length: number of blocks - 1
    vector<float> thetas = {}; // length: number of blocks
    
    vector<float> lower_sums = {};
    vector<float> upper_sums = {};
    vector<float> diagonals = {};
    vector<float> factors = {};
    vector<float> lower_diagonals = {};
    vector<float> upper_diagonals = {};
    
    float prev_rho = -1;
    float prev_theta = -1;
    Node_ptr prev_node = nullptr;
    
    int dim = 0;
    vector<float> temp = {};
    vector<float> null_emit_probs = {};
    vector<float> mut_emit_probs = {};
    int sample_index = -1;
    vector<float> trace_back_probs = {};
    vector<vector<float>> forward_probs = {};
    vector<float> emissions = vector<float>(4);
    
    float recomb_cdf(float s, float t);
    
    float recomb_quantile(float s, float q, float lb, float ub);
    
    float recomb_prob(float s, float t1, float t2);
    
    float standard_recomb_cdf(float rho, float s, float t);
    
    float psmc_cdf(float rho, float s, float t);
    
    float psmc_prob(float rho, float s, float t1, float t2);
    
    float get_exp_quantile(float p);
    
    vector<float> generate_grid(float lb, float ub);
    
    float random();
    
    float get_prop(float lb1, float ub1, float lb2, float ub2);
    
    void set_dimensions();
    
    void compute_null_emit_probs(float theta, Node_ptr query_node);
    
    void compute_mut_emit_probs(float theta, float bin_size, set<float> &mut_set, Node_ptr query_node);
    
    void compute_diagonals(float rho);
    
    void compute_lower_diagonals(float rho);
    
    void compute_upper_diagonals(float rho);
    
    void compute_lower_sums();
    
    void compute_upper_sums();
    
    void compute_factors();
    
    void compute_emissions(set<float> &mut_set, Branch branch, Node_ptr node);
    
    void compute_trace_back_probs(float rho, Interval *interval, vector<Interval *> &intervals);
    
    void sanity_check(Recombination &r);
    
    void transfer_intervals(Recombination &r, Branch &prev_branch, Branch &next_branch);
    
    void generate_intervals(Branch &next_branch, float lb, float ub);
    
    vector<Interval *> get_state_space(int x);
    
    int get_interval_index(Interval *interval, vector<Interval *> &intervals);
    
    int get_prev_breakpoint(int x);
    
    Interval *sample_curr_interval(int x);
    
    Interval *sample_prev_interval(Interval *interval, int x);
    
    Interval *sample_source_interval(Interval *interval, int x);
    
    Interval *sample_recomb_interval(Interval *interval, int x);
    
    int trace_back_helper(Interval *interval, int x);
    
    float sample_time(float lb, float ub);
    
    float sample_time(float lb, float ub, float t);
    
    Node_ptr sample_joining_node(Interval *interval);
    
    Node_ptr sample_joining_node(Interval *interval, Node_ptr n);
    
    void write_mean_time(string filename);
};

#endif /* TSP_smc_hpp */
