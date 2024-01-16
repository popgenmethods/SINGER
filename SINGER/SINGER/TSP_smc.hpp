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
    
    double cut_time = 0;
    double lower_bound = 0;
    double gap = 0;
    int min_num = 1;
    double epsilon = 1e-7;
    set<double> check_points = {};
    shared_ptr<Emission> eh;
    static int counter;
    
    TSP_smc();
    
    ~TSP_smc();
    
    void set_gap(double q);
    
    void set_emission(shared_ptr<Emission> e);
    
    void set_check_points(set<double> &p);
    
    void reserve_memory(int length);
    
    void start(Branch &branch, double t);
    
    void set_interval_constraint(Recombination &r);
    
    void set_point_constraint(Recombination &r);
    
    Interval *search_point_interval(Recombination &r);
    
    void transfer(Recombination &r, Branch &prev_branch, Branch &next_branch);
    
    void recombine(Branch &prev_branch, Branch &next_branch);
    
    void forward(double rho);
    
    void null_emit(double theta, Node_ptr query_node);
    
    void mut_emit(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node);
    
    map<double, Node_ptr > sample_joining_nodes(int start_index, vector<double> &coordinates);
    
// private:

    int curr_index = 0;
    Branch curr_branch = Branch();
    vector<Interval *> curr_intervals = {};
    map<int, vector<Interval *>>  state_spaces = {{INT_MAX, {}}};
    map<Interval *, Interval *> source_interval = {};
    
    vector<double> rhos = {}; // length: number of blocks - 1
    vector<double> thetas = {}; // length: number of blocks
    
    vector<double> lower_sums = {};
    vector<double> upper_sums = {};
    vector<double> diagonals = {};
    vector<double> factors = {};
    vector<double> lower_diagonals = {};
    vector<double> upper_diagonals = {};
    
    double prev_rho = -1;
    double prev_theta = -1;
    Node_ptr prev_node = nullptr;
    
    int dim = 0;
    vector<double> temp = {};
    vector<double> null_emit_probs = {};
    vector<double> mut_emit_probs = {};
    int sample_index = -1;
    vector<double> trace_back_probs = {};
    vector<vector<double>> forward_probs = {};
    vector<double> emissions = vector<double>(4);
    
    double recomb_cdf(double s, double t);
    
    double recomb_quantile(double s, double q, double lb, double ub);
    
    double recomb_prob(double s, double t1, double t2);
    
    double standard_recomb_cdf(double rho, double s, double t);
    
    double psmc_cdf(double rho, double s, double t);
    
    double psmc_prob(double rho, double s, double t1, double t2);
    
    double get_exp_quantile(double p);
    
    vector<double> generate_grid(double lb, double ub);
    
    double random();
    
    double get_prop(double lb1, double ub1, double lb2, double ub2);
    
    void set_dimensions();
    
    void compute_null_emit_probs(double theta, Node_ptr query_node);
    
    void compute_mut_emit_probs(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node);
    
    void compute_diagonals(double rho);
    
    void compute_lower_diagonals(double rho);
    
    void compute_upper_diagonals(double rho);
    
    void compute_lower_sums();
    
    void compute_upper_sums();
    
    void compute_factors();
    
    void compute_emissions(set<double> &mut_set, Branch branch, Node_ptr node);
    
    void compute_trace_back_probs(double rho, Interval *interval, vector<Interval *> &intervals);
    
    void sanity_check(Recombination &r);
    
    void transfer_intervals(Recombination &r, Branch &prev_branch, Branch &next_branch);
    
    void generate_intervals(Branch &next_branch, double lb, double ub);
    
    vector<Interval *> get_state_space(int x);
    
    int get_interval_index(Interval *interval, vector<Interval *> &intervals);
    
    int get_prev_breakpoint(int x);
    
    Interval *sample_curr_interval(int x);
    
    Interval *sample_prev_interval(Interval *interval, int x);
    
    Interval *sample_source_interval(Interval *interval, int x);
    
    Interval *sample_recomb_interval(Interval *interval, int x);
    
    int trace_back_helper(Interval *interval, int x);
    
    double sample_time(double lb, double ub);
    
    double sample_time(double lb, double ub, double t);
    
    double exp_median(double lb, double ub);
    
    Node_ptr sample_joining_node(Interval *interval);
    
    Node_ptr sample_joining_node(Interval *interval, Node_ptr n);
    
    void write_mean_time(string filename);
};

#endif /* TSP_smc_hpp */
