//
//  Sampler.hpp
//  SINGER
//
//  Created by Yun Deng on 3/31/23.
//

#ifndef Sampler_hpp
#define Sampler_hpp

#include <stdio.h>
#include <chrono>
#include <sstream>
#include <zlib.h>
#include "ARG.hpp"
#include "Threader_smc.hpp"
#include "Binary_emission.hpp"
#include "Emission.hpp"
#include "Normalizer.hpp"
#include "Scaler.hpp"
#include "Rate_map.hpp"
#include "gzstream.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <numeric>
#include <unordered_set>
#include <random>

class Sampler {
    
public:
    
    double rho_unit = 4e-3;
    double Ne = 1;
    Rate_map recomb_map;
    Rate_map mut_map;
    double mut_rate = 0;
    double recomb_rate = 0;
    string input_prefix = "";
    string output_prefix = "";
    string log_prefix = "";
    double start = 0;
    double end = 0;
    double sequence_length = 0;
    int num_samples = 0;
    ARG arg;
    bool fast_mode = false;
    double bsp_c = 0.01;
    double tsp_q = 0.05;
    int random_seed = 0;
    double penalty = 0.01;
    double polar = 0.99;
    int sample_index = 0;
    set<Node_ptr, compare_node> sample_nodes = {};
    vector<Node_ptr> ordered_sample_nodes = {};
    unordered_map<double, set<Node_ptr>> carriers = {};
    unordered_map<Node_ptr, set<double>> mutation_sets = {};
    
    Sampler();
    
    Sampler(double pop_size, double r, double m);
    
    Sampler(double pop_size, Rate_map &rm, Rate_map &mm);
    
    void set_pop_size(double n);
    
    void set_precision(double c, double q);
    
    void set_input_file_prefix(string f);
    
    void set_output_file_prefix(string f);
    
    void set_log_file_prefix(string f);
    
    void set_sequence_length(double x);
    
    void set_num_samples(int n);
    
    void naive_read_vcf(string prefix, double start_pos, double end_pos);
    
    void guide_read_vcf(string prefix, double start, double end);
    
    void load_vcf(string prefix, double start, double end);
    
    void optimal_ordering();
    
    Node_ptr build_node(int index, double time);
    
    void build_all_nodes();
    
    void build_singleton_arg();
    
    void build_void_arg();
    
    void iterative_start();
    
    void fast_iterative_start();
    
    // void recombination_climb(int num_iters, int spacing);
    
    // void mutation_climb(int num_iters, int spacing);
    
    // void fast_recombination_climb(int num_iters, int spacing);
    
    // void fast_mutation_climb(int num_iters, int spacing);
    
    // void terminal_sample(int num_iters);
    
    void internal_sample(int num_iters, int spacing);
    
    // void fast_terminal_sample(int num_iters);
    
    void fast_internal_sample(int num_iters, int spacing);
    
    void resume_internal_sample(int num_iters, int spacing);
    
    void debug_resume_internal_sample(int num_iters, int spacing);
    
    void resume_fast_internal_sample(int num_iters, int spacing);
    
    void debug_resume_fast_internal_sample(int num_iters, int spacing);
    
    void normalize();
    
    void rescale();
    
    void start_log();
    
    void write_iterative_start();
    
    void write_sample();
    
    void write_cut(tuple<double, Branch, double> cut_point);
    
    void load_resume_arg();
    
    vector<string> read_last_line(string filename);
    
    void read_resume_point(string filename);
    
    void retract_log(int k);
};

#endif /* Sampler_hpp */
