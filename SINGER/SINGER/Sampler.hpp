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
#include "ARG.hpp"
#include "Parsimony_pruner.hpp"
#include "Threader_smc.hpp"
#include "Binary_emission.hpp"
#include "Emission.hpp"
#include "Normalizer.hpp"

class Sampler {
    
public:
    
    float rho_unit = 4e-3;
    float Ne = 1;
    float mut_rate = 0;
    float recomb_rate = 0;
    string input_prefix = "";
    string output_prefix = "";
    string log_prefix = "";
    float sequence_length = 0;
    int num_samples = 0;
    ARG arg;
    double bsp_c;
    double tsp_q;
    int random_seed = 0;
    shared_ptr<Emission> eh = make_shared<Binary_emission>();
    int sample_index = 0;
    set<Node_ptr, compare_node> sample_nodes = {};
    vector<Node_ptr> ordered_sample_nodes = {};
    unordered_map<float, set<Node_ptr>> carriers = {};
    
    Sampler(float pop_size, float r, float m);
    
    void set_pop_size(float n);
    
    void set_precision(float c, float q);
    
    void set_input_file_prefix(string f);
    
    void set_output_file_prefix(string f);
    
    void set_log_file_prefix(string f);
    
    void set_sequence_length(float x);
    
    void set_num_samples(int n);
    
    void load_vcf(string vcf_file);
    
    void optimal_ordering();
    
    Node_ptr build_node(int index, float time);
    
    void build_all_nodes();
    
    void build_singleton_arg();
    
    void build_void_arg();
    
    void iterative_start();
    
    void fast_iterative_start();
    
    void multiple_iterative_start(int num_iters);
    
    void multiple_fast_iterative_start(int num_iters);
    
    void recombination_climb(int num_iters, int spacing);
    
    void fast_recombination_climb(int num_iters, int spacing);
    
    void fast_mutation_climb(int num_iters, int spacing);
    
    void terminal_sample(int num_iters);
    
    void internal_sample(int num_iters, int spacing);
    
    void fast_terminal_sample(int num_iters);
    
    void fast_internal_sample(int num_iters, int spacing);
    
    void resume_internal_sample(int num_iters, int spacing, int resume_point);
    
    void resume_fast_internal_sample(int num_iters, int spacing, int resume_point, int seed);
    
    void normalize();
};

#endif /* Sampler_hpp */
