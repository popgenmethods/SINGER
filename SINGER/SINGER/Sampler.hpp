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
    
    Sampler(float pop_size, float r, float m);
    
    string get_time();
    
    void set_pop_size(float n);
    
    void set_precision(float c, float q);
    
    void set_input_file_prefix(string f);
    
    void set_output_file_prefix(string f);
    
    void set_log_file_prefix(string f);
    
    void set_sequence_length(float x);
    
    void set_num_samples(int n);
    
    Node *build_node(int index, float time);
    
    void build_singleton_arg();
    
    void iterative_start();
    
    void fast_iterative_start();
    
    void terminal_sample(int num_iters);
    
    void sample(int num_iters, int spacing);
    
};

#endif /* Sampler_hpp */
