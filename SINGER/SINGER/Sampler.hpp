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
#include "Threader.hpp"

class Sampler {
    
public:
    
    float rho_unit = 4e-3;
    float Ne = 1;
    float mut_rate = 0;
    float recomb_rate = 0;
    // Binary_emission eh = Binary_emission();
    string input_prefix = "";
    string output_prefix = "";
    string log_prefix = "";
    float sequence_length = 0;
    int num_samples = 0;
    ARG arg;
    double bsp_c;
    double tsp_q;
    int random_seed = 0;
    
    Sampler(float pop_size, float r, float m);
    
    string get_time();
    
    void set_pop_size(float n);
    
    void set_precision(float c, float q);
    
    void set_input_file_prefix(string f);
    
    void set_output_file_prefix(string f);
    
    void set_log_file_prefix(string f);
    
    void set_sequence_length(float x);
    
    void set_num_samples(int n);
    
    void build_singleton_arg();
};

#endif /* Sampler_hpp */
