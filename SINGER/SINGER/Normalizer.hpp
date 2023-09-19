//
//  Normalizer.hpp
//  SINGER
//
//  Created by Yun Deng on 7/13/23.
//

#ifndef Normalizer_hpp
#define Normalizer_hpp

#include <stdio.h>
#include "random_utils.hpp"
#include "ARG.hpp"

class Normalizer {
    
public:
    
    int num_windows = 100;
    float ls = 0;
    float max_time = 20;
    float cap = 10;
    
    vector<Node_ptr> all_root_nodes = {};
    vector<float> all_root_spans = {};
    vector<Node_ptr> all_nodes = {};
    vector<float> all_spans = {};
    set<float> mutation_ages = {};
    vector<float> old_grid = {};
    vector<float> new_grid = {};
    vector<float> observed_mutation_counts = {};
    vector<float> expected_mutation_counts = {};
    vector<float> observed_recombination_counts = {};
    vector<float> recombination_density = {};
    
    Normalizer();
    
    ~Normalizer();
    
    void get_root_span(ARG &a);
    
    void get_node_span(ARG &a);
    
    void randomize_mutation_ages(ARG &a);
    
    void randomized_normalize(ARG &a);
    
    void normalize(ARG &a, float theta);
    
    void partition_arg(ARG &a);
    
    void normalize_recombinations(ARG &a);
    
    void count_mutations(ARG &a);
    
    void count_recombinations(ARG &a);
    
    void add_mutation(float lb, float ub);
    
    void add_recombination(float lb, float ub);
    
    float sample_recombination_time(float lb, float ub);
    
    void sample_recombinations(ARG &a);
};

#endif /* Normalizer_hpp */
