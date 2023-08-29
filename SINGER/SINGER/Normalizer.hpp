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
    
    vector<Node_ptr> all_root_nodes = {};
    vector<float> all_root_spans = {};
    vector<Node_ptr> all_nodes = {};
    vector<float> all_spans = {};
    set<float> mutation_ages = {};
    vector<float> old_grid = {};
    vector<float> new_grid = {};
    vector<float> observed_mutation_counts = {};
    vector<float> expected_mutation_counts = {};
    vector<float> observed_branch_length = {};
    vector<float> expected_branch_length = {};
    
    Normalizer();
    
    ~Normalizer();
    
    void get_root_span(ARG &a);
    
    void get_node_span(ARG &a);
    
    void randomize_mutation_ages(ARG &a);
    
    void randomized_normalize(ARG &a);
    
    void normalize(ARG &a, float theta);
    
    void partition_arg(ARG &a);
    
    void count_mutations(ARG &a);
    
    void calculate_branch_length(ARG &a, float Ne);
    
    void add_mutation(float lb, float ub);
};

#endif /* Normalizer_hpp */
