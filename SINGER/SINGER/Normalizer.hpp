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
    
    // map<Node_ptr, float, compare_node> root_span = {};
    vector<Node_ptr> all_root_nodes = {};
    vector<float> all_root_spans = {};
    vector<Node_ptr> all_nodes = {};
    vector<float> all_spans = {};
    set<float> mutation_ages = {};
    map<float, float> mutation_counts = {};
    
    Normalizer();
    
    ~Normalizer();
    
    void get_root_span(ARG &a);
    
    void get_node_span(ARG &a);
    
    void count_mutations(ARG &a);
    
    void randomize_mutation_ages(ARG &a);
    
    void randomized_normalize(ARG &a);
    
    // void normalize(ARG &a);
    
    // void add_mutation(float lb, float ub);
};

#endif /* Normalizer_hpp */
