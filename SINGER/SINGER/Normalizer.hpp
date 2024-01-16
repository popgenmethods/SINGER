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
    double ls = 0;
    double max_time = 20;
    
    vector<Node_ptr> all_root_nodes = {};
    vector<double> all_root_spans = {};
    vector<Node_ptr> all_nodes = {};
    vector<double> all_spans = {};
    set<double> mutation_ages = {};
    vector<double> old_grid = {};
    vector<double> new_grid = {};
    vector<double> observed_mutation_counts = {};
    vector<double> expected_mutation_counts = {};
    vector<double> observed_recombination_counts = {};
    vector<double> recombination_density = {};
    
    Normalizer();
    
    ~Normalizer();
    
    void compute_max_time(ARG &a);
    
    void get_root_span(ARG &a);
    
    void get_node_span(ARG &a);
    
    void randomize_mutation_ages(ARG &a);
    
    void randomized_normalize(ARG &a);
    
    void normalize(ARG &a, double theta);
    
    void partition_arg(ARG &a);
    
    void normalize_recombinations(ARG &a);
    
    void count_mutations(ARG &a);
    
    void count_recombinations(ARG &a);
    
    void add_mutation(double lb, double ub);
    
    void add_recombination(double lb, double ub);
    
    double sample_recombination_time(double lb, double ub);
    
    void sample_recombinations(ARG &a);
};

#endif /* Normalizer_hpp */
