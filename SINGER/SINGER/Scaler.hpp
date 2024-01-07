//
//  Scaler.hpp
//  SINGER
//
//  Created by Yun Deng on 1/3/24.
//

#ifndef Scaler_hpp
#define Scaler_hpp

#include <stdio.h>
#include "random_utils.hpp"
#include "ARG.hpp"

class Scaler {
    
public:
    
    int num_windows = 100;
    vector<Node_ptr> sorted_nodes = {};
    vector<double> node_deltas = {};
    vector<double> rates = {};
    vector<double> accumulated_arg_length = {};
    vector<double> old_grid = {0};
    vector<double> new_grid = {0};
    vector<double> expected_arg_length = {};
    vector<double> observed_arg_length = {};
    vector<double> scaling_factors = {};
    
    Scaler();
    
    void compute_deltas(ARG &a);
    
    void compute_old_grid();
    
    void compute_new_grid(double theta);
    
    void map_mutations(ARG &a);
    
    void add_mutation(double w, double lb, double ub);
    
    void rescale(ARG &a, double theta);
    
};

#endif /* Scaler_hpp */
