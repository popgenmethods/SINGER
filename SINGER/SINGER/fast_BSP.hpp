//
//  fast_BSP.hpp
//  SINGER
//
//  Created by Yun Deng on 4/21/23.
//

#ifndef fast_BSP_hpp
#define fast_BSP_hpp

#include <stdio.h>
#include "Tree.hpp"
#include "Emission.hpp"
#include "Interval.hpp"
#include "Coalescent_calculator.hpp"

class fast_BSP {
    
public:
    
    // basic setup
    float cut_time = 0.0;
    float cutoff = 0;
    shared_ptr<Emission> eh;
    set<float> check_points = {};
    
    // hmm running results
    vector<float> rhos = {};
    vector<float> recomb_sums = {}; // length: number of blocks - 1
    vector<float> weight_sums = {}; // length: number of blocks
    
    // hmm states
    int curr_index = 0;
    map<int, vector<Interval *>>  state_spaces = {{INT_MAX, {}}};
    vector<Interval *> curr_intervals = {};
    
    // coalescent computation
    shared_ptr<Coalescent_calculator> cc;
    
    // transfer at recombinations
    map<Interval *, vector<Interval *>> source_intervals = {};
    map<Interval *, vector<float>> source_weights = {};
    map<Interval_info, vector<Interval *>> temp_intervals = {};
    map<Interval_info, vector<float>> temp_weights = {};
    
    // cache:
    float prev_rho = -1;
    float prev_theta = -1;
    Node *prev_node = nullptr;
    
    // vector computation:
    int dim = 0;
    float recomb_sum = 0;
    float weight_sum = 0;
    vector<float> temp = {};
    vector<float> recomb_probs = {};
    vector<float> recomb_weights = {};
    vector<float> null_emit_probs = {};
    vector<float> mut_emit_probs = {};
    int sample_index = -1;
    vector<float> trace_back_probs = {};
    vector<vector<float>> forward_probs = {};
    
    // states after pruning:
    bool states_change = false;
    set<Branch> valid_branches = {};
    
    
    
    
};

#endif /* fast_BSP_hpp */
