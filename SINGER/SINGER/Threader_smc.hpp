//
//  Threader_smc.hpp
//  SINGER
//
//  Created by Yun Deng on 4/4/23.
//

#ifndef Threader_smc_hpp
#define Threader_smc_hpp

#include <stdio.h>
#include <chrono>
#include <sstream>
#include "ARG.hpp"
#include "BSP.hpp"
#include "succint_BSP.hpp"
#include "BSP_smc.hpp"
#include "fast_BSP.hpp"
#include "fast_BSP_smc.hpp"
#include "reduced_BSP.hpp"
#include "sub_BSP.hpp"
#include "fast_BSP.hpp"
#include "approx_BSP.hpp"
#include "TSP_smc.hpp"
#include "TSP.hpp"
#include "Trace_pruner.hpp"

class Threader_smc {
    
public:
    
    // Threader_smc(double c, double q, shared_ptr<Emission> e);
    
    Threader_smc(double c, double q);
    
    ~Threader_smc();
    
    void thread(ARG &a, Node_ptr n);
    
    void internal_rethread(ARG &a, tuple<double, Branch, double> cut_point);
    
    void terminal_rethread(ARG &a, tuple<double, Branch, double> cut_point);
    
    void fast_thread(ARG &a, Node_ptr n);
    
    void fast_internal_rethread(ARG &a, tuple<double, Branch, double> cut_point);
    
    void fast_terminal_rethread(ARG &a, tuple<double, Branch, double> cut_point);
    
// private:
    
    double cut_time = 0;
    double start = 0;
    double end = 0;
    int start_index = 0;
    int end_index = 0;
    Trace_pruner pruner = Trace_pruner();
    approx_BSP bsp = approx_BSP();
    fast_BSP fbsp = fast_BSP();
    TSP tsp = TSP();
    double gap;
    double cutoff;
    shared_ptr<Binary_emission> be = make_shared<Binary_emission>();
    shared_ptr<Polar_emission> pe = make_shared<Polar_emission>();
    map<double, Branch> new_joining_branches = {};
    map<double, Branch> added_branches = {};
    
    void get_boundary(ARG &a);
    
    void set_check_points(ARG &a);
    
    void run_pruner(ARG &a);
    
    void run_BSP(ARG &a);
    
    void run_fast_BSP(ARG &a);
    
    void run_TSP(ARG &a);
    
    void sample_joining_branches(ARG &a);
    
    void sample_fast_joining_branches(ARG &a);
    
    void sample_joining_points(ARG &a);
    
    double acceptance_ratio(ARG &a);
    
    double random();
    
    vector<double> expected_diff(double m);
    
    vector<double> observed_diff(ARG &a);
    
};

#endif /* Threader_smc_hpp */
