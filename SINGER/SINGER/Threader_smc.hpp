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
#include "BSP_smc.hpp"
#include "TSP_smc.hpp"

class Threader_smc {
    
public:
    
    Threader_smc(float c, float q, shared_ptr<Emission> e);
    
    ~Threader_smc();
    
    void thread_TSP(ARG &a, Node *n, map<int, Branch> joining_branches);
    
    void thread(ARG &a, Node *n);
    
    void internal_rethread(ARG &a, tuple<int, Branch, float> cut_point);
    
    void terminal_rethread(ARG &old_arg, tuple<int, Branch, float> cut_point);
    
// private:
    
    float cut_time = 0;
    float start = 0;
    float end = 0;
    int start_index = 0;
    int end_index = 0;
    BSP_smc bsp = BSP_smc();
    TSP_smc tsp = TSP_smc();
    float gap;
    float cutoff;
    shared_ptr<Emission> eh;
    std::chrono::system_clock::time_point time = std::chrono::system_clock::now();
    map<float, Branch> new_joining_branches = {};
    map<float, Branch> added_branches = {};
    
    string get_time();
    
    void get_boundary(ARG &a);
    
    void set_check_points(ARG &a);
    
    void run_BSP(ARG &a);
    
    void run_TSP(ARG &a);
    
    void sample_joining_points(ARG &a);
    
    float random();
    
};

#endif /* Threader_smc_hpp */
