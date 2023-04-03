//
//  Interval.hpp
//  SINGER
//
//  Created by Yun Deng on 4/9/22.
//

#ifndef Interval_hpp
#define Interval_hpp

#include <stdio.h>
#include <map>
#include <numeric>
#include "Recombination.hpp"

using namespace std;

class Interval {
    
private:
    
    vector<float> probs;
    
public:
    
    Branch branch;
    float lb = 0;
    float ub = 0;
    float weight = 0.0;
    float time = 0.0;
    float start_pos = 0;
    vector<float> source_weights = {};
    vector<Interval *> source_intervals = {};
    Node *node = nullptr;
    float reduction = 1.0;
    
    Interval(Branch b, float tl, float tu, float init_pos, float init_prob);
    
    void new_prob(float p);
    
    void add_prob(float p);
    
    float get_prob();
    
    float get_prob_at(int x);
    
    void update_prob(float p);
    
    void assign_weight(float w);
    
    void assign_time(float t);
    
    int last_pos();
    
    float get_recomb_prob(float rho, float base_time);
    
    void fill_time();
    
    void rescale(float a);
    
    void multiply(float a);
    
    bool full_branch(float t);

    void set_source(vector<Interval *> intervals, vector<float> weights);
    
    void set_node(Node *n);
    
    Interval *sample_source();
};

struct compare_interval {
    
    bool operator()(const Interval *i1, const Interval *i2) const {
        if (i1->branch != i2->branch) {
            return i1->branch < i2->branch;
        }
        return i1->time < i2->time;
    }
};

class Interval_info {
    
public:
    
    Branch branch;
    float lb;
    float ub;
    
    Interval_info();
    
    Interval_info(Branch b, float tl, float tu);
    
};

bool operator==(const Interval_info& i1, const Interval_info& i2);

bool operator!=(const Interval_info& i1, const Interval_info& i2);

bool operator<(const Interval_info& i1, const Interval_info& i2);

#endif /* Interval_hpp */
