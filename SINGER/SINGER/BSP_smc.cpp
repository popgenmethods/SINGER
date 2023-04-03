//
//  BSP_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 4/2/23.
//

#include "BSP_smc.hpp"

BSP_smc::BSP_smc() {}

BSP_smc::~BSP_smc() {
    for (auto x : state_spaces) {
        for (Interval *i : x.second) {
            delete(i);
        }
    }
}

void BSP_smc::start(set<Branch> branches, float x, float t) {
    cut_time = t;
    start_pos = x;
    end_pos = x;
    curr_pos = x;
    float lb;
    float ub;
    float p;
    Interval *new_interval = nullptr;
    coalescence_times = {cut_time};
    for (Branch b : branches) {
        if (b.upper_node->time > cut_time) {
            coalescence_times.insert(b.upper_node->time);
        }
    }
    calculate_coalescence_stats();
    for (Branch b : branches) {
        if (b.upper_node->time > cut_time) {
            lb = max(b.lower_node->time, cut_time);
            ub = b.upper_node->time;
            p = coalescence_probs.at(ub) - coalescence_probs.at(lb);
            new_interval = new Interval(b, lb, ub, curr_pos, p);
            curr_intervals.push_back(new_interval);
        }
    }
    fill_interval_info();
    state_spaces.insert({curr_pos, curr_intervals});
}

void BSP_smc::set_cutoff(float x) {
    cutoff = x;
}

void BSP_smc::set_check_points(set<float> p) {
    check_points = p;
}

void BSP_smc::forward(float rho) {
    curr_pos += 1;
    end_pos += 1;
    rhos.push_back(rho);
    float w;
    float p;
    float recomb_prob;
    float recomb_sum = 0.0;
    float weight_sum = 0.0;
    for (Interval *i : curr_intervals) {
        recomb_prob = i->get_recomb_prob(rho, cut_time);
        p = i->get_prob()*recomb_prob;
        recomb_sum += p;
        p = i->get_prob() - p;
        i->new_prob(p);
        if (i->full_branch(cut_time)) {
            weight_sum += recomb_prob*i->weight;
        }
    }
    for (Interval *i : curr_intervals) {
        if (i->full_branch(cut_time)) {
            w = i->weight*i->get_recomb_prob(rho, cut_time);
            p = w*recomb_sum/weight_sum;
            i->add_prob(p);
        }
    }
    recomb_sums.push_back(recomb_sum);
    weight_sums.push_back(weight_sum);
}


void BSP_smc::transfer(Recombination r) {
    rhos.push_back(0);
    recomb_sums.push_back(0);
    curr_pos += 1;
    end_pos += 1;
    sanity_check(r);
    vector<Interval *> intervals = curr_intervals;
    curr_intervals.clear();
    for (Interval *i : intervals) {
        process_interval(r, i);
    }
    update_coalescence_times(r);
    calculate_coalescence_stats();
    // add_new_branches(r);
    // get_new_intervals(r);
    state_spaces.insert({curr_pos, curr_intervals});
}
