//
//  fast_coalescent_calculator.cpp
//  SINGER
//
//  Created by Yun Deng on 6/15/23.
//

#include "fast_coalescent_calculator.hpp"

fast_coalescent_calculator::fast_coalescent_calculator(double t) {
    cut_time = t;
    coalescence_times.insert(cut_time);
    coalescence_times.insert(numeric_limits<double>::infinity());
}

fast_coalescent_calculator::~fast_coalescent_calculator() {}

void fast_coalescent_calculator::start(set<Branch> &branches) {
    for (const Branch &b : branches) {
        if (b.lower_node->time > cut_time) {
            coalescence_times.insert(b.lower_node->time);
        }
    }
    compute_first_moment();
}

void fast_coalescent_calculator::update(Recombination &r) {
    double t_old = r.deleted_node->time;
    double t_new = r.inserted_node->time;
    if (t_old > cut_time) {
        coalescence_times.erase(t_old);
        first_moment = 0;
    }
    if (t_new > cut_time) {
        coalescence_times.insert(t_new);
        first_moment = 0;
    }
    if (first_moment == 0) {
        compute_first_moment();
    }
}

void fast_coalescent_calculator::compute_first_moment() {
    first_moment = 0;
    double integral = 0;
    double num_lineages = coalescence_times.size() - 1;
    auto it = coalescence_times.begin();
    double prev_time = cut_time;
    double next_time = cut_time;
    double prev_prob, next_prob;
    double increment = 0;
    while (next(it) != coalescence_times.end()) {
        prev_prob = exp(-integral);
        prev_time = *it;
        next_time = *next(it);
        integral += num_lineages*(next_time - prev_time);
        next_prob = exp(-integral);
        increment = (prev_prob - next_prob)/num_lineages;
        first_moment += increment;
        num_lineages -= 1;
        it++;
    }
}

pair<double, double> fast_coalescent_calculator::compute_time_weights(double x, double y) {
    if (x == y) { // no need to compute when point mass
        return {x, 0};
    }
    double integral = get_integral(x);
    auto it = coalescence_times.upper_bound(x);
    it--;
    // auto u_it = coalescence_times.upper_bound(y);
    double prev_time, next_time;
    double p = 0, q = 0;
    double prev_prob = exp(-integral), next_prob = prev_prob;
    double num_lineages = get_num_lineages(x);
    while (*it < y) {
        prev_prob = exp(-integral);
        prev_time = max(*it, x);
        next_time = min(*next(it), y);
        integral += num_lineages*(next_time - prev_time);
        next_prob = exp(-integral);
        p += (prev_prob - next_prob)/num_lineages;
        if (!isinf(next_time)) {
            q += ((prev_time - cut_time)*prev_prob - (next_time - cut_time)*next_prob)/num_lineages + (prev_prob - next_prob)/num_lineages/num_lineages;
        } else {
            q += (prev_time - cut_time)*prev_prob/num_lineages + (prev_prob - next_prob)/num_lineages/num_lineages;
        }
        num_lineages -= 1;
        it++;
    }
    double t = q/p + cut_time;
    double w = q/first_moment;
    if (y - x < 0.001) {
        t = 0.5*(x + y);
        w = (t - cut_time)*p/first_moment;
    }
    assert(w >= 0);
    assert(t >= x and t <= y);
    return {t, w};
}

double fast_coalescent_calculator::prob(double x, double y) {
    double integral = get_integral(x);
    auto it = coalescence_times.upper_bound(x);
    it--;
    // auto u_it = coalescence_times.upper_bound(y);
    double prev_time = *it, next_time;
    double p = 0;
    double prev_prob = exp(-integral), next_prob = prev_prob;
    double num_lineages = get_num_lineages(x);
    while (*it < y) {
        prev_prob = exp(-integral);
        prev_time = max(*it, x);
        next_time = min(*next(it), y);
        integral += num_lineages*(next_time - prev_time);
        next_prob = exp(-integral);
        p += (prev_prob - next_prob)/num_lineages;
        num_lineages -= 1;
        it++;
    }
    return p;
}

double fast_coalescent_calculator::get_num_lineages(double x) {
    auto u_it = coalescence_times.upper_bound(x);
    // double d0 = distance(u_it, coalescence_times.end());
    double d = coalescence_times.size() - distance(coalescence_times.begin(), u_it);
    // assert(d == d0);
    return d;
}

double fast_coalescent_calculator::get_integral(double x) {
    double integral = 0;
    double num_lineages = coalescence_times.size() - 1;
    auto it = coalescence_times.begin();
    double prev_time = cut_time;
    double next_time = cut_time;
    while (*it < x) {
        prev_time = *it;
        next_time = min(*next(it), x);
        integral += num_lineages*(next_time - prev_time);
        num_lineages -= 1;
        it++;
    }
    return integral;
}
