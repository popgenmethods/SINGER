//
//  Coalescent_calculator.cpp
//  SINGER
//
//  Created by Yun Deng on 4/17/23.
//

#include "Coalescent_calculator.hpp"

Coalescent_calculator::Coalescent_calculator(float t) {
    cut_time = t;
}

Coalescent_calculator::~Coalescent_calculator() {}

void Coalescent_calculator::start(set<Branch> &inserted_branches) {
    for (Branch b : inserted_branches) {
        branches.insert(b);
    }
    compute();
}

void Coalescent_calculator::update(set<Branch> &deleted_branches, set<Branch> &inserted_branches) {
    for (Branch b : deleted_branches) {
        assert(branches.count(b) > 0);
        branches.erase(b);
    }
    for (Branch b : inserted_branches) {
        branches.insert(b);
    }
    compute();
}

void Coalescent_calculator::compute() {
    compute_rate_changes();
    compute_rates();
    compute_probs_quantiles();
}

float Coalescent_calculator::weight(float lb, float ub) {
    float p = prob(ub) - prob(lb);
    return p;
}

float Coalescent_calculator::time(float lb, float ub) {
    float lq = prob(lb);
    float uq = prob(ub);
    float t;
    if (isinf(ub)) {
        return lb + log(2);
    }
    if (ub - lb < 1e-3 or uq - lq < 1e-3) {
        t = 0.5*(lb + ub);
    } else {
        t = quantile(0.5*(lq + uq));
    }
    return t;
}

void Coalescent_calculator::compute_rate_changes() {
    rate_changes.clear();
    float lb = 0;
    float ub = 0;
    for (Branch b : branches) {
        lb = max(cut_time, b.lower_node->time);
        ub = b.upper_node->time;
        if (rate_changes.count(lb) > 0) {
            rate_changes[lb] += 1;
        } else {
            rate_changes[lb] = 1;
        }
        if (rate_changes.count(ub) > 0) {
            rate_changes[ub] -= 1;
        } else {
            rate_changes[ub] = -1;
        }
    }
}

void Coalescent_calculator::compute_rates() {
    rates.clear();
    int curr_rate = 0;
    for (auto x : rate_changes) {
        curr_rate += x.second;
        rates[x.first] = curr_rate;
    }
}

void Coalescent_calculator::compute_probs_quantiles() {
    probs.clear();
    quantiles.clear();
    int curr_rate = 0;
    float prev_time = 0.0;
    float next_time = 0.0;
    float prev_prob = 1.0;
    float next_prob = 1.0;
    float cum_prob = 0.0;
    for (auto it = rates.begin(); it != prev(rates.end()); it++) {
        curr_rate = it->second;
        prev_time = it->first;
        next_time = next(it)->first;
        next_prob = prev_prob*exp(-curr_rate*(next_time - prev_time));
        cum_prob += (prev_prob - next_prob)/curr_rate;
        probs[next_time] = cum_prob;
        quantiles[cum_prob] = next_time;
        prev_prob = next_prob;
    }
    probs[cut_time] = 0;
    quantiles[0] = cut_time;
}


float Coalescent_calculator::prob(float x) {
    if (probs.count(x) > 0) {
        return probs[x];
    }
    auto u_it = probs.upper_bound(x);
    auto l_it = probs.upper_bound(x);
    l_it--;
    int rate = rates[l_it->first];
    float delta_t = u_it->first - l_it->first;
    float delta_p = u_it->second - l_it->second;
    float new_delta_t = x - l_it->first;
    float new_delta_p = delta_p*(1 - exp(-rate*new_delta_t))/(1 - exp(-rate*delta_t));
    return l_it->second + new_delta_p;
}

float Coalescent_calculator::quantile(float p) {
    if (quantiles.count(p) > 0) {
        return quantiles[p];
    }
    auto u_it = quantiles.upper_bound(p);
    auto l_it = quantiles.upper_bound(p);
    l_it--;
    int rate = rates[l_it->second];
    float delta_t = u_it->second - l_it->second;
    float delta_p = u_it->first - l_it->first;
    float new_delta_p = p - l_it->first;
    float new_delta_t = 1 - new_delta_p/delta_p*(1 - exp(-rate*delta_t));
    new_delta_t = -log(new_delta_t)/rate;
    assert(!isnan(new_delta_t));
    return l_it->second + new_delta_t;
}
