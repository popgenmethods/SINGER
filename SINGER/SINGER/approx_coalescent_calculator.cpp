//
//  approx_coalescent_calculator.cpp
//  SINGER
//
//  Created by Yun Deng on 6/25/23.
//

#include "approx_coalescent_calculator.hpp"

approx_coalescent_calculator::approx_coalescent_calculator(float t) {
    cut_time = t;
}

void approx_coalescent_calculator::start(Tree &tree) {
    for (auto &x : tree.parents) {
        if (x.second->time > cut_time) {
            coalescence_times.insert(x.second->time);
        }
    }
}

void approx_coalescent_calculator::update(Recombination &r) {
    if (r.deleted_node->time > cut_time) {
        coalescence_times.erase(r.deleted_node->time);
    } else {
        coalescence_times.erase(r.inserted_node->time);
    }
}

float approx_coalescent_calculator::prob(float x, float y) {
    int n0 = (int) coalescence_times.size();
    if (n0 == 1) {
        return exp(-x) - exp(-y);
    }
    float p = 1.0/(n0 + (1 - n0)*exp(-0.5*y)) - 1.0/(n0 + (1 - n0)*exp(-0.5*x));
    p *= -2.0/(n0 - 1);
    return p;
}
