//
//  approx_coalescent_calculator.cpp
//  SINGER
//
//  Created by Yun Deng on 6/25/23.
//

#include "approx_coalescent_calculator.hpp"

approx_coalescent_calculator::approx_coalescent_calculator(double t) {
    cut_time = t;
}

approx_coalescent_calculator::~approx_coalescent_calculator() {}

void approx_coalescent_calculator::start(set<Branch> &branches) {
    for (auto &x : branches) {
        if (x.upper_node->time > cut_time and x.lower_node->time <= cut_time) {
            n0 += 1;
        }
    }
}

void approx_coalescent_calculator::start(Tree &tree) {
    for (auto &x : tree.parents) {
        if (x.second->time > cut_time and x.first->time <= cut_time) {
            n0 += 1;
        }
    }
}

void approx_coalescent_calculator::update(Recombination &r) {
    if (r.deleted_node->time > cut_time) {
        n0 -= 1;
    }
    if (r.inserted_node->time > cut_time) {
        n0 += 1;
    }
}

pair<double, double> approx_coalescent_calculator::compute_time_weights(double x, double y) {
    if (x == y) { // no need to compute when point mass
        return {x, 0};
    }
    double t, w, p, q;
    p = prob(x, y);
    t = find_median(x, y);
    q = p*(t - cut_time);
    w = q*n0/2;
    if (y - x < 0.001) {
        t = 0.5*(x + y);
    }
    assert(t >= x and t <= y);
    assert(q >= 0);
    return {t, w};
}

void approx_coalescent_calculator::compute_first_moment() {
    first_moment = 2.0/n0;
}

double approx_coalescent_calculator::prob(double x, double y) {
    double p1 = prob_integral(x);
    double p2 = prob_integral(y);
    return max(p2 - p1, 0.0);
}

double approx_coalescent_calculator::prob_integral(double x) {
    x -= cut_time;
    if (n0 == 1) {
        return 1 - exp(-x);
    }
    double v = n0 + (1 - n0)*exp(-0.5*x);
    double p = -2*log(v) -2*n0/v;
    p = p/(1 - n0)/(1 - n0);
    assert(!isnan(p));
    return p;
}

double approx_coalescent_calculator::find_median(double x, double y) {
    if (x == y) {
        return x;
    }
    if (y - x <= 0.01) {
        return 0.5*(x + y);
    }
    if (n0 == 1) {
        double q = 0.5*(exp(-x) + exp(-y));
        double t = -log(q);
        assert(t >= x and t <= y);
        return t;
    }
    x -= cut_time;
    y -= cut_time;
    double z = 0;
    double nx = n0 + (1 - n0)*exp(-0.5*x);
    double ny = n0 + (1 - n0)*exp(-0.5*y);
    double nm = sqrt(nx*ny);
    z = 2*log(n0 - 1) - 2*log(n0 - nm);
    if (isinf(y)) {
        z = x + 1;
    }
    if (z < x or z > y) {
        z = 0.5*(x + y);
    }
    z += cut_time;
    return z;
}
