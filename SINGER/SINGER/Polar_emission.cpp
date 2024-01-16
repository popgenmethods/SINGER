//
//  Polar_emission.cpp
//  SINGER
//
//  Created by Yun Deng on 6/14/23.
//

#include "Polar_emission.hpp"

Polar_emission::Polar_emission() {}

Polar_emission::~Polar_emission() {}

/*
double Polar_emission::null_emit(Branch &branch, double time, double theta, Node_ptr node) {
    double emit_prob = 1;
    double old_prob = 1;
    double ll = time - branch.lower_node->time;
    double lu = branch.upper_node->time - time;
    double l0 = time - node->time;
    emit_prob = null_prob(theta, ll, lu, l0);
    if (!isinf(lu)) {
        old_prob = null_prob(theta*(ll + lu));
    } else {
        old_prob = 1;
    }
    emit_prob /= old_prob;
    return emit_prob;
}
 */

double Polar_emission::null_emit(Branch &branch, double time, double theta, Node_ptr node) {
    double emit_prob = 1;
    double ll = time - branch.lower_node->time;
    double lu = branch.upper_node->time - time;
    double l0 = time - node->time;
    if (!isinf(lu)) {
        emit_prob = null_prob(theta*l0);
    } else {
       emit_prob = null_prob(theta*(ll + l0));
    }
    return emit_prob;
}

/*
double Polar_emission::mut_emit(Branch &branch, double time, double theta, double bin_size, set<double> &mut_set, Node_ptr node) {
    double emit_prob = 1;
    double old_prob = 1;
    double ll = time - branch.lower_node->time;
    double lu = branch.upper_node->time - time;
    double l0 = time - node->time;
    for (double m : mut_set) {
        get_diff(m, branch, node);
        emit_prob *= mut_prob(theta, bin_size, ll, lu, l0, diff[0], diff[1], diff[2]);
        old_prob *= mut_prob(theta*(ll + lu), bin_size, diff[3]);
    }
    emit_prob *= null_prob(theta, ll, lu, l0);
    old_prob *= null_prob(theta*(ll + lu));
    emit_prob /= old_prob;
    emit_prob *= root_reward;
    emit_prob = max(emit_prob, 1e-20f);
    assert(emit_prob > 0);
    return emit_prob;
}
 */

double Polar_emission::mut_emit(Branch &branch, double time, double theta, double bin_size, set<double> &mut_set, Node_ptr node) {
    double emit_prob = 1;
    double old_prob = 1;
    double ll = time - branch.lower_node->time;
    double lu = branch.upper_node->time - time;
    double l0 = time - node->time;
    for (double m : mut_set) {
        get_diff(m, branch, node);
        emit_prob *= mut_prob(theta, bin_size, ll, lu, l0, diff[0], diff[1], diff[2]);
        old_prob *= mut_prob(theta*(ll + lu), bin_size, diff[3]);
    }
    if (!isinf(lu)) {
       emit_prob *= null_prob(theta*l0);
    } else {
       emit_prob *= null_prob(theta*(l0 + ll));
    }
    emit_prob /= old_prob;
    emit_prob *= root_reward;
    emit_prob = max(emit_prob, 1e-20);
    assert(emit_prob > 0);
    return emit_prob;
}

double Polar_emission::emit(Branch &branch, double time, double theta, double bin_size, vector<double> &emissions, Node_ptr node) {
    double emit_prob = 1;
    double old_prob = 1;
    double ll = time - branch.lower_node->time;
    double lu = branch.upper_node->time - time;
    double l0 = time - node->time;
    emit_prob = mut_prob(theta, bin_size, ll, lu, l0, emissions[0], emissions[1], emissions[2]);
    old_prob = mut_prob(theta*(ll + lu), bin_size, emissions[3]);
    emit_prob *= null_prob(theta, ll, lu, l0);
    old_prob *= null_prob(theta*(ll + lu));
    emit_prob /= old_prob;
    return emit_prob;
}

double Polar_emission::mut_prob(double theta, double bin_size, double ll, double lu, double l0, int sl, int su, int s0) {
    double prob = 1;
    prob *= mut_prob(ll*theta, bin_size, sl);
    prob *= mut_prob(lu*theta, bin_size, su);
    prob *= mut_prob(l0*theta, bin_size, s0);
    if (s0 >= 1) {
        prob *= penalty;
    }
    return prob;
}

double Polar_emission::null_prob(double theta, double ll, double lu, double l0) {
    double prob = 1;
    prob *= null_prob(ll*theta);
    if (!isinf(lu)) {
        prob *= null_prob(lu*theta);
    }
    prob *= null_prob(l0*theta);
    return prob;
}

double Polar_emission::mut_prob(double theta, double bin_size, int s) {
    if (isinf(theta)) {
        return 1.0;
    }
    double unit_theta = theta/bin_size;
    return pow(unit_theta, abs(s));
}

double Polar_emission::null_prob(double theta) {
    if (isinf(theta)) {
        return 1;
    }
    return exp(-theta);
}

void Polar_emission::get_diff(double m, Branch branch, Node_ptr node) {
    double sl = 0;
    double su = 0;
    double s0 = 0;
    double sm = 0;
    sl = branch.lower_node->get_state(m);
    su = branch.upper_node->get_state(m);
    s0 = node->get_state(m);
    if (sl + su + s0 > 1.5) {
        sm = 1;
    } else {
        sm = 0;
    }
    if (branch.upper_node->index == -1) {
        if (sm == 0 and sl == 1) {
            root_reward = ancestral_prob/(1 - ancestral_prob);
        } else {
            root_reward = 1;
        }
    } else {
        root_reward = 1;
    }
    // remember: lower - upper, to keep the directionality of mutations
    diff[0] = sl - sm;
    diff[1] = sm - su;
    diff[2] = s0 - sm;
    diff[3] = sl - su;
}
