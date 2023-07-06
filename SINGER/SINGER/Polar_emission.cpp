//
//  Polar_emission.cpp
//  SINGER
//
//  Created by Yun Deng on 6/14/23.
//

#include "Polar_emission.hpp"

Polar_emission::Polar_emission() {}

Polar_emission::~Polar_emission() {}

float Polar_emission::null_emit(Branch &branch, float time, float theta, Node_ptr node) {
    float emit_prob = 1;
    float old_prob = 1;
    float ll = time - branch.lower_node->time;
    float lu = branch.upper_node->time - time;
    float l0 = time - node->time;
    emit_prob = null_prob(theta, ll, lu, l0);
    if (!isinf(lu)) {
        old_prob = null_prob(theta*(ll + lu));
    } else {
        old_prob = 1;
    }
    // old_prob = null_prob(theta*(ll + lu));
    emit_prob /= old_prob;
    return emit_prob;
}

float Polar_emission::mut_emit(Branch &branch, float time, float theta, float bin_size, set<float> &mut_set, Node_ptr node) {
    float emit_prob = 1;
    float old_prob = 1;
    float ll = time - branch.lower_node->time;
    float lu = branch.upper_node->time - time;
    float l0 = time - node->time;
    for (float m : mut_set) {
        get_diff(m, branch, node);
        emit_prob *= mut_prob(theta, bin_size, ll, lu, l0, diff[0], diff[1], diff[2]);
        old_prob *= mut_prob(theta*(ll + lu), bin_size, diff[3]);
    }
    emit_prob *= null_prob(theta, ll, lu, l0);
    old_prob *= null_prob(theta*(ll + lu));
    emit_prob /= old_prob;
    assert(emit_prob != 0);
    return emit_prob;
}

float Polar_emission::emit(Branch &branch, float time, float theta, float bin_size, vector<float> &emissions, Node_ptr node) {
    float emit_prob = 1;
    float old_prob = 1;
    float ll = time - branch.lower_node->time;
    float lu = branch.upper_node->time - time;
    float l0 = time - node->time;
    emit_prob = mut_prob(theta, bin_size, ll, lu, l0, emissions[0], emissions[1], emissions[2]);
    old_prob = mut_prob(theta*(ll + lu), bin_size, emissions[3]);
    emit_prob *= null_prob(theta, ll, lu, l0);
    old_prob *= null_prob(theta*(ll + lu));
    emit_prob /= old_prob;
    // assert(emit_prob != 0);
    return emit_prob;
}

float Polar_emission::mut_prob(float theta, float bin_size, float ll, float lu, float l0, int sl, int su, int s0) {
    float prob = 1;
    prob *= mut_prob(ll*theta, bin_size, sl);
    prob *= mut_prob(lu*theta, bin_size, su);
    prob *= mut_prob(l0*theta, bin_size, s0);
    return prob;
}

float Polar_emission::null_prob(float theta, float ll, float lu, float l0) {
    float prob = 1;
    prob *= null_prob(ll*theta);
    if (!isinf(lu)) {
        prob *= null_prob(lu*theta);
    }
    prob *= null_prob(l0*theta);
    return prob;
}

float Polar_emission::mut_prob(float theta, float bin_size, int s) {
    if (isinf(theta)) {
        return 1.0;
    }
    float unit_theta = theta/bin_size;
    if (s < 0) {
        unit_theta *= reverse_penalty;
    }
    return pow(unit_theta, abs(s));
}

float Polar_emission::null_prob(float theta) {
    if (isinf(theta)) {
        return 1;
    }
    return exp(-theta);
}

void Polar_emission::get_diff(float m, Branch branch, Node_ptr node) {
    float sl = 0;
    float su = 0;
    float s0 = 0;
    float sm = 0;
    sl = branch.lower_node->get_state(m);
    su = branch.upper_node->get_state(m);
    s0 = node->get_state(m);
    if (sl + su + s0 > 1.5) {
        sm = 1;
    } else {
        sm = 0;
    }
    // remember: lower - upper, to keep the directionality of mutations
    diff[0] = sl - sm;
    diff[1] = sm - su;
    diff[2] = s0 - sm;
    diff[3] = sl - su;
}
