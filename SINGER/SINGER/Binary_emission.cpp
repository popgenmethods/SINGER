//
//  Binary_emission.cpp
//  SINGER
//
//  Created by Yun Deng on 4/6/23.
//

#include "Binary_emission.hpp"

Binary_emission::Binary_emission() {}

Binary_emission::~Binary_emission() {}

float Binary_emission::null_emit(Branch branch, float time, float theta, Node *node) {
    float emit_prob = 1;
    float lu = branch.upper_node->time - time;
    float ll = time - branch.lower_node->time;
    float l0 = time - node->time;
    emit_prob = poisson_prob(theta, 1, lu, ll, l0, 0, 0, 0);
    return emit_prob;
}

float Binary_emission::mut_emit(Branch branch, float time, float theta, float bin_size, set<float> mut_set, Node *node) {
    return 0;
}

float Binary_emission::poisson_prob(float theta, float bin_size, float ll, float lu, float l0, float sl, float su, float s0) {
    float prob = 1;
    vector<int> diff = get_diff(sl, su, s0);
    prob *= poisson_prob(theta, bin_size, diff[0]);
    prob *= poisson_prob(theta, bin_size, diff[1]);
    prob *= poisson_prob(theta, bin_size, diff[2]);
    return prob;
}

float Binary_emission::poisson_prob(float theta, float bin_size, int s) {
    return exp(-theta)*pow(theta/bin_size, s);
}

vector<int> Binary_emission::get_diff(float sl, float su, float s0) {
    float sm = 0;
    if (sl + su + s0 > 1.5) {
        sm = 1;
    } else {
        sm = 0;
    }
    vector<int> diff(3);
    diff[0] = abs(sm - sl);
    diff[1] = abs(sm - su);
    diff[2] = abs(sm - s0);
    return diff;
}
