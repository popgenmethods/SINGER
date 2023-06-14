//
//  RSP_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#include "RSP_smc.hpp"

RSP_smc::RSP_smc() {
}

float RSP_smc::sample_start_time(Branch b, int density, float join_time, float cut_time) {
    float lb = b.lower_node->time;
    float ub = b.upper_node->time;
    lb = max(cut_time, lb);
    ub = min(join_time, ub);
    assert(lb < ub);
    float p;
    float t;
    float w;
    float start_time = 0;
    vector<float> start_times = {};
    vector<float> weights = {};
    float weight_sum = 0;
    for (int i = 0; i < density; i++) {
        // t = random()*(ub - lb) + lb;
        t = random_time(lb, ub);
        w = recomb_pdf(t, ub);
        // w = recomb_pdf(t, join_time)*exp(- ub + t);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
    }
    p = random();
    weight_sum = weight_sum*p;
    for (int i = 0; i < density; i++) {
        weight_sum -= weights[i];
        if (weight_sum <= 0) {
            start_time = start_times[i];
            break;
        }
    }
    assert(start_time > cut_time);
    return start_time;
}

pair<Branch, float> RSP_smc::sample_start_time(Branch b1, Branch b2, int density, float join_time, float cut_time) {
    float lb1 = max(b1.lower_node->time, cut_time);
    float ub1 = min(b1.upper_node->time, join_time);
    float lb2 = max(b2.lower_node->time, cut_time);
    float ub2 = min(b2.upper_node->time, join_time);
    assert(lb1 < ub1);
    assert(lb2 < ub2);
    float q = (ub1 - lb1)/(ub1 + ub2 - lb1 - lb2);
    int n1 = round(density*q);
    int n2 = density - n1;
    float p;
    float t;
    float w;
    Branch source_branch;
    float start_time = 0;
    vector<float> start_times(0);
    vector<float> weights(0);
    vector<int> branch_indices(0);
    float weight_sum = 0;
    for (int i = 0; i < n1; i++) {
        // t = random()*(ub1 - lb1) + lb1;
        t = random_time(lb1, ub1);
        w = recomb_pdf(t, ub1);
        // w = recomb_pdf(t, join_time)*exp(- ub1 + t);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
        branch_indices.push_back(1);
    }
    for (int i = 0; i < n2; i++) {
        // t = random()*(ub2 - lb2) + lb2;
        t = random_time(lb2, ub2);
        w = recomb_pdf(t, ub2);
        // w = recomb_pdf(t, join_time)*exp(- ub2 + t);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
        branch_indices.push_back(2);
    }
    p = random();
    weight_sum = weight_sum*p;
    for (int i = 0; i < density; i++) {
        weight_sum -= weights[i];
        if (weight_sum <= 0) {
            start_time = start_times[i];
            if (branch_indices[i] == 1) {
                source_branch = b1;
            } else {
                source_branch = b2;
            }
            break;
        }
    }
    assert(weight_sum <= 0);
    assert(start_time > cut_time);
    return {source_branch, start_time};
}

void RSP_smc::sample_recombination(Recombination &r, float cut_time, Tree &tree) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    get_coalescence_rate(tree, r, cut_time);
    vector<Branch> source_candidates;
    for (Branch b : r.deleted_branches) {
        if (b.upper_node == r.deleted_node and b.lower_node->time < r.inserted_node->time) {
            Branch candidate_recombined_branch = Branch(b.lower_node, r.inserted_node);
            if (r.create(candidate_recombined_branch)) {
                source_candidates.push_back(b);
            }
        }
    }
    if (source_candidates.size() == 1) {
        r.source_branch = source_candidates[0];
        r.start_time = sample_start_time(r.source_branch, 20, r.inserted_node->time, cut_time);
    } else if (source_candidates.size() == 2) {
        pair<Branch, float> breakpoint = sample_start_time(source_candidates[0], source_candidates[1], 40, r.inserted_node->time, cut_time);
        r.source_branch = breakpoint.first;
        r.start_time = breakpoint.second;
    } else {
        cout << r.pos << " " << source_candidates.size() << endl;
        cerr << "no candidates in smc sampling" << endl;
        exit(1);
    }
    r.find_target_branch();
    r.find_recomb_info();
    assert(r.target_branch != Branch());
    assert(r.merging_branch != Branch());
    assert(r.start_time <= r.inserted_node->time);
    assert(r.start_time > cut_time);
}

// private methods:

float RSP_smc::recomb_pdf(float s, float t) {
    float pdf = 1.0;
    float curr_time = s;
    float next_coalescence_time = s;
    map<float, int>::iterator rate_it = coalescence_rates.upper_bound(s);
    rate_it--;
    int rate;
    while (next_coalescence_time < t) {
        rate = rate_it->second;
        rate_it++;
        next_coalescence_time = rate_it->first;
        pdf *= exp(-rate*(-curr_time + min(t, next_coalescence_time)));
        curr_time = next_coalescence_time;
    }
    return pdf;
}

void RSP_smc::get_coalescence_rate(Tree &tree, Recombination &r, float cut_time) {
    coalescence_rates.clear();
    vector<float> coalescence_times = {cut_time};
    for (Branch b : tree.branches) {
        if (b.lower_node->time > cut_time and b.lower_node != r.deleted_node) {
            coalescence_times.push_back(b.lower_node->time);
        }
    }
    coalescence_times.push_back(numeric_limits<float>::infinity());
    sort(coalescence_times.begin(), coalescence_times.end());
    int n = (int) coalescence_times.size();
    for (int i = 0; i < n; i++) {
        coalescence_rates[coalescence_times[i]] = n-i-1;
    }
}

float RSP_smc::random() {
    float p = uniform_random();
    p = 0.1 + 0.9*p;
    return p;
}

float RSP_smc::random_time(float lb, float ub) {
    /*
    float t = 0;
    float p = random();
    float q = p + (1 - p)*exp(lb - ub);
    float delta_t = -log(q);
    t = ub - delta_t;
    assert(t > lb);
    return t;
     */
    float t = random()*(ub - lb) + lb;
    /*
    if (t <= lb) {
        t = ub - 1e-12;
    }
    assert(t > lb);
     */
    return t;
}
