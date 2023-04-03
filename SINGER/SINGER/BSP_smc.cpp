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

void BSP_smc::null_emit(float theta, Node *base_node) {
    float ws = 0.0f;
    float emit_prob = 1;
    for (Interval *i : curr_intervals) {
        emit_prob = eh->null_emit(i->branch, i->time, theta, base_node);
        i->multiply(emit_prob);
        ws += i->get_prob();
    }
    for (Interval *i : curr_intervals) {
        i->rescale(ws);
    }
}

void BSP_smc::mut_emit(float theta, float mut_pos, Node *base_node) {
    float emit_prob = 1;
    float ws = 0.0f;
    for (Interval *i : curr_intervals) {
        emit_prob = eh->mut_emit(i->branch, i->time, theta, mut_pos, base_node);
        i->multiply(emit_prob);
        ws += i->get_prob();
    }
    for (Interval *i : curr_intervals) {
        i->rescale(ws);
    }
}

map<float, Branch> BSP_smc::sample_joining_branches() {
    map<float, Branch> joining_branches = {};
    int x = end_pos;
    Interval *interval = sample_terminal_interval(x);
    Branch b;
    while (x >= start_pos) {
        x = trace_back_helper(interval, x);
        b = interval->branch;
        joining_branches.insert({x, b});
        if (x == start_pos) {
            break;
        } else if (x == interval->start_pos) {
            x -= 1;
            interval = interval->sample_source();
            b = interval->branch;
        } else {
            x -= 1;
            interval = sample_prev_interval(x);
            b = interval->branch;
        }
    }
    joining_branches.insert({INT_MAX, joining_branches.rbegin()->second});
    return joining_branches;
}

// private methods:

void BSP_smc::transfer_helper(Interval_info interval_info, Interval *i, float w) {
    assert(!isnan(w));
    if (interval_info.branch == Branch()) {
        return;
    }
    if (w == 0) {
        return;
    }
    if (transfer_weights.count(interval_info) > 0) {
        transfer_weights.at(interval_info).push_back(w);
        transfer_intervals.at(interval_info).push_back(i);
    } else {
        transfer_weights.insert({interval_info, {w}});
        transfer_intervals.insert({interval_info, {i}});
    }
}

void BSP_smc::add_new_branches(Recombination r) { // add recombined branch and merging branch, if legal
    Interval_info ii;
    if (r.merging_branch != Branch() and r.merging_branch.upper_node->time > cut_time) {
        ii = Interval_info(r.merging_branch, max(cut_time, r.merging_branch.lower_node->time), r.merging_branch.upper_node->time);
        transfer_weights.insert({ii, {}});
        transfer_intervals.insert({ii, {}});
    }
    if (r.recombined_branch != Branch() and r.recombined_branch.upper_node->time > cut_time) {
        ii = Interval_info(r.recombined_branch, max(cut_time, r.recombined_branch.lower_node->time), r.recombined_branch.upper_node->time);
        transfer_weights.insert({ii, {}});
        transfer_intervals.insert({ii, {}});
    }
}

void BSP_smc::fill_interval_info() {
    float t;
    float p;
    for (Interval *i : curr_intervals) {
        p = get_prop(i->lb, i->ub);
        t = get_median(i->lb, i->ub);
        i->assign_time(t);
        i->assign_weight(p);
    }
    weight_sums.push_back(0.0);
    sort(curr_intervals.begin(), curr_intervals.end(), compare_interval());
}

void BSP_smc::sanity_check(Recombination r) {
    for (Interval *i : curr_intervals) {
        if (i->lb == i->ub and i->lb == r.inserted_node->time and i->branch != r.target_branch) {
            i->update_prob(0.0);
        }
    }
}

void BSP_smc::generate_intervals(Recombination r) {
    Branch b;
    float lb;
    float ub;
    float p;
    vector<Interval *> source_intervals;
    vector<float> source_weights;
    Interval *new_interval = nullptr;
    for (auto x : transfer_weights) {
        source_weights = x.second;
        source_intervals = transfer_intervals.at(x.first);
        b = x.first.branch;
        lb = x.first.lb;
        ub = x.first.ub;
        p = accumulate(source_weights.begin(), source_weights.end(), 0.0);
        assert(!isnan(p));
        if (lb == max(cut_time, b.lower_node->time) and ub == b.upper_node->time) { // full intervals
            new_interval = new Interval(b, lb, ub, curr_pos, p);
            new_interval->set_source(source_intervals, source_weights);
            curr_intervals.push_back(new_interval);
        } else if (p >= cutoff) { // partial intervals
            new_interval = new Interval(b, lb, ub, curr_pos, p);
            new_interval->set_source(source_intervals, source_weights);
            curr_intervals.push_back(new_interval);
        }
    }
    fill_interval_info();
    transfer_weights.clear();
    transfer_intervals.clear();
}

float BSP_smc::get_prop(float lb, float ub) {
    assert(lb <= ub);
    float x = get_prob(ub) - get_prob(lb);
    assert(!isnan(x));
    return x;
}

float BSP_smc::get_overwrite_prob(Recombination r, float lb, float ub) {
    if (check_points.count(curr_pos) > 0) {
        return 0.0;
    }
    float join_time = r.inserted_node->time;
    float p1 = get_prop(lb, ub);
    float p2 = get_prop(max(cut_time, r.start_time), join_time);
    if (p1 == 0 and p2 == 0) {
        return 1.0;
    }
    float overwrite_prob = p2/(p1 + p2);
    assert(!isnan(overwrite_prob));
    return overwrite_prob;
}

float BSP_smc::get_prob(float x) {
    if (coalescence_probs.count(x) > 0) {
        return coalescence_probs.at(x);
    }
    map<float, float>::iterator u_it = coalescence_probs.upper_bound(x);
    map<float, float>::iterator l_it = coalescence_probs.upper_bound(x);
    l_it--;
    int rate = coalescence_rates.at(l_it->first);
    float delta_t = u_it->first - l_it->first;
    float delta_p = u_it->second - l_it->second;
    float new_delta_t = x - l_it->first;
    float new_delta_p = delta_p*(1 - exp(-rate*new_delta_t))/(1 - exp(-rate*delta_t));
    return l_it->second + new_delta_p;
}

float BSP_smc::get_quantile(float p) {
    if (coalescence_quantiles.count(p) > 0) {
        return coalescence_quantiles.at(p);
    }
    map<float, float>::iterator u_it = coalescence_quantiles.upper_bound(p);
    map<float, float>::iterator l_it = coalescence_quantiles.upper_bound(p);
    l_it--;
    int rate = coalescence_rates.at(l_it->second);
    float delta_t = u_it->second - l_it->second;
    float delta_p = u_it->first - l_it->first;
    float new_delta_p = p - l_it->first;
    float new_delta_t = 1 - new_delta_p/delta_p*(1 - exp(-rate*delta_t));
    new_delta_t = -log(new_delta_t)/rate;
    assert(l_it->second + new_delta_t != numeric_limits<float>::infinity());
    return l_it->second + new_delta_t;
}

float BSP_smc::get_median(float lb, float ub) {
    float lq = get_prob(lb);
    float uq = get_prob(ub);
    float t;
    if (ub == numeric_limits<float>::infinity()) {
        return lb + log(2);
    }
    if (ub - lb < 1e-3) {
        t = 0.5*(lb + ub);
    } else if (uq - lq < 1e-3) {
        t = 0.5*(lb + ub);
    } else {
        t = get_quantile(0.5*(lq + uq));
    }
    assert(t >= lb and t <= ub);
    return t;
}
