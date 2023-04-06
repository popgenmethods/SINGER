//
//  TSP_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 4/5/23.
//

#include "TSP_smc.hpp"

TSP_smc::TSP_smc() {
}

TSP_smc::~TSP_smc() {
    for (auto x : state_spaces) {
        for (Interval *i : x.second) {
            delete(i);
        }
    }
}

void TSP_smc::set_coordinates(vector<float> c) {
    coordinates = c;
}

void TSP_smc::set_gap(float q) {
    gap = q;
}

void TSP_smc::set_emission(shared_ptr<Emission> e) {
    eh = e;
}

void TSP_smc::set_check_points(set<float> p) {
    check_points = p;
}

void TSP_smc::start(Branch branch, float t) {
    cut_time = t;
    curr_index = 0;
    if (branch == Branch()) {
        query_time = cut_time;
    } else {
        query_time = max(cut_time, branch.lower_node->time);
    }
    generate_intervals(branch, branch.lower_node->time, branch.upper_node->time);
    for (Interval *i : curr_intervals) {
        i->update_prob(exp(-i->lb) - exp(-i->ub));
    }
    for (Interval *i : curr_intervals) {
        assert(i->lb >= cut_time);
    }
    state_spaces[0] = curr_intervals;
    set_quantities();
}

void TSP_smc::transfer(Recombination &r, Branch prev_branch, Branch next_branch) {
    rhos.push_back(0);
    prev_rho = -1;
    sanity_check(r);
    if (next_branch == Branch()) {
        query_time = cut_time;
    } else {
        query_time = max(cut_time, next_branch.lower_node->time);
    }
    if (next_branch == Branch()) {
        ;
    } else if (prev_branch == r.source_branch and next_branch == r.merging_branch) {
        set_interval_constraint(r);
    } else if (prev_branch == r.target_branch and next_branch == r.recombined_branch) {
        set_point_constraint(r);
    }
    curr_intervals.clear();
    if (next_branch == Branch()) {
        ;
    } else if (prev_branch == r.source_branch and next_branch == r.merging_branch) { // switch to a point mass
        float t = r.deleted_node->time;
        generate_intervals(next_branch, next_branch.lower_node->time, t);
        generate_intervals(next_branch, t, t);
        curr_intervals.back()->update_prob(1.0);
        curr_intervals.back()->set_node(r.deleted_node);
        generate_intervals(next_branch, t, next_branch.upper_node->time);
    } else if (prev_branch == r.target_branch and next_branch == r.recombined_branch) { // switch from a point mass
        generate_intervals(next_branch, next_branch.lower_node->time, r.start_time);
        generate_intervals(next_branch, r.start_time, next_branch.upper_node->time);
        for (Interval *i : curr_intervals) {
            if (i->time >= r.start_time) {
                i->update_prob(1);
            }
        }
    } else {
        float lb;
        float ub;
        lb = next_branch.lower_node->time;
        ub = max(prev_branch.lower_node->time, next_branch.lower_node->time);
        generate_intervals(next_branch, lb, ub);
        transfer_intervals(r, prev_branch, next_branch);
        if (curr_intervals.size() > 0) {
            lb = min(curr_intervals.back()->ub, next_branch.upper_node->time);
        } else {
            lb = cut_time;
        }
        ub = next_branch.upper_node->time;
        generate_intervals(next_branch, lb, ub);
    }
    state_spaces[curr_index] = curr_intervals;
    for (Interval *i : curr_intervals) {
        assert(i->lb >= cut_time);
        assert(i->time >= i->lb and i->time <= i->ub);
    }
    set_quantities();
}

void TSP_smc::recombine(Branch prev_branch, Branch next_branch) {
    assert(next_branch != Branch());
    vector<Interval *> prev_intervals = curr_intervals;
    curr_intervals.clear();
    rhos.push_back(0);
    prev_rho = -1;
    curr_index += 1;
    query_time = max(cut_time, next_branch.lower_node->time);
    generate_intervals(next_branch, next_branch.lower_node->time, next_branch.upper_node->time);
    state_spaces[curr_index] = curr_intervals;
    set_quantities();
    float new_prob;
    float base;
    for (int i = 0; i < prev_intervals.size(); i++) {
        base = recomb_prob(prev_intervals[i]->time, curr_intervals.front()->lb, curr_intervals.back()->ub);
        for (int j = 0; j < curr_intervals.size(); j++) {
            if (base == 0) {
                new_prob = 1;
            } else {
                new_prob = recomb_prob(prev_intervals[i]->time, curr_intervals[j]->lb, curr_intervals[j]->ub)*prev_intervals[i]->get_prob()/base;
            }
            curr_intervals[j]->add_prob(new_prob);
        }
    }
    for (Interval *i : curr_intervals) {
        assert(i->lb >= cut_time);
        assert(i->time >= i->lb and i->time <= i->ub);
    }
}

float TSP_smc::get_exp_quantile(float p) {
    assert(p >= 0 and p <= 1);
    if (p < 1e-6) {
        return 0;
    }
    if (1 - p < 1e-6) {
        return numeric_limits<float>::infinity();
    }
    return -log(1 - p);
}

vector<float> TSP_smc::generate_grid(float lb, float ub) {
    assert(lb < ub);
    vector<float> points = {lb};
    float lq = 1 - exp(-lb);
    float uq = 1 - exp(-ub);
    float q = uq - lq;
    int n = ceil(q/gap);
    float l;
    for (int i = 1; i < n; i++) {
        l = get_exp_quantile(lq + i*q/n);
        points.push_back(l);
    }
    points.push_back(ub);
    return points;
}

float TSP_smc::recomb_cdf(float s, float t) {
    if (t == numeric_limits<float>::infinity()) {
        return 1;
    }
    if (t == 0) {
        return 0;
    }
    float cdf = 0;
    float l = s - cut_time;
    if (s > t) {
        cdf = t - 1 + exp(cut_time - t) - cut_time;
    } else {
        cdf = s + exp(cut_time - t) - exp(s - t) - cut_time;
    }
    cdf = cdf/l;
    assert(!isnan(cdf));
    return cdf;
}

float TSP_smc::recomb_quantile(float s, float q, float lb, float ub) {
    float mid_q = 0;
    float mid = lb;
    float p = 0;
    float w = 0;
    float gap = min(0.5*(ub - lb), 1e-5);
    while (ub - mid > gap) {
        p = recomb_cdf(s, mid);
        if (p >= q) {
            lb = lb;
            ub = mid;
        } else {
            lb = mid;
            ub = ub;
        }
        w = random()*0.02 + 0.49;
        mid_q = w*(exp(-lb) + exp(-ub));
        mid = -log(mid_q);
    }
    return mid;
}

void TSP_smc::forward(float rho) {
    rhos.push_back(rho);
    curr_index += 1;
    if (curr_intervals.size() == 0) {
        return;
    }
    compute_diagonals(rho);
    compute_lower_diagonals(rho);
    compute_upper_diagonals(rho);
    compute_lower_sums();
    compute_upper_sums();
    prev_rho = rho;
    float new_prob;
    for (int i = 0; i < curr_intervals.size(); i++) {
        new_prob = lower_sums[i] + diagonals[i]*curr_intervals[i]->get_prob() + lower_diagonals[i]*upper_sums[i];
        assert(!isnan(new_prob) and new_prob >= 0);
        curr_intervals[i]->new_prob(new_prob);
    }
}

void TSP_smc::null_emit(float theta, Node *query_node) {
    float ws = 0.0;
    float emit_prob = 0.0;
    for (Interval *i : curr_intervals) {
        emit_prob = eh->null_emit(i->branch, i->time, theta, query_node);
        i->multiply(emit_prob);
        ws += i->get_prob();
    }
    if (ws > 0) {
        for (Interval *i : curr_intervals) {
            i->rescale(ws);
        }
    } else {
        float p = 1.0/curr_intervals.size();
        for (Interval *i : curr_intervals) {
            i->update_prob(p);
        }
    }
}

void TSP_smc::mut_emit(float theta, float mut_pos, Node *query_node) {
    float ws = 0.0f;
    float emit_prob = 0.0;
    for (Interval *i : curr_intervals) {
        emit_prob = eh->mut_emit(i->branch, i->time, theta, mut_pos, query_node);
        i->multiply(emit_prob);
        ws += i->get_prob();
    }
    if (ws > 0) {
        for (Interval *i : curr_intervals) {
            i->rescale(ws);
        }
    } else {
        float p = 1.0/curr_intervals.size();
        for (Interval *i : curr_intervals) {
            i->update_prob(p);
        }
    }
}

float TSP_smc::recomb_prob(float s, float t1, float t2) {
    assert(t1 <= t2);
    assert(t1 >= cut_time and s >= cut_time);
    return recomb_cdf(s, t2) - recomb_cdf(s, t1);
}

float TSP_smc::psmc_cdf(float rho, float s, float t) {
    float l;
    float pre_factor;
    if (t <= s) {
        l = 2*t - query_time - cut_time;
    } else {
        l = 2*s - query_time - cut_time;
    }
    if (l == 0) {
        pre_factor = rho;
    } else {
        pre_factor = (1 - exp(-rho*l))/l;
    }
    float integral;
    float cdf;
    if (t == cut_time and t == query_time) {
        return 0;
    } else if (t <= s) {
        integral = 2*t + exp(-t)*(exp(cut_time) + exp(query_time)) - cut_time - query_time - 2;
    } else {
        integral = 2*s + exp(cut_time - t) + exp(query_time - t) - 2*exp(s-t) - cut_time - query_time;
    }
    cdf = pre_factor*integral;
    return cdf;
}

float TSP_smc::standard_recomb_cdf(float rho, float s, float t) {
    float l = 2*s;
    float integral;
    float cdf;
    if (t <= s) {
        integral = t + 0.5*exp(-2*t) - 0.5;
    } else {
        integral = s + 0.5*exp(-2*s) - 0.5 + (1 - exp(s-t))*(1 - exp(-2*s));
    }
    cdf = integral/l;
    return cdf;
}

float TSP_smc::psmc_prob(float rho, float s, float t1, float t2) {
    assert(s != numeric_limits<float>::infinity());
    assert(t1 <= t2);
    assert(t1 >= query_time and s >= query_time);
    float prob;
    float base;
    float l = 2*s - query_time - cut_time;
    if (t1 == s and t2 == s) {
        base = exp(-rho*l);
    } else if (t1 < s and t2 > s) {
        base = exp(-rho*l);
    } else {
        base = 0;
    }
    float gap;
    float uq = 0;
    float lq = 0;
    if (t2 - t1 > 0) {
        uq = psmc_cdf(rho, s, t2);
        lq = psmc_cdf(rho, s, t1);
        gap = uq - lq;
    } else {
        gap = 0;
    }
    gap = max(gap, 0.0f);
    prob = base + gap;
    assert(!isnan(prob) and prob >= 0 and prob <= 1);
    return base + gap;
}

void TSP_smc::generate_intervals(Branch next_branch, float lb, float ub) {
    assert(lb >= next_branch.lower_node->time and ub <= next_branch.upper_node->time);
    lb = max(cut_time, lb);
    ub = max(cut_time, ub);
    if (lb == ub) {
        if (lb == next_branch.lower_node->time or lb == next_branch.upper_node->time or lb == cut_time) {
            return;
        }
        else {
            Interval *new_interval = new Interval(next_branch, lb, ub, curr_index, 0);
            new_interval->fill_time();
            curr_intervals.push_back(new_interval);
            return;
        }
    }
    vector<float> points = generate_grid(lb, ub);
    float l;
    float u;
    for (int i = 0; i < points.size() - 1; i++) {
        l = points[i];
        u = points[i+1];
        Interval *new_interval = new Interval(next_branch, l, u, curr_index, 0);
        new_interval->fill_time();
        curr_intervals.push_back(new_interval);
    }
}

void TSP_smc::transfer_intervals(Recombination &r, Branch prev_branch, Branch next_branch) {
    float lb;
    float ub;
    float p;
    vector<Interval *> prev_intervals = get_state_space(curr_index);
    for (Interval *i : prev_intervals) {
        lb = max(i->lb, next_branch.lower_node->time);
        ub = min(i->ub, next_branch.upper_node->time);
        if (prev_branch == r.source_branch) {
            ub = min(ub, r.start_time);
            if (lb == r.start_time) {
                continue;
            }
        }
        if (lb == ub and ub == next_branch.upper_node->time) {
            continue;
        }
        if (lb == ub and ub == next_branch.lower_node->time) {
            continue;
        }
        if (lb == ub or ub - lb > 1e-6) {
            float w = get_prop(lb, ub, i->lb, i->ub);
            assert(!isnan(w) and w >= 0);
            p = w*i->get_prob();
            Interval *new_interval = new Interval(next_branch, lb, ub, curr_index, p);
            new_interval->set_source({i}, {1});
            new_interval->fill_time();
            if (i->node != nullptr) {
                new_interval->set_node(i->node);
            }
            curr_intervals.push_back(new_interval);
        }
    }
}

void TSP_smc::set_quantities() {
    int n = (int) curr_intervals.size();
    diagonals = vector<float>(n);
    lower_diagonals = vector<float>(n);
    upper_diagonals = vector<float>(n);
    lower_sums = vector<float>(n);
    upper_sums = vector<float>(n);
}

float TSP_smc::random() {
    float p = (float) rand()/RAND_MAX;
    return p;
}

float TSP_smc::get_prop(float lb1, float ub1, float lb2, float ub2) {
    float p;
    if (ub2 - lb2 < 1e-6) {
        p = 1;
    } else {
        float p1 = exp(-lb1) - exp(-ub1);
        float p2 = exp(-lb2) - exp(-ub2);
        assert(!isnan(p1));
        assert(p2 != 0);
        p = p1/p2;
    }
    assert(p>=0 and p<=1);
    return p;
}

void TSP_smc::compute_diagonals(float rho) {
    if (rho == prev_rho) {
        return;
    }
    float t;
    float base;
    float lb = curr_intervals.front()->lb;
    float ub = curr_intervals.back()->ub;
    Interval *curr_interval = nullptr;
    float diag = 0;
    for (int i = 0; i < diagonals.size(); i++) {
        curr_interval = curr_intervals[i];
        t = curr_interval->time;
        base = psmc_prob(rho, curr_interval->time, lb, ub);
        assert(base > 0 and base <= 1);
        diag = psmc_prob(rho, curr_interval->time, curr_interval->lb, curr_interval->ub);
        diagonals[i] = diag/base;
        assert(diagonals[i] > 0);
    }
}

void TSP_smc::compute_lower_diagonals(float rho) {
    if (rho == prev_rho) {
        return;
    }
    float t;
    float base;
    int n = (int) lower_diagonals.size();
    lower_diagonals[n-1] = 0;
    float lb = max(cut_time, curr_intervals.front()->lb);
    float ub = curr_intervals.back()->ub;
    for (int i = 0; i < lower_diagonals.size() - 1; i++) {
        t = curr_intervals[i+1]->time;
        base = psmc_prob(rho, t, lb, ub);
        lower_diagonals[i] = psmc_prob(rho, t, curr_intervals[i]->lb, curr_intervals[i]->ub)/base;
        assert(!isnan(lower_diagonals[i]) and lower_diagonals[i] >= 0);
    }
}

void TSP_smc::compute_upper_diagonals(float rho) {
    if (rho == prev_rho) {
        return;
    }
    upper_diagonals[0] = 0;
    float lb = max(cut_time, curr_intervals.front()->lb);
    float ub = curr_intervals.back()->ub;
    float t;
    float base;
    for (int i = 1; i < upper_diagonals.size(); i++) {
        t = curr_intervals[i-1]->time;
        base = psmc_prob(rho, t, lb, ub);
        upper_diagonals[i] = psmc_prob(rho, t, curr_intervals[i]->lb, curr_intervals[i]->ub)/base;
        assert(!isnan(upper_diagonals[i]) and upper_diagonals[i] >= 0);
    }
}

void TSP_smc::compute_lower_sums() {
    lower_sums[0] = 0;
    float factor;
    for (int i = 1; i < curr_intervals.size(); i++) {
        if (lower_sums[i-1] == 0) {
            lower_sums[i] = upper_diagonals[i]*curr_intervals[i-1]->get_prob();
        } else {
            factor = (exp(-curr_intervals[i]->lb) - exp(-curr_intervals[i]->ub))/(exp(-curr_intervals[i-1]->lb) - exp(-curr_intervals[i-1]->ub));
            lower_sums[i] = upper_diagonals[i]*curr_intervals[i-1]->get_prob() + factor*lower_sums[i-1];
        }
        assert(!isnan(lower_sums[i]));
    }
}

void TSP_smc::compute_upper_sums() {
    int n = (int) curr_intervals.size();
    upper_sums[n-1] = 0;
    for (int i = n-2; i >=0 ; i--) {
        upper_sums[i] = upper_sums[i+1] + curr_intervals[i+1]->get_prob();
        assert(!isnan(upper_sums[i]));
    }
}

void TSP_smc::sanity_check(Recombination r) {
    for (Interval *i : curr_intervals) {
        if (i->lb == i->ub and i->lb == r.inserted_node->time and i->branch != r.target_branch) {
            i->update_prob(0.0);
        }
    }
}

vector<Interval *> TSP_smc::get_state_space(int x) {
    map<int, vector<Interval *>>::iterator state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->second;
}

int TSP_smc::get_prev_breakpoint(int x) {
    map<int, vector<Interval *>>::iterator state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->first;
}

void TSP_smc::set_interval_constraint(Recombination &r) {
    for (Interval *i : curr_intervals) {
        if (i->ub < r.start_time) {
            i->update_prob(0.0);
        } else {
            i->lb = max(r.start_time, i->lb);
            i->fill_time();
        }
    }
}

void TSP_smc::set_point_constraint(Recombination &r) {
    float t = r.inserted_node->time;
    if (t == numeric_limits<float>::infinity()) {
        return;
    }
    Interval *point_interval = search_point_interval(r);
    point_interval->node = r.inserted_node;
    for (Interval *i : curr_intervals) {
        if (i == point_interval) {
            i->update_prob(1.0);
        } else {
            i->update_prob(0.0);
        }
    }
}

Interval *TSP_smc::search_point_interval(Recombination r) {
    float t = r.inserted_node->time;
    Interval *point_interval = nullptr;
    for (Interval *i : curr_intervals) {
        if (i->lb < t and i->ub > t) {
            point_interval = i;
        }
    }
    for (Interval *i : curr_intervals) {
        if (i->ub == i->lb and i->lb == t) {
            point_interval = i;
        }
    }
    if (point_interval != nullptr) {
        return point_interval;
    }
    vector<Interval *> candidate_point_intervals = {};
    for (Interval *i : curr_intervals) {
        if (i->lb <= t and i->ub >= t) {
            candidate_point_intervals.push_back(i);
        }
    }
    assert(candidate_point_intervals.size() == 2);
    Interval *test_interval = candidate_point_intervals[0];
    while (test_interval->source_intervals.size() > 0) {
        test_interval = test_interval->source_intervals[0];
        if (test_interval->branch.upper_node == r.inserted_node or test_interval->branch.lower_node == r.inserted_node) {
            return candidate_point_intervals[1];
        }
    }
    return candidate_point_intervals[0];
}
