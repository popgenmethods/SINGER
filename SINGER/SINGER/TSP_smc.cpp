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
}

void TSP_smc::set_gap(float q) {
    gap = q;
}

void TSP_smc::set_emission(shared_ptr<Emission> e) {
    eh = e;
}

void TSP_smc::set_check_points(set<float> &p) {
    check_points = p;
}

void TSP_smc::reserve_memory(int length) {
    forward_probs.reserve(length);
    rhos.reserve(length);
}

void TSP_smc::start(Branch &branch, float t) {
    cut_time = t;
    curr_index = 0;
    curr_branch = branch;
    lower_bound = max(cut_time, branch.lower_node->time);
    generate_intervals(branch, branch.lower_node->time, branch.upper_node->time);
    set_dimensions();
    compute_factors();
    for (int i = 0; i < curr_intervals.size(); i++) {
        forward_probs[0][i] = exp(-curr_intervals[i].lb) - exp(-curr_intervals[i].ub);
    }
    state_spaces[0] = curr_intervals;
}

void TSP_smc::transfer(Recombination &r, Branch &prev_branch, Branch &next_branch) {
    rhos.push_back(0);
    prev_rho = -1;
    prev_theta = -1;
    prev_node = nullptr;
    curr_index += 1;
    curr_branch = next_branch;
    sanity_check(r);
    lower_bound = max(cut_time, next_branch.lower_node->time);
    if (prev_branch == r.source_branch and next_branch == r.merging_branch) {
        set_interval_constraint(r);
    } else if (prev_branch == r.target_branch and next_branch == r.recombined_branch) {
        set_point_constraint(r);
    }
    curr_intervals.clear();
    if (prev_branch == r.source_branch and next_branch == r.merging_branch) { // switch to a point mass
        float t = r.deleted_node->time;
        generate_intervals(next_branch, next_branch.lower_node->time, t);
        generate_intervals(next_branch, t, t);
        forward_probs[curr_index].back() = 1.0f; // set point mass prob marker
        curr_intervals.back().node = r.deleted_node; // set point mass node marker
        generate_intervals(next_branch, t, next_branch.upper_node->time);
    } else if (prev_branch == r.target_branch and next_branch == r.recombined_branch) { // switch from a point mass
        generate_intervals(next_branch, next_branch.lower_node->time, r.start_time);
        generate_intervals(next_branch, r.start_time, next_branch.upper_node->time);
        for (int i = 0; i < curr_intervals.size(); i++) {
            if (curr_intervals[i].time >= r.start_time) {
                forward_probs[curr_index][i] = 1.0;
            }
        }
    } else {
        float lb;
        float ub;
        lb = next_branch.lower_node->time;
        ub = max(prev_branch.lower_node->time, next_branch.lower_node->time);
        generate_intervals(next_branch, lb, ub);
        transfer_intervals(r, prev_branch, next_branch);
        lb = min(curr_intervals.back().ub, next_branch.upper_node->time);
        ub = next_branch.upper_node->time;
        generate_intervals(next_branch, lb, ub);
    }
    state_spaces[curr_index] = curr_intervals;
    set_dimensions();
}

void TSP_smc::recombine(Branch &prev_branch, Branch &next_branch) {
    assert(next_branch != Branch());
    vector<Interval> prev_intervals = curr_intervals;
    curr_intervals.clear();
    rhos.push_back(0);
    prev_rho = -1;
    curr_index += 1;
    lower_bound = max(cut_time, next_branch.lower_node->time);
    generate_intervals(next_branch, next_branch.lower_node->time, next_branch.upper_node->time);
    state_spaces[curr_index] = curr_intervals;
    set_dimensions();
    float new_prob;
    float base;
    for (int i = 0; i < prev_intervals.size(); i++) {
        base = recomb_prob(prev_intervals[i].time, curr_intervals.front().lb, curr_intervals.back().ub);
        for (int j = 0; j < curr_intervals.size(); j++) {
            if (base == 0) {
                new_prob = 1;
            } else {
                new_prob = recomb_prob(prev_intervals[i].time, curr_intervals[j].lb, curr_intervals[j].ub)*forward_probs[curr_index-1][i]/base;
            }
            forward_probs[curr_index].push_back(j);
        }
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
    compute_diagonals(rho);
    compute_lower_diagonals(rho);
    compute_upper_diagonals(rho);
    compute_lower_sums();
    compute_upper_sums();
    curr_index += 1;
    prev_rho = rho;
    forward_probs.push_back(lower_sums);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] += diagonals[i]*forward_probs[curr_index-1][i] + lower_diagonals[i]*upper_sums[i];
    }
}

void TSP_smc::null_emit(float theta, Node *query_node) {
    compute_null_emit_probs(theta, query_node);
    prev_theta = theta;
    prev_node = query_node;
    float ws = 0;
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] *= null_emit_probs[i];
        ws += forward_probs[curr_index][i];
    }
    assert(ws > 0);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] /= ws;
    }
}

void TSP_smc::mut_emit(float theta, float bin_size, set<float> &mut_set, Node *query_node) {
    compute_mut_emit_probs(theta, bin_size, mut_set, query_node);
    float ws = 0;
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] *= mut_emit_probs[i];
        ws += forward_probs[curr_index][i];
    }
    assert(ws > 0);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] /= ws;
    }
}

map<float, Node *> TSP_smc::sample_joining_nodes(int start_index, vector<float> &coordinates) {
    prev_rho = -1;
    map<float, Node *> joining_nodes = {};
    int x = curr_index;
    float pos = coordinates[x + start_index + 1];
    Interval interval = sample_curr_interval(x);
    Node *n = sample_joining_node(interval);
    joining_nodes[pos] = n;
    while (x >= 0) {
        x = trace_back_helper(interval, x);
        pos = coordinates[x];
        joining_nodes[pos] = n;
        if (x == 0) {
            break;
        } if (x == interval.start_pos) {
            if (source_interval.count(interval) > 0) {
                x -= 1;
                interval = sample_source_interval(interval);
            } else {
                x -= 1;
                interval = sample_recomb_interval(interval, x);
                n = sample_joining_node(interval);
            }
        } else {
            x -= 1;
            interval = sample_prev_interval(interval, x);
            n = sample_joining_node(interval);
        }
        prev_rho = -1;
    }
    return joining_nodes;
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
        l = 2*t - lower_bound - cut_time;
    } else {
        l = 2*s - lower_bound - cut_time;
    }
    if (l == 0) {
        pre_factor = rho;
    } else {
        pre_factor = (1 - exp(-rho*l))/l;
    }
    float integral;
    float cdf;
    if (t == cut_time and t == lower_bound) {
        return 0;
    } else if (t <= s) {
        integral = 2*t + exp(-t)*(exp(cut_time) + exp(lower_bound)) - cut_time - lower_bound - 2;
    } else {
        integral = 2*s + exp(cut_time - t) + exp(lower_bound - t) - 2*exp(s-t) - cut_time - lower_bound;
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
    assert(t1 >= lower_bound and s >= lower_bound);
    float prob;
    float base;
    float l = 2*s - lower_bound - cut_time;
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

void TSP_smc::generate_intervals(Branch &next_branch, float lb, float ub) {
    Interval new_interval = Interval();
    lb = max(cut_time, lb);
    ub = max(cut_time, ub);
    if (lb == ub) {
        if (lb == max(cut_time, next_branch.lower_node->time) or lb == next_branch.upper_node->time) {
            return;
        }
        else {
            new_interval = Interval(next_branch, lb, ub, curr_index);
            new_interval.fill_time();
            new_interval.hmm_index = (int) curr_intervals.size();
            curr_intervals.push_back(new_interval);
            temp.push_back(0);
            return;
        }
    }
    vector<float> points = generate_grid(lb, ub);
    float l;
    float u;
    for (int i = 0; i < points.size() - 1; i++) {
        l = points[i];
        u = points[i+1];
        new_interval = Interval(next_branch, l, u, curr_index);
        new_interval.fill_time();
        new_interval.hmm_index = (int) curr_intervals.size();
        curr_intervals.push_back(new_interval);
        temp.push_back(0);
    }
    forward_probs.push_back(temp);
    temp.clear();
}

void TSP_smc::transfer_intervals(Recombination &r, Branch &prev_branch, Branch &next_branch) {
    float lb;
    float ub;
    float p;
    vector<Interval> prev_intervals = get_state_space(curr_index);
    Interval interval = Interval();
    Interval new_interval = Interval();
    for (int i = 0; i < prev_intervals.size(); i++) {
        interval = prev_intervals[i];
        lb = max(interval.lb, next_branch.lower_node->time);
        ub = min(interval.ub, next_branch.upper_node->time);
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
            float w = get_prop(lb, ub, interval.lb, interval.ub);
            p = w*forward_probs[curr_index-1][i];
            new_interval = Interval(next_branch, lb, ub, curr_index);
            new_interval.fill_time();
            new_interval.node = interval.node;
            new_interval.hmm_index = (int) curr_intervals.size();
            curr_intervals.push_back(new_interval);
            forward_probs[curr_index].push_back(p);
        }
    }
}

void TSP_smc::set_dimensions() {
    dim = (int) curr_intervals.size();
    diagonals.resize(dim); diagonals.assign(dim, 0);
    lower_diagonals.resize(dim); lower_diagonals.assign(dim, 0);
    upper_diagonals.resize(dim); upper_diagonals.assign(dim, 0);
    lower_sums.resize(dim); lower_sums.assign(dim, 0);
    upper_sums.resize(dim); upper_sums.assign(dim, 0);
    null_emit_probs.resize(dim); null_emit_probs.assign(dim, 0);
    mut_emit_probs.resize(dim); mut_emit_probs.assign(dim, 0);
    factors.resize(dim); factors.assign(dim, 0);
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
        p = p1/p2;
    }
    return p;
}

void TSP_smc::compute_null_emit_probs(float theta, Node *query_node) {
    if (theta == prev_theta and query_node == prev_node) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        null_emit_probs[i] = eh->null_emit(curr_branch, curr_intervals[i].time, theta, query_node);
    }
}

void TSP_smc::compute_mut_emit_probs(float theta, float bin_size, set<float> &mut_set, Node *query_node) {
    for (int i = 0; i < dim; i++) {
        mut_emit_probs[i] = eh->mut_emit(curr_branch, curr_intervals[i].time, theta, bin_size, mut_set, query_node);
    }
}

void TSP_smc::compute_diagonals(float rho) {
    if (rho == prev_rho) {
        return;
    }
    float t;
    float base;
    float lb = curr_intervals.front().lb;
    float ub = curr_intervals.back().ub;
    Interval curr_interval = Interval();
    float diag = 0;
    for (int i = 0; i < dim; i++) {
        curr_interval = curr_intervals[i];
        t = curr_interval.time;
        base = psmc_prob(rho, curr_interval.time, lb, ub);
        diag = psmc_prob(rho, curr_interval.time, curr_interval.lb, curr_interval.ub);
        diagonals[i] = diag/base;
    }
}

void TSP_smc::compute_lower_diagonals(float rho) {
    if (rho == prev_rho) {
        return;
    }
    float t;
    float base;
    lower_diagonals[dim-1] = 0;
    float lb = max(cut_time, curr_intervals.front().lb);
    float ub = curr_intervals.back().ub;
    for (int i = 0; i < dim - 1; i++) {
        t = curr_intervals[i+1].time;
        base = psmc_prob(rho, t, lb, ub);
        lower_diagonals[i] = psmc_prob(rho, t, curr_intervals[i].lb, curr_intervals[i].ub)/base;
    }
}

void TSP_smc::compute_upper_diagonals(float rho) {
    if (rho == prev_rho) {
        return;
    }
    upper_diagonals[0] = 0;
    float lb = max(cut_time, curr_intervals.front().lb);
    float ub = curr_intervals.back().ub;
    float t;
    float base;
    for (int i = 1; i < dim; i++) {
        t = curr_intervals[i-1].time;
        base = psmc_prob(rho, t, lb, ub);
        upper_diagonals[i] = psmc_prob(rho, t, curr_intervals[i].lb, curr_intervals[i].ub)/base;
    }
}

void TSP_smc::compute_lower_sums() {
    for (int i = 0; i < dim; i++) {
        lower_sums[i] = upper_diagonals[i]*forward_probs[curr_index][i-1] + factors[i]*lower_sums[i-1];
    }
}

void TSP_smc::compute_upper_sums() {
    partial_sum(forward_probs[curr_index].rbegin(), forward_probs[curr_index].rend()-1, upper_sums.rbegin()+1);
}

void TSP_smc::compute_factors() {
    factors[0] = 0;
    for (int i = 1; i < dim; i++) {
        if (curr_intervals[i-1].lb == curr_intervals[i-1].ub) {
            factors[i] = 0;
        } else {
            factors[i] = (exp(-curr_intervals[i].lb) - exp(-curr_intervals[i].ub))/(exp(-curr_intervals[i-1].lb) - exp(-curr_intervals[i-1].ub));
        }
    }
}

void TSP_smc::compute_trace_back_probs(float rho, Interval &interval, vector<Interval> &intervals) {
    if (rho == prev_rho) {
        return;
    }
    for (int i = 0; i < trace_back_probs.size(); i++) {
        trace_back_probs[i] = psmc_prob(rho, intervals[i].time, interval.lb, interval.ub);
    }
}

void TSP_smc::sanity_check(Recombination &r) {
    for (int i = 0; i < curr_intervals.size(); i++) {
        if (curr_intervals[i].lb == curr_intervals[i].ub and curr_intervals[i].lb == r.inserted_node->time and curr_intervals[i].branch != r.target_branch) {
            forward_probs[curr_index-1][i] = 0; // curr_index - 1 because the index has already moved forward
        }
    }
}

vector<Interval> &TSP_smc::get_state_space(int x) {
    map<int, vector<Interval>>::iterator state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->second;
}

int TSP_smc::get_prev_breakpoint(int x) {
    map<int, vector<Interval>>::iterator state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->first;
}

Interval &TSP_smc::sample_curr_interval(int x) {
    vector<Interval> intervals = get_state_space(x);
    float ws = accumulate(forward_probs[x].begin(), forward_probs[x].end(), 0.0f);
    float q = random();
    float w = ws*q;
    for (int i = 0; i < intervals.size(); i++) {
        w -= forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            return curr_intervals[i];
        }
    }
    cerr << "Sample curr interval failed" << endl;
    exit(1);
}

Interval &TSP_smc::sample_prev_interval(Interval &interval, int x) {
    vector<Interval> intervals = get_state_space(x);
    lower_bound = intervals.front().lb;
    float ws = 0;
    float rho = rhos[x];
    compute_trace_back_probs(rho, interval, intervals);
    for (int i = 0; i < intervals.size(); i++) {
        if (intervals[i] != interval) {
            ws += trace_back_probs[i]*forward_probs[x][i];
        }
    }
    float q = random();
    float w = ws*q;
    for (int i = 0; i < intervals.size(); i++) {
        if (intervals[i] != interval) {
            w -= trace_back_probs[i]*forward_probs[x][i];
            if (w <= 0) {
                sample_index = i;
                return curr_intervals[i];
            }
        }
    }
    cerr << "sampling failed" << endl;
    exit(1);
}

Interval &TSP_smc::sample_source_interval(Interval &interval) {
    assert(source_interval.count(interval) > 0);
    return source_interval[interval];
}

Interval &TSP_smc::sample_recomb_interval(Interval &interval, int x) {
    vector<Interval> &intervals = get_state_space(x);
    float ws = 0;
    Interval prev_interval = Interval();
    for (int i = 0; i < intervals.size(); i++) {
        prev_interval = curr_intervals[i];
        ws += max(epsilon, recomb_prob(prev_interval.time, interval.lb, interval.ub))*forward_probs[x][i];
    }
    float q = random();
    float w = ws*q;
    for (int i = 0; i< intervals.size(); i++) {
        prev_interval = curr_intervals[i];
        w -= max(epsilon, recomb_prob(prev_interval.time, interval.lb, interval.ub))*forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            return curr_intervals[i];
        }
    }
    cerr << "sampling failed" << endl;
    exit(1);
}

int TSP_smc::trace_back_helper(Interval &interval, int x) {
    int y = get_prev_breakpoint(x);
    float non_recomb_prob = 0;
    float all_prob = 0;
    float q = random();
    float p = 1;
    float shrinkage;
    float rho;
    vector<Interval> intervals = get_state_space(x);
    lower_bound = intervals.front().lb;
    trace_back_probs = vector<float>(intervals.size());
    while (p > q and x > y) {
        rho = rhos[x-1];
        compute_trace_back_probs(rho, interval, intervals);
        prev_rho = rho;
        all_prob = inner_product(trace_back_probs.begin(), trace_back_probs.end(), forward_probs[x-1].begin(), 0.0);
        non_recomb_prob = trace_back_probs[sample_index]*forward_probs[x-1][sample_index];
        shrinkage = non_recomb_prob/all_prob;
        p *= shrinkage;
        if (p <= q) {
            return x;
        }
        x -= 1;
    }
    return x;
}

void TSP_smc::set_interval_constraint(Recombination &r) {
    for (int i = 0; i < curr_intervals.size(); i++) {
        if (curr_intervals[i].ub < r.start_time) {
            forward_probs[curr_index-1][i] = 0; // curr_index - 1 because the index has already moved forward
        } else {
            curr_intervals[i].lb = max(r.start_time, curr_intervals[i].lb);
            curr_intervals[i].fill_time();
        }
    }
}

void TSP_smc::set_point_constraint(Recombination &r) {
    Interval point_interval = search_point_interval(r);
    for (Interval &i : curr_intervals) {
        if (i == point_interval) {
            forward_probs[curr_index-1][i.hmm_index] = 1; // curr_index - 1 because the index has already moved forward
            i.node = r.inserted_node;
        } else {
            forward_probs[curr_index-1][i.hmm_index] = 0;
        }
    }
}

Interval TSP_smc::search_point_interval(Recombination &r) {
    float t = r.inserted_node->time;
    Interval point_interval = Interval();
    for (Interval &i : curr_intervals) {
        if (i.ub > t and i.lb < t) {
            point_interval = i;
        }
    }
    for (Interval &i : curr_intervals) {
        if (i.ub == i.lb and i.lb == t) {
            point_interval = i;
        }
    }
    if (point_interval != Interval()) {
        return point_interval;
    }
    vector<Interval> candidate_point_intervals = {};
    for (Interval &i : curr_intervals) {
        if (i.lb <= t and i.ub >= t) {
            candidate_point_intervals.push_back(i);
        }
    }
    assert(candidate_point_intervals.size() == 2);
    Interval test_interval = candidate_point_intervals[0];
    while (source_interval.count(test_interval) > 0) {
        test_interval = source_interval[test_interval];
        if (test_interval.branch.upper_node == r.inserted_node or test_interval.branch.lower_node == r.inserted_node) {
            return candidate_point_intervals[1];
        }
    }
    return candidate_point_intervals[0];
}

float TSP_smc::sample_time(float lb, float ub) {
    assert(lb <= ub);
    if (lb == ub) {
        return lb;
    }
    float delta = exp(-lb) - exp(-ub);
    float q;
    float p;
    float t = 0;
    if (lb > 5 and ub == numeric_limits<float>::infinity()) {
        t = lb + random()*log(2);
    }
    while (t <= lb or t >= ub) {
        q = random();
        p = 1 - exp(-lb) + q*delta;
        t = -log(1 - p);
    }
    assert(t != numeric_limits<float>::infinity());
    assert(lb == ub or (t > lb and t < ub));
    assert(t != 0);
    return t;
}

Node *TSP_smc::sample_joining_node(Interval &interval) {
    Node *n = nullptr;
    float t;
    if (interval.node != nullptr) {
        n = interval.node;
        interval.node = nullptr;
    } else {
        t = sample_time(interval.lb, interval.ub);
        n = new Node(t);
    }
    assert(n != nullptr);
    return n;
}
