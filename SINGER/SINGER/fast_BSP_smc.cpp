//
//  fast_BSP_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 4/21/23.
//

#include "fast_BSP_smc.hpp"

fast_BSP_smc::fast_BSP_smc() {}

fast_BSP_smc::~fast_BSP_smc() {
    for (auto &x : state_spaces) {
        for (Interval *interval : x.second) {
            delete interval;
        }
    }
    vector<vector<float>>().swap(forward_probs);
    map<Interval *, vector<float>>().swap(source_weights);
    map<Interval *, vector<Interval *>>().swap(source_intervals);
    map<int, vector<Interval *>>().swap(state_spaces);
}

void fast_BSP_smc::reserve_memory(int length) {
    forward_probs.reserve(length);
}

void fast_BSP_smc::start(set<Branch> &branches, float t) {
    cut_time = t;
    curr_index = 0;
    valid_branches = branches;
    float lb = 0;
    float ub = 0;
    float p = 0;
    Interval *new_interval = nullptr;
    cc = make_shared<Coalescent_calculator>(cut_time);
    cc->compute(valid_branches);
    for (Branch b : branches) {
        if (b.upper_node->time > cut_time) {
            lb = max(b.lower_node->time, cut_time);
            ub = b.upper_node->time;
            p = cc->weight(lb, ub);
            new_interval = new Interval(b, lb, ub, curr_index);
            curr_intervals.push_back(new_interval);
            temp_probs.push_back(p);
        }
    }
    forward_probs.push_back(temp_probs);
    compute_interval_info();
    weight_sums.push_back(0.0);
    set_dimensions();
    state_spaces[curr_index] = curr_intervals;
    temp_probs.clear();
}

void fast_BSP_smc::set_cutoff(float x) {
    cutoff = x;
}

void fast_BSP_smc::set_emission(shared_ptr<Emission> e) {
    eh = e;
}

void fast_BSP_smc::set_check_points(set<float> &p) {
    check_points = p;
}

void fast_BSP_smc::forward(float rho) {
    if (states_change) {
        update(rho);
    } else {
        regular_forward(rho);
    }
}


void fast_BSP_smc::transfer(Recombination &r) {
    rhos.push_back(0);
    prev_rho = -1;
    prev_theta = -1;
    recomb_sums.push_back(0);
    weight_sums.push_back(0);
    sanity_check(r);
    curr_index += 1;
    transfer_weights.clear();
    transfer_intervals.clear();
    temp_probs.clear();
    temp_intervals.clear();
    covered_branches.clear();
    for (int i = 0; i < curr_intervals.size(); i++) {
        process_interval(r, i);
    }
    generate_intervals(r);
    set_dimensions();
    state_spaces[curr_index] = curr_intervals;
    states_change = false;
}

void fast_BSP_smc::regular_forward(float rho) {
    rhos.push_back(rho);
    compute_recomb_probs(rho);
    compute_recomb_weights(rho);
    prev_rho = rho;
    curr_index += 1;
    recomb_sum = inner_product(recomb_probs.begin(), recomb_probs.end(), forward_probs[curr_index - 1].begin(), 0.0);
    forward_probs.push_back(recomb_probs);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] = forward_probs[curr_index - 1][i]*(1 - recomb_probs[i]) + recomb_sum*recomb_weights[i];
    }
    recomb_sums.push_back(recomb_sum);
    weight_sums.push_back(weight_sum);
}

void fast_BSP_smc::update(float rho) {
    float lb, ub;
    Interval *prev_interval, *new_interval;
    Branch prev_branch;
    rhos.push_back(rho);
    compute_recomb_probs(rho);
    prev_rho = -1;
    prev_theta = -1;
    curr_index += 1;
    recomb_sum = inner_product(recomb_probs.begin(), recomb_probs.end(), forward_probs[curr_index - 1].begin(), 0.0);
    temp_probs.clear();
    temp_intervals.clear();
    covered_branches.clear();
    for (int i = 0; i < dim; i++) {
        prev_interval = curr_intervals[i];
        prev_branch = prev_interval->branch;
        if (valid_branches.count(prev_branch) > 0) {
            new_interval = duplicate_interval(prev_interval);
            temp_intervals.push_back(new_interval);
            temp_probs.push_back(forward_probs[curr_index - 1][i]);
            covered_branches.insert(prev_branch);
        }
    }
    for (Branch b : valid_branches) {
        if (covered_branches.count(b) == 0) {
            lb = max(cut_time, b.lower_node->time);
            ub = b.upper_node->time;
            new_interval = new Interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp_probs.push_back(0);
        }
    }
    forward_probs.push_back(temp_probs);
    curr_intervals = temp_intervals;
    state_spaces[curr_index] = curr_intervals;
    compute_interval_info();
    set_dimensions();
    compute_recomb_probs(rho);
    compute_recomb_weights(rho);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] += recomb_sum*recomb_weights[i];
    }
    recomb_sums.push_back(recomb_sum);
    weight_sums.push_back(weight_sum);
    states_change = false;
}

float fast_BSP_smc::get_recomb_prob(float rho, float t) {
    float p = rho*(t - cut_time)*exp(-rho*(t - cut_time));
    return p;
}

void fast_BSP_smc::null_emit(float theta, Node *query_node) {
    compute_null_emit_prob(theta, query_node);
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

void fast_BSP_smc::mut_emit(float theta, float bin_size, set<float> &mut_set, Node *query_node) {
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

map<float, Branch> fast_BSP_smc::sample_joining_branches(int start_index, vector<float> &coordinates) {
    prev_rho = -1;
    map<float, Branch> joining_branches = {};
    int x = curr_index;
    float pos = coordinates[x + start_index + 1];
    Interval *interval = sample_curr_interval(x);
    Branch b = interval->branch;
    joining_branches[pos] = b;
    while (x >= 0) {
        x = trace_back_helper(interval, x);
        b = interval->branch;
        pos = coordinates[x + start_index];
        joining_branches[pos] = b;
        if (x == 0) {
            break;
        } else if (x == interval->start_pos) {
            x -= 1;
            if (rhos[x] == 0) {
                interval = sample_source_interval(interval, x);
                b = interval->branch;
            } else {
                interval = sample_connection_interval(interval, x);
                b = interval->branch;
            }
        } else {
            x -= 1;
            interval = sample_prev_interval(x);
            b = interval->branch;
        }
    }
    simplify(joining_branches);
    return joining_branches;
}

void fast_BSP_smc::write_forward_probs(string filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open the file." << std::endl;
    }

    // Write the vector<vector<float>> to the file
    for (const auto &row : forward_probs) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ' '; // separate values with a space (or use another delimiter if you prefer)
            }
        }
        file << '\n'; // add a newline character for each row
    }

    // Close the file
    file.close();
}

// private methods:

void fast_BSP_smc::update_states(set<Branch> &deletions, set<Branch> &insertions) {
    for (Branch b : deletions) {
        assert(valid_branches.count(b) > 0);
        valid_branches.erase(b);
    }
    for (Branch b : insertions) {
        valid_branches.insert(b);
    }
    states_change = true;
}

void fast_BSP_smc::set_dimensions() {
    dim = (int) curr_intervals.size();
    recomb_probs.resize(dim); recomb_probs.assign(dim, 0);
    recomb_weights.resize(dim); recomb_weights.assign(dim, 0);
    null_emit_probs.resize(dim); null_emit_probs.assign(dim, 0);
    mut_emit_probs.resize(dim); mut_emit_probs.assign(dim, 0);
}

void fast_BSP_smc::compute_recomb_probs(float rho) {
    if (prev_rho == rho) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        recomb_probs[i] = get_recomb_prob(rho, curr_intervals[i]->time);
    }
}

void fast_BSP_smc::compute_recomb_weights(float rho) {
    if (prev_rho == rho) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        if (curr_intervals[i]->full(cut_time)) {
            recomb_weights[i] = recomb_probs[i]*curr_intervals[i]->weight;
        }
    }
    weight_sum = accumulate(recomb_weights.begin(), recomb_weights.end(), 0.0);
    for (int i = 0; i < dim; i++) {
        recomb_weights[i] /= weight_sum;
    }
}

void fast_BSP_smc::compute_null_emit_prob(float theta, Node *query_node) {
    if (theta == prev_theta and query_node == prev_node) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        null_emit_probs[i] = eh->null_emit(curr_intervals[i]->branch, curr_intervals[i]->time, theta, query_node);
    }
}

void fast_BSP_smc::compute_mut_emit_probs(float theta, float bin_size, set<float> &mut_set, Node *query_node) {
    for (int i = 0; i < dim; i++) {
        mut_emit_probs[i] = eh->mut_emit(curr_intervals[i]->branch, curr_intervals[i]->time, theta, bin_size, mut_set, query_node);
    }
}

void fast_BSP_smc::transfer_helper(Interval_info next_interval, Interval *prev_interval, float w) {
    if (valid_branches.count(next_interval.branch) == 0) {
        return;
    }
    if (transfer_weights.count(next_interval) > 0) {
        transfer_weights[next_interval].push_back(w);
        transfer_intervals[next_interval].push_back(prev_interval);
    } else {
        transfer_weights[next_interval] = {w};
        transfer_intervals[next_interval] = {prev_interval};
    }
}

void fast_BSP_smc::transfer_helper(Interval_info next_interval) {
    if (valid_branches.count(next_interval.branch) == 0) {
        return;
    }
    if (transfer_weights.count(next_interval) == 0) {
        transfer_weights[next_interval] = {};
        transfer_intervals[next_interval] = {};
    }
}

float fast_BSP_smc::compute_transfer_prob() {
    float p = 0;
    for (auto &x : transfer_weights) {
        p += accumulate(x.second.begin(), x.second.end(), 0.0f);
    }
    return p;
}

Interval *fast_BSP_smc::duplicate_interval(Interval *interval) {
    Interval *new_interval = new Interval(interval->branch, interval->lb, interval->ub, curr_index);
    // source_intervals[new_interval] = {interval};
    // source_weights[new_interval] = {1};
    return new_interval;
}

void fast_BSP_smc::add_new_branches(Recombination &r) { // add recombined branch and merging branch, if legal
    Interval_info next_interval;
    float lb = 0;
    float ub = 0;
    if (r.merging_branch != Branch() and r.merging_branch.upper_node->time > cut_time) {
        lb = max(cut_time, r.merging_branch.lower_node->time);
        ub = r.merging_branch.upper_node->time;
        next_interval = Interval_info(r.merging_branch, lb, ub);
        transfer_helper(next_interval);
    }
    if (r.recombined_branch != Branch() and r.recombined_branch.upper_node->time > cut_time) {
        lb = max(cut_time, r.recombined_branch.lower_node->time);
        ub = r.recombined_branch.upper_node->time;
        next_interval = Interval_info(r.recombined_branch, lb, ub);
        transfer_helper(next_interval);
    }
}

void fast_BSP_smc::compute_interval_info() {
    cc->compute(valid_branches);
    float t;
    float p;
    for (Interval *i : curr_intervals) {
        p = cc->weight(i->lb, i->ub);
        t = cc->time(i->lb, i->ub);
        i->assign_weight(p);
        i->assign_time(t);
    }
}

void fast_BSP_smc::sanity_check(Recombination &r) {
}

void fast_BSP_smc::generate_intervals(Recombination &r) {
    float p0 = compute_transfer_prob();
    Branch b;
    float lb;
    float ub;
    float p;
    vector<Interval *> intervals;
    vector<float> weights;
    Interval_info interval;
    Interval *new_interval;
    for (auto x : transfer_weights) {
        interval = x.first;
        weights = x.second;
        intervals = transfer_intervals[x.first];
        b = interval.branch;
        lb = interval.lb;
        ub = interval.ub;
        p = accumulate(weights.begin(), weights.end(), 0.0);
        assert(!isnan(p));
        if (lb == max(cut_time, b.lower_node->time) and ub == b.upper_node->time) { // full intervals
            new_interval = new Interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp_probs.push_back(p);
            if (weights.size() > 0) {
                source_weights[new_interval] = weights;
                source_intervals[new_interval] = intervals;
            }
            covered_branches.insert(b);
        } else if (p >= cutoff*p0) { // partial intervals
            new_interval = new Interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp_probs.push_back(p);
            if (weights.size() > 0) {
                source_weights[new_interval] = weights;
                source_intervals[new_interval] = intervals;
            }
        }
    }
    for (Branch b : valid_branches) {
        if (covered_branches.count(b) == 0) {
            lb = max(cut_time, b.lower_node->time);
            ub = b.upper_node->time;
            new_interval = new Interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp_probs.push_back(0);
        }
    }
    forward_probs.push_back(temp_probs);
    curr_intervals = temp_intervals;
    compute_interval_info();
}

float fast_BSP_smc::get_overwrite_prob(Recombination &r, float lb, float ub) {
    if (check_points.count(curr_index) > 0) {
        return 0.0;
    }
    float join_time = r.inserted_node->time;
    float p1 = cc->weight(lb, ub);
    float p2 = cc->weight(max(cut_time, r.start_time), join_time);
    if (p1 == 0 and p2 == 0) {
        return 1.0;
    }
    float overwrite_prob = p2/(p1 + p2);
    return overwrite_prob;
}

/*
void fast_BSP_smc::process_interval(Recombination &r, int i) {
    Branch prev_branch = curr_intervals[i]->branch;
    if (prev_branch == r.source_branch) {
        process_source_interval(r, i);
    } else if (prev_branch == r.target_branch) {
        process_target_interval(r, i);
    } else {
        process_other_interval(r, i);
    }
}
 */

void fast_BSP_smc::process_interval(Recombination &r, int i) {
    Interval *prev_interval = curr_intervals[i];
    Interval *next_interval;
    Branch prev_branch = prev_interval->branch;
    /*
    if (!r.affect(prev_branch)) {
        next_interval = duplicate_interval(prev_interval);
        temp_intervals.push_back(next_interval);
        temp_probs.push_back(forward_probs[curr_index - 1][i]);
        source_intervals[next_interval] = {prev_interval};
        source_weights[next_interval] = {1};
        return;
    }
     */
    if (prev_branch == r.source_branch) {
        process_source_interval(r, i);
    } else if (prev_branch == r.target_branch) {
        process_target_interval(r, i);
    } else {
        process_other_interval(r, i);
    }
}

void fast_BSP_smc::process_source_interval(Recombination &r, int i) {
    float w1, w2, lb, ub = 0;
    Interval *prev_interval = curr_intervals[i];
    float p = forward_probs[curr_index - 1][i];
    float point_time = r.source_branch.upper_node->time;
    float break_time = r.start_time;
    Branch next_branch;
    Interval_info next_interval;
    if (prev_interval->ub <= break_time) {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        next_branch = r.recombined_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else if (prev_interval->lb >= break_time) {
        lb = point_time;
        ub = point_time;
        next_branch = r.merging_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else {
        w1 = cc->weight(prev_interval->lb, break_time);
        w2 = cc->weight(break_time, prev_interval->ub);
        if (w1 == 0 and w2 == 0) {
            w1 = 1;
            w2 = 0;
        } else {
            w1 = w1/(w1 + w2);
            w2 = 1 - w1;
        }
        lb = prev_interval->lb;
        ub = break_time;
        next_branch = r.recombined_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w1*p);
        lb = point_time;
        ub = point_time;
        next_branch = r.merging_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w2*p);
    }
}

void fast_BSP_smc::process_target_interval(Recombination &r, int i) {
    float w0, w1, w2, lb, ub = 0;
    Interval *prev_interval = curr_intervals[i];
    float p = forward_probs[curr_index - 1][i];
    float join_time = r.inserted_node->time;
    Branch next_branch;
    Interval_info next_interval;
    if (prev_interval->lb == prev_interval->ub and prev_interval->lb == join_time) {
        lb = max(cut_time, r.start_time);
        ub = r.recombined_branch.upper_node->time;
        next_branch = r.recombined_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else if (prev_interval->lb >= join_time) {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        next_branch = r.upper_transfer_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else if (prev_interval->ub <= join_time) {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        next_branch = r.lower_transfer_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else {
        w0 = get_overwrite_prob(r, prev_interval->lb, prev_interval->ub);
        w1 = cc->weight(prev_interval->lb, join_time);
        w2 = cc->weight(join_time, prev_interval->ub);
        if (w1 + w2 == 0) {
            w1 = 0;
            w2 = 0;
            w0 = 1.0;
        } else {
            w1 = w1/(w1 + w2);
            w2 = 1 - w1;
            w1 *= 1 - w0;
            w2 *= 1 - w0;
        }
        lb = prev_interval->lb;
        ub = join_time;
        next_branch = r.lower_transfer_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w1*p);
        lb = join_time;
        ub = prev_interval->ub;
        next_branch = r.upper_transfer_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w2*p);
        lb = max(r.start_time, cut_time);
        ub = r.recombined_branch.upper_node->time;;
        next_branch = r.recombined_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w0*p);
    }
}

void fast_BSP_smc::process_other_interval(Recombination &r, int i) {
    float lb, ub = 0;
    Interval *prev_interval = curr_intervals[i];
    float p = forward_probs[curr_index - 1][i];
    Branch next_branch;
    Interval_info next_interval;
    if (r.affect(prev_interval->branch)) {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        next_branch = r.merging_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        next_branch = prev_interval->branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    }
}

float fast_BSP_smc::random() {
    return (float) rand()/RAND_MAX;
}

int fast_BSP_smc::get_prev_breakpoint(int x) {
    auto state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->first;
}

vector<Interval *> &fast_BSP_smc::get_state_space(int x) {
    auto state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->second;
}

int fast_BSP_smc::get_interval_index(Interval *interval, vector<Interval *> &intervals) {
    auto it = find(intervals.begin(), intervals.end(), interval);
    // assert(it != intervals.end());
    int index = (int) distance(intervals.begin(), it);
    // assert(intervals[index] == interval);
    return index;
}
 
void fast_BSP_smc::simplify(map<float, Branch> &joining_branches) {
    map<float, Branch> simplified_joining_branches = {};
    Branch curr_branch = joining_branches.begin()->second;
    simplified_joining_branches[joining_branches.begin()->first] = curr_branch;
    for (auto x : joining_branches) {
        if (x.second != curr_branch) {
            simplified_joining_branches.insert(x);
            curr_branch = x.second;
        }
    }
    simplified_joining_branches[joining_branches.rbegin()->first] = joining_branches.rbegin()->second;
    joining_branches = simplified_joining_branches;
}

Interval *fast_BSP_smc::sample_curr_interval(int x) {
    vector<Interval *> &intervals = get_state_space(x);
    float ws = accumulate(forward_probs[x].begin(), forward_probs[x].end(), 0.0);
    float q = random();
    float w = ws*q;
    for (int i = 0; i < intervals.size(); i++) {
        w -= forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            return intervals[i];
        }
    }
    cerr << "bsp curr sampling failed" << endl;
    exit(1);
}

Interval *fast_BSP_smc::sample_prev_interval(int x) {
    vector<Interval *> &intervals = get_state_space(x);
    float rho = rhos[x];
    float ws = recomb_sums[x];
    float q = random();
    float w = ws*q;
    for (int i = 0; i < intervals.size(); i++) {
        w -= get_recomb_prob(rho, intervals[i]->time)*forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            return intervals[i];
        }
    }
    cerr << "bsp prev sampling failed" << endl;
    exit(1);
}

Interval *fast_BSP_smc::sample_source_interval(Interval *interval, int x) {
    assert(source_intervals.count(interval) > 0);
    vector<Interval *> &prev_intervals = get_state_space(x);
    vector<Interval *> intervals = source_intervals[interval];
    vector<float> weights = source_weights[interval];
    float q = random();
    float ws = accumulate(weights.begin(), weights.end(), 0.0);
    float w = ws*q;
    for (int i = 0; i < weights.size(); i++) {
        w -= weights[i];
        if (w <= 0) {
            sample_index = get_interval_index(intervals[i], prev_intervals);
            return intervals[i];
        }
    }
    cerr << "fbsp sample source interval failed" << endl;
    exit(1);
}

Interval *fast_BSP_smc::sample_connection_interval(Interval *interval, int x) {
    vector<Interval *> &prev_intervals = get_state_space(x);
    int n = (int) prev_intervals.size();
    float source_recomb_prob, target_proportion;
    vector<float> weights = vector<float>(n);
    recomb_sum = recomb_sums[x];
    weight_sum = weight_sums[x + 1];
    if (interval->full(cut_time)) {
        target_proportion = interval->weight*get_recomb_prob(rhos[x], interval->time)/weight_sum;
    } else {
        target_proportion = 0;
    }
    for (int i = 0; i < n; i++) {
        source_recomb_prob = forward_probs[x][i]*get_recomb_prob(rhos[x], interval->time); // probability that goes out from i-state
        if (ei(interval, prev_intervals[i])) {
            weights[i] += forward_probs[x][i] - source_recomb_prob;
        }
        weights[i] += source_recomb_prob*target_proportion;
    }
    float q = random();
    float ws = accumulate(weights.begin(), weights.end(), 0.0);
    float w = ws*q;
    for (int i = 0; i < weights.size(); i++) {
        w -= weights[i];
        if (w <= 0) {
            sample_index = i;
            // assert(ei(prev_intervals[i], interval));
            return prev_intervals[i];
        }
    }
    cerr << "fbsp sample connection interval failed" << endl;
    exit(1);
    /*
    for (int i = 0; i < prev_intervals.size(); i++) {
        if (ei(prev_intervals[i], interval)) {
            sample_index = i;
            return prev_intervals[i];
        }
    }
    cerr << "fbsp sample connection interval failed" << endl;
    exit(1);
     */
}

int fast_BSP_smc::trace_back_helper(Interval *interval, int x) {
    // int y = get_prev_breakpoint(x);
    int y = interval->start_pos;
    if (!interval->full(cut_time)) {
        // return interval->start_pos;
        return y;
    }
    float p = random();
    float q = 1;
    float shrinkage = 0;
    float recomb_prob = 0;
    float non_recomb_prob = 0;
    float all_prob = 0;
    while (x > y) {
        recomb_sum = recomb_sums[x - 1];
        weight_sum = weight_sums[x];
        if (recomb_sum == 0) {
            shrinkage = 1;
        } else {
            recomb_prob = get_recomb_prob(rhos[x - 1], interval->time);
            non_recomb_prob = (1 - recomb_prob)*forward_probs[x - 1][sample_index];
            all_prob = non_recomb_prob + recomb_sum*interval->weight*recomb_prob/weight_sum;
            shrinkage = non_recomb_prob/all_prob;
            assert(shrinkage >= 0 and shrinkage <= 1);
        }
        q *= shrinkage;
        if (p >= q) {
            return x;
        }
        x -= 1;
    }
    return y;
}

void fast_BSP_smc::check_intervals() {
    vector<Interval *> intervals = curr_intervals;
    set<Branch> full_branches = {};
    for (Interval *i : intervals) {
        if (i->full(cut_time)) {
            full_branches.insert(i->branch);
        }
        assert(valid_branches.count(i->branch) > 0);
    }
    for (Branch b : full_branches) {
        assert(valid_branches.count(b) > 0);
    }
    for (Branch b : valid_branches) {
        assert(full_branches.count(b) > 0);
    }
    sort(intervals.begin(), intervals.end(), compare_interval());
    compare_interval comp;
    for (int i = 0; i < intervals.size() - 1; i++) {
        if (!comp(intervals[i], intervals[i+1]) && !comp(intervals[i+1], intervals[i])) {
            assert(false);
        }
    }
}

void fast_BSP_smc::check_recomb_sums() {
    float ws, rs= 0;
    float rho;
    for (int i = 0; i < rhos.size() - 1; i++) {
        vector<Interval *> intervals = get_state_space(i);
        vector<float> r = vector<float>(intervals.size());
        ws = 0.0;
        rho = rhos[i];
        for (int j = 0; j < intervals.size(); j++) {
            r[j] = get_recomb_prob(rho, intervals[j]->time);
        }
        ws = inner_product(r.begin(), r.end(), forward_probs[i].begin(), 0.0);
        rs = recomb_sums[i];
        assert(ws == rs);
    }
}
