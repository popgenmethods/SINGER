//
//  reduced_BSP.cpp
//  SINGER
//
//  Created by Yun Deng on 6/7/23.
//

#include "reduced_BSP.hpp"

reduced_BSP::reduced_BSP() {}

reduced_BSP::~reduced_BSP() {
    vector<vector<double>>().swap(forward_probs);
    map<int, vector<Interval_ptr>>().swap(state_spaces);
}

void reduced_BSP::reserve_memory(int length) {
    forward_probs.reserve(length);
}

void reduced_BSP::start(set<Branch> &start_branches, set<Interval_info> &start_intervals, double t) {
    cut_time = t;
    curr_index = 0;
    set<Interval_info> empty_set = {};
    update_states(empty_set, start_intervals);
    for (const Branch &b : start_branches) {
        if (b.upper_node->time > cut_time) {
            all_branches.insert(b);
        }
    }
    double lb = 0;
    double ub = 0;
    double p = 0;
    Interval_ptr new_interval = nullptr;
    cc = make_shared<Coalescent_calculator>(cut_time);
    cc->compute(all_branches);
    for (const Branch &b : all_branches) {
        if (b.upper_node->time > cut_time) {
            lb = max(b.lower_node->time, cut_time);
            ub = b.upper_node->time;
            p = cc->weight(lb, ub);
            new_interval = create_interval(b, lb, ub, curr_index);
            curr_intervals.push_back(new_interval);
            if (reduced_branches.count(b) > 0) {
                temp.push_back(p);
            } else {
                temp.push_back(0);
            }
        }
    }
    forward_probs.push_back(temp);
    weight_sums.push_back(0.0);
    reduced_sums.push_back(0.0);
    set_dimensions();
    compute_interval_info();
    state_spaces[curr_index] = curr_intervals;
    temp.clear();
}

void reduced_BSP::start(set<Branch> &branches, double t) {
    cut_time = t;
    curr_index = 0;
    for (Branch b : branches) {
        if (b.upper_node->time > cut_time) {
            all_branches.insert(b);
        }
    }
    double lb = 0;
    double ub = 0;
    double p = 0;
    Interval_ptr new_interval = nullptr;
    cc = make_shared<Coalescent_calculator>(cut_time);
    cc->compute(all_branches);
    for (Branch b : branches) {
        if (b.upper_node->time > cut_time) {
            lb = max(b.lower_node->time, cut_time);
            ub = b.upper_node->time;
            p = cc->weight(lb, ub);
            new_interval = create_interval(b, lb, ub, curr_index);
            curr_intervals.push_back(new_interval);
            temp.push_back(p);
        }
    }
    forward_probs.push_back(temp);
    weight_sums.push_back(0.0);
    reduced_sums.push_back(0.0);
    set_dimensions();
    compute_interval_info();
    state_spaces[curr_index] = curr_intervals;
    temp.clear();
}

void reduced_BSP::set_cutoff(double x) {
    cutoff = x;
}

void reduced_BSP::set_emission(shared_ptr<Emission> e) {
    eh = e;
}

void reduced_BSP::set_check_points(set<double> &p) {
    check_points = p;
}

void reduced_BSP::forward(double rho) {
    rhos.push_back(rho);
    compute_recomb_probs(rho);
    compute_recomb_weights(rho);
    compute_reduced_sums();
    prev_rho = rho;
    curr_index += 1;
    recomb_sum = inner_product(recomb_probs.begin(), recomb_probs.end(), forward_probs[curr_index - 1].begin(), 0.0);
    forward_probs.push_back(recomb_probs);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] = forward_probs[curr_index - 1][i]*(1 - recomb_probs[i]) + recomb_sum*recomb_weights[i];
    }
    for (int i = 0; i < dim; i++) {
        Interval_ptr interval = curr_intervals[i];
        if (reduced_branches.count(interval->branch) == 0) {
            forward_probs[curr_index][i] = 0;
        }
    }
    recomb_sums.push_back(recomb_sum);
    weight_sums.push_back(weight_sum);
    reduced_sums.push_back(reduced_sum);
}


void reduced_BSP::transfer(Recombination &r) {
    rhos.push_back(0);
    prev_rho = -1;
    prev_theta = -1;
    recomb_sums.push_back(0);
    weight_sums.push_back(0);
    reduced_sums.push_back(0);
    sanity_check(r);
    curr_index += 1;
    transfer_weights.clear();
    transfer_intervals.clear();
    temp.clear();
    temp_intervals.clear();
    update_states(r.deleted_branches, r.inserted_branches);
    for (int i = 0; i < curr_intervals.size(); i++) {
        process_interval(r, i);
    }
    add_new_branches(r);
    generate_intervals(r);
    set_dimensions();
    compute_interval_info();
    state_spaces[curr_index] = curr_intervals;
}

double reduced_BSP::get_recomb_prob(double rho, double t) {
    double p = rho*(t - cut_time)*exp(-rho*(t - cut_time));
    return p;
}

void reduced_BSP::null_emit(double theta, Node_ptr query_node) {
    compute_null_emit_prob(theta, query_node);
    prev_theta = theta;
    prev_node = query_node;
    double ws = 0;
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] *= null_emit_probs[i];
        // ws += forward_probs[curr_index][i];
    }
    ws = accumulate(forward_probs[curr_index].begin(), forward_probs[curr_index].end(), 0.0f);
    assert(ws > 0);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] /= ws;
    }
}

void reduced_BSP::mut_emit(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node) {
    compute_mut_emit_probs(theta, bin_size, mut_set, query_node);
    double ws = 0;
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] *= mut_emit_probs[i];
        ws += forward_probs[curr_index][i];
    }
    assert(ws > 0);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] /= ws;
    }
}

map<double, Branch> reduced_BSP::sample_joining_branches(int start_index, vector<double> &coordinates) {
    prev_rho = -1;
    map<double, Branch> joining_branches = {};
    int x = curr_index;
    int y = 0;
    double pos = coordinates[x + start_index + 1];
    Interval_ptr interval = sample_curr_interval(x);
    Branch b = interval->branch;
    joining_branches[pos] = b;
    while (x >= 0) {
        vector<Interval_ptr> &intervals = get_state_space(x);
        assert(intervals[sample_index] == interval);
        x = trace_back_helper(interval, x);
        b = interval->branch;
        pos = coordinates[x + start_index];
        joining_branches[pos] = b;
        y = get_prev_breakpoint(x);
        if (x == 0) {
            break;
        } else if (x == y) {
            x -= 1;
            interval = sample_source_interval(interval, x);
            b = interval->branch;
        } else {
            x -= 1;
            interval = sample_prev_interval(x);
            b = interval->branch;
        }
    }
    simplify(joining_branches);
    return joining_branches;
}

void reduced_BSP::write_forward_probs(string filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open the file." << std::endl;
    }

    // Write the vector<vector<double>> to the file
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

void reduced_BSP::update_states(set<Branch> &deletions, set<Branch> &insertions) {
    for (Branch b : deletions) {
        if (b.upper_node->time > cut_time) {
            assert(all_branches.count(b) > 0);
            all_branches.erase(b);
        }
        states_change = true;
    }
    for (Branch b : insertions) {
        if (b.upper_node->time > cut_time) {
            all_branches.insert(b);
        }
        states_change = true;
    }
}

void reduced_BSP::update_states(set<Interval_info> &deletions, set<Interval_info> &insertions) {
    for (const Interval_info &ii : deletions) {
        const Branch &b = ii.branch;
        set<Interval_info> &intervals = reduced_intervals[b];
        assert(intervals.count(ii) > 0);
        intervals.erase(ii);
        if (intervals.size() == 0) {
            reduced_intervals.erase(b);
            reduced_branches.erase(b);
        }
    }
    for (const Interval_info &ii : insertions) {
        const Branch b = ii.branch;
        if (reduced_branches.count(b) == 0) {
            reduced_branches.insert(b);
            reduced_intervals[b];
        }
        set<Interval_info> &intervals = reduced_intervals[b];
        intervals.insert(ii);
    }
    assert(reduced_branches.size() == reduced_intervals.size() and reduced_branches.size() > 0);
}

void reduced_BSP::set_dimensions() {
    dim = (int) curr_intervals.size();
    time_points.resize(dim); time_points.assign(dim, 0);
    raw_weights.resize(dim); raw_weights.assign(dim, 0);
    recomb_probs.resize(dim); recomb_probs.assign(dim, 0);
    recomb_weights.resize(dim); recomb_weights.assign(dim, 0);
    null_emit_probs.resize(dim); null_emit_probs.assign(dim, 0);
    mut_emit_probs.resize(dim); mut_emit_probs.assign(dim, 0);
}

void reduced_BSP::compute_recomb_probs(double rho) {
    if (prev_rho == rho) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        recomb_probs[i] = get_recomb_prob(rho, time_points[i]);
    }
}

void reduced_BSP::compute_recomb_weights(double rho) {
    if (prev_rho == rho) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        if (curr_intervals[i]->full(cut_time)) {
            recomb_weights[i] = recomb_probs[i]*raw_weights[i];
        }
    }
    weight_sum = accumulate(recomb_weights.begin(), recomb_weights.end(), 0.0);
    for (int i = 0; i < dim; i++) {
        recomb_weights[i] /= weight_sum;
    }
}

void reduced_BSP::compute_reduced_sums() {
    reduced_sum = 0;
    for (int i = 0; i < curr_intervals.size(); i++) {
        Interval_ptr interval = curr_intervals[i];
        const Branch &b = interval->branch;
        if (interval->full(cut_time) and reduced_branches.count(b) > 0) {
            reduced_sum += recomb_weights[i];
        }
    }
    assert(reduced_sum <= 1.001);
    reduced_sum = min(reduced_sum, 1.0);
}

void reduced_BSP::compute_null_emit_prob(double theta, Node_ptr query_node) {
    if (theta == prev_theta and query_node == prev_node) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        null_emit_probs[i] = eh->null_emit(curr_intervals[i]->branch, time_points[i], theta, query_node);
    }
}

void reduced_BSP::compute_mut_emit_probs(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node) {
    for (int i = 0; i < dim; i++) {
        mut_emit_probs[i] = eh->mut_emit(curr_intervals[i]->branch, time_points[i], theta, bin_size, mut_set, query_node);
    }
}

void reduced_BSP::transfer_helper(Interval_info &next_interval, Interval_ptr &prev_interval, double w) {
    transfer_weights[next_interval].push_back(w);
    transfer_intervals[next_interval].push_back(prev_interval);
}

void reduced_BSP::transfer_helper(Interval_info &next_interval) {
    transfer_weights[next_interval];
    transfer_intervals[next_interval];
}

void reduced_BSP::add_new_branches(Recombination &r) { // add recombined branch and merging branch, if legal
    Interval_info next_interval;
    double lb = 0;
    double ub = 0;
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

void reduced_BSP::compute_interval_info() {
    if (states_change) {
        cc->compute(all_branches);
    }
    states_change = false;
    double t;
    double p;
    for (int i = 0; i < curr_intervals.size(); i++) {
        Interval_ptr interval = curr_intervals[i];
        p = cc->weight(interval->lb, interval->ub);
        t = cc->time(interval->lb, interval->ub);
        raw_weights[i] = p;
        time_points[i] = t;
    }
    times[curr_index] = time_points;
    weights[curr_index] = raw_weights;
}

void reduced_BSP::sanity_check(Recombination &r) {
    for (int i = 0; i < curr_intervals.size(); i++) {
        Interval_ptr interval = curr_intervals[i];
        if (interval->lb == interval->ub and interval->lb == r.inserted_node->time and interval->branch != r.target_branch) {
            forward_probs[curr_index][i] = 0;
        }
    }
}

void reduced_BSP::get_full_branches(Recombination &r) {
    for (auto x : transfer_weights) {
        const Branch &b = x.first.branch;
        vector<double> &weights = x.second;
        double p = accumulate(weights.begin(), weights.end(), 0.0);
        if (p > 0 and x.first.lb == max(cut_time, b.lower_node->time) and x.first.ub == b.upper_node->time) {
            full_branches.insert(b);
        }
    }
    for (int i = 0; i < curr_intervals.size(); i++) {
        Interval_ptr interval = curr_intervals[i];
        double p = forward_probs[curr_index - 1][i];
        if (interval->full(cut_time) and !r.affect(interval->branch)) {
            if (p > 0) {
                full_branches.insert(interval->branch);
            }
        }
    }
}

void reduced_BSP::generate_intervals(Recombination &r) {
    Branch b;
    double lb;
    double ub;
    double p;
    vector<Interval_ptr > intervals;
    vector<double> weights;
    Interval_info interval;
    Interval_ptr new_interval = nullptr;
    auto y = transfer_intervals.begin();
    for (auto x = transfer_weights.begin(); x != transfer_weights.end(); ++x, ++y) {
        interval = x->first;
        const auto &weights = x->second;
        const auto &intervals = y->second;
        b = interval.branch;
        lb = interval.lb;
        ub = interval.ub;
        p = accumulate(weights.begin(), weights.end(), 0.0);
        assert(!isnan(p));
        if (lb == max(cut_time, b.lower_node->time) and ub == b.upper_node->time) { // full intervals
            new_interval = create_interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                new_interval->source_weights = move(weights);
                new_interval->intervals = move(intervals);
            }
        } else if (reduced_branches.count(b) > 0 or p >= cutoff) { // partial intervals
            new_interval = create_interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                new_interval->source_weights = move(weights);
                new_interval->intervals = move(intervals);
            }
        }
    }
    forward_probs.push_back(temp);
    curr_intervals = temp_intervals;
    for (int i = 0; i < curr_intervals.size(); i++) {
        if (reduced_branches.count(curr_intervals[i]->branch) == 0) {
            forward_probs[curr_index][i] = 0;
        }
    }
}

/*
void reduced_BSP::generate_intervals(Recombination &r) {
    full_branches.clear();
    get_full_branches(r);
    Branch b;
    double lb;
    double ub;
    double p;
    vector<Interval_ptr > intervals;
    vector<double> weights;
    Interval_info interval;
    Interval_ptr new_interval = nullptr;
    auto y = transfer_intervals.begin();
    for (auto x = transfer_weights.begin(); x != transfer_weights.end(); ++x, ++y) {
        interval = x->first;
        const auto &weights = x->second;
        const auto &intervals = y->second;
        b = interval.branch;
        lb = interval.lb;
        ub = interval.ub;
        p = accumulate(weights.begin(), weights.end(), 0.0);
        assert(!isnan(p));
        if (lb == max(cut_time, b.lower_node->time) and ub == b.upper_node->time) { // full intervals
            new_interval = create_interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                new_interval->source_weights = move(weights);
                new_interval->intervals = move(intervals);
            }
        } else if (reduced_branches.count(b) > 0) { // partial intervals
            new_interval = create_interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                new_interval->source_weights = move(weights);
                new_interval->intervals = move(intervals);
            }
        }
    }
    forward_probs.push_back(temp);
    curr_intervals = temp_intervals;
    for (int i = 0; i < curr_intervals.size(); i++) {
        if (reduced_branches.count(curr_intervals[i]->branch) == 0) {
            forward_probs[curr_index][i] = 0;
        }
    }
}
 */

double reduced_BSP::get_overwrite_prob(Recombination &r, double lb, double ub) {
    if (check_points.count(r.pos) > 0) {
        return 0.0;
    }
    double join_time = r.inserted_node->time;
    double p1 = cc->weight(lb, ub);
    double p2 = cc->weight(max(cut_time, r.start_time), join_time);
    if (p1 == 0 and p2 == 0) {
        return 1.0;
    }
    double overwrite_prob = p2/(p1 + p2);
    assert(!isnan(overwrite_prob));
    return overwrite_prob;
}

void reduced_BSP::process_interval(Recombination &r, int i) {
    Branch prev_branch = curr_intervals[i]->branch;
    if (prev_branch == r.source_branch) {
        process_source_interval(r, i);
    } else if (prev_branch == r.target_branch) {
        process_target_interval(r, i);
    } else {
        process_other_interval(r, i);
    }
}

void reduced_BSP::process_source_interval(Recombination &r, int i) {
    double w1, w2, lb, ub = 0;
    Interval_ptr prev_interval = curr_intervals[i];
    double p = forward_probs[curr_index - 1][i];
    double point_time = r.source_branch.upper_node->time;
    double break_time = r.start_time;
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

void reduced_BSP::process_target_interval(Recombination &r, int i) {
    double w0, w1, w2, lb, ub = 0;
    Interval_ptr prev_interval = curr_intervals[i];
    double p = forward_probs[curr_index - 1][i];
    double join_time = r.inserted_node->time;
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

void reduced_BSP::process_other_interval(Recombination &r, int i) {
    double lb, ub = 0;
    Interval_ptr prev_interval = curr_intervals[i];
    double p = forward_probs[curr_index - 1][i];
    if (!r.affect(prev_interval->branch)) {
        if (p > cutoff or prev_interval->full(cut_time)) {
            temp_intervals.push_back(prev_interval);
            temp.push_back(p);
        }
    } else {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        Branch &next_branch = r.merging_branch;
        Interval_info next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    }
}

double reduced_BSP::random() {
    double p = uniform_random();
    return p;
}

int reduced_BSP::get_prev_breakpoint(int x) {
    auto state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->first;
}

vector<Interval_ptr> &reduced_BSP::get_state_space(int x) {
    auto state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->second;
}

vector<double> &reduced_BSP::get_time_points(int x) {
    auto time_it = times.upper_bound(x);
    time_it--;
    return time_it->second;
}

vector<double> &reduced_BSP::get_raw_weights(int x) {
    auto weight_it = weights.upper_bound(x);
    weight_it--;
    return weight_it->second;
}

int reduced_BSP::get_interval_index(Interval_ptr interval, vector<Interval_ptr > &intervals) {
    auto it = find(intervals.begin(), intervals.end(), interval);
    int index = (int) distance(intervals.begin(), it);
    return index;
}
 
void reduced_BSP::simplify(map<double, Branch> &joining_branches) {
    map<double, Branch> simplified_joining_branches = {};
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

Interval_ptr reduced_BSP::sample_curr_interval(int x) {
    vector<Interval_ptr > &intervals = get_state_space(x);
    double ws = accumulate(forward_probs[x].begin(), forward_probs[x].end(), 0.0);
    double q = random();
    double w = ws*q;
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

Interval_ptr reduced_BSP::sample_prev_interval(int x) {
    vector<Interval_ptr > &intervals = get_state_space(x);
    vector<double> &prev_times = get_time_points(x);
    double rho = rhos[x];
    double ws = recomb_sums[x];
    double q = random();
    double w = ws*q;
    for (int i = 0; i < intervals.size(); i++) {
        w -= get_recomb_prob(rho, prev_times[i])*forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            return intervals[i];
        }
    }
    cerr << "bsp prev sampling failed" << endl;
    exit(1);
}

Interval_ptr reduced_BSP::sample_source_interval(Interval_ptr interval, int x) {
    vector<Interval_ptr> &intervals = interval->intervals;
    vector<double> &weights = interval->source_weights;
    vector<Interval_ptr> &prev_intervals = get_state_space(x);
    if (x == interval->start_pos - 1) {
        double q = random();
        double ws = accumulate(weights.begin(), weights.end(), 0.0);
        double w = ws*q;
        for (int i = 0; i < weights.size(); i++) {
            w -= weights[i];
            if (w <= 0) {
                sample_index = get_interval_index(intervals[i], prev_intervals);
                return intervals[i];
            }
        }
        cerr << "sampling failed" << endl;
        exit(1);
    } else {
        sample_index = get_interval_index(interval, prev_intervals);
        assert(prev_intervals[sample_index] == interval);
        return interval;
    }
}

int reduced_BSP::trace_back_helper(Interval_ptr interval, int x) {
    int y = get_prev_breakpoint(x);
    if (!interval->full(cut_time)) {
        return y;
    }
    vector<double> &ts = get_time_points(x);
    vector<double> &ws = get_raw_weights(x);
    assert(sample_index < ts.size());
    double t = ts[sample_index];
    double w = ws[sample_index];
    double p = random();
    double q = 1;
    double shrinkage = 0;
    double recomb_prob = 0;
    double non_recomb_prob = 0;
    double all_prob = 0;
    while (x > y) {
        recomb_sum = recomb_sums[x - 1];
        weight_sum = weight_sums[x];
        if (recomb_sum == 0) {
            shrinkage = 1;
        } else {
            recomb_prob = get_recomb_prob(rhos[x - 1], t);
            non_recomb_prob = (1 - recomb_prob)*forward_probs[x - 1][sample_index];
            all_prob = non_recomb_prob + recomb_sum*w*recomb_prob/weight_sum;
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

/*
int reduced_BSP::trace_back_helper(Interval_ptr interval, int x) {
    int y = get_prev_breakpoint(x);
    if (!interval->full(cut_time)) {
        return y;
    }
    vector<double> &ts = get_time_points(x);
    vector<double> &ws = get_raw_weights(x);
    assert(sample_index < ts.size());
    double t = ts[sample_index];
    double w = ws[sample_index];
    double p = random();
    double q = 1;
    double shrinkage = 0;
    double recomb_prob = 0;
    double non_recomb_prob = 0;
    double all_prob = 0;
    while (x > y) {
        recomb_sum = recomb_sums[x - 1];
        weight_sum = weight_sums[x];
        reduced_sum = reduced_sums[x];
        if (recomb_sum == 0) {
            shrinkage = 1;
        } else {
            recomb_prob = get_recomb_prob(rhos[x - 1], t);
            non_recomb_prob = (1 - recomb_prob)*forward_probs[x - 1][sample_index];
            all_prob = non_recomb_prob + recomb_sum*w*recomb_prob/weight_sum + recomb_sum*forward_probs[x - 1][sample_index]*(1 - reduced_sum);
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
 */

double reduced_BSP::avg_num_states() {
    int span = 0;
    double count = 0;
    auto x = state_spaces.begin();
    ++x;
    while (x->first != INT_MAX) {
        count += x->second.size()*(x->first - prev(x)->first);
        span = x->first;
        ++x;
    }
    double avg = (double) count/span;
    return avg;
}

