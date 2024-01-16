//
//  BSP_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 4/2/23.
//

#include "BSP_smc.hpp"

BSP_smc::BSP_smc() {}

BSP_smc::~BSP_smc() {
    for (auto &x : state_spaces) {
        for (Interval *interval : x.second) {
            delete interval;
        }
    }
    vector<vector<double>>().swap(forward_probs);
    // map<Interval *, vector<double>>().swap(source_weights);
    // map<Interval *, vector<Interval *>>().swap(source_intervals);
    map<int, vector<Interval *>>().swap(state_spaces);
}

void BSP_smc::reserve_memory(int length) {
    forward_probs.reserve(length);
}

void BSP_smc::start(set<Branch> &branches, double t) {
    cut_time = t;
    curr_index = 0;
    for (Branch b : branches) {
        if (b.upper_node->time > cut_time) {
            valid_branches.insert(b);
        }
    }
    double lb = 0;
    double ub = 0;
    double p = 0;
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
            temp.push_back(p);
        }
    }
    forward_probs.push_back(temp);
    compute_interval_info();
    weight_sums.push_back(0.0);
    set_dimensions();
    state_spaces[curr_index] = curr_intervals;
    temp.clear();
}

void BSP_smc::set_cutoff(double x) {
    cutoff = x;
}

void BSP_smc::set_emission(shared_ptr<Emission> e) {
    eh = e;
}

void BSP_smc::set_check_points(set<double> &p) {
    check_points = p;
}

void BSP_smc::forward(double rho) {
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


void BSP_smc::transfer(Recombination &r) {
    rhos.push_back(0);
    prev_rho = -1;
    prev_theta = -1;
    recomb_sums.push_back(0);
    weight_sums.push_back(0);
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
    state_spaces[curr_index] = curr_intervals;
}

double BSP_smc::get_recomb_prob(double rho, double t) {
    double p = rho*(t - cut_time)*exp(-rho*(t - cut_time));
    return p;
}

void BSP_smc::null_emit(double theta, Node_ptr query_node) {
    compute_null_emit_prob(theta, query_node);
    prev_theta = theta;
    prev_node = query_node;
    double ws = 0;
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] *= null_emit_probs[i];
        ws += forward_probs[curr_index][i];
    }
    assert(ws > 0);
    for (int i = 0; i < dim; i++) {
        forward_probs[curr_index][i] /= ws;
    }
}

void BSP_smc::mut_emit(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node) {
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

map<double, Branch> BSP_smc::sample_joining_branches(int start_index, vector<double> &coordinates) {
    prev_rho = -1;
    map<double, Branch> joining_branches = {};
    int x = curr_index;
    double pos = coordinates[x + start_index + 1];
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

void BSP_smc::check_recomb_sums() {
    double ws, rs= 0;
    double rho;
    for (int i = 0; i < rhos.size() - 1; i++) {
        vector<Interval *> intervals = get_state_space(i);
        vector<double> r = vector<double>(intervals.size());
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

void BSP_smc::write_forward_probs(string filename) {
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

double BSP_smc::avg_num_states() {
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

// private methods:

void BSP_smc::update_states(set<Branch> &deletions, set<Branch> &insertions) {
    for (Branch b : deletions) {
        if (b.upper_node->time > cut_time) {
            assert(valid_branches.count(b) > 0);
            valid_branches.erase(b);
        }
        states_change = true;
    }
    for (Branch b : insertions) {
        if (b.upper_node->time > cut_time) {
            valid_branches.insert(b);
        }
        states_change = true;
    }
    // states_change = true;
}

void BSP_smc::set_dimensions() {
    dim = (int) curr_intervals.size();
    recomb_probs.resize(dim); recomb_probs.assign(dim, 0);
    recomb_weights.resize(dim); recomb_weights.assign(dim, 0);
    null_emit_probs.resize(dim); null_emit_probs.assign(dim, 0);
    mut_emit_probs.resize(dim); mut_emit_probs.assign(dim, 0);
}

void BSP_smc::compute_recomb_probs(double rho) {
    if (prev_rho == rho) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        recomb_probs[i] = get_recomb_prob(rho, curr_intervals[i]->time);
    }
}

void BSP_smc::compute_recomb_weights(double rho) {
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

void BSP_smc::compute_null_emit_prob(double theta, Node_ptr query_node) {
    if (theta == prev_theta and query_node == prev_node) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        null_emit_probs[i] = eh->null_emit(curr_intervals[i]->branch, curr_intervals[i]->time, theta, query_node);
    }
}

void BSP_smc::compute_mut_emit_probs(double theta, double bin_size, set<double> &mut_set, Node_ptr query_node) {
    for (int i = 0; i < dim; i++) {
        mut_emit_probs[i] = eh->mut_emit(curr_intervals[i]->branch, curr_intervals[i]->time, theta, bin_size, mut_set, query_node);
    }
}

void BSP_smc::transfer_helper(Interval_info next_interval, Interval *prev_interval, double w) {
    transfer_weights[next_interval].push_back(w);
    transfer_intervals[next_interval].push_back(prev_interval);
}

void BSP_smc::transfer_helper(Interval_info next_interval) {
    transfer_weights[next_interval];
    transfer_intervals[next_interval];
}

void BSP_smc::add_new_branches(Recombination &r) { // add recombined branch and merging branch, if legal
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

void BSP_smc::compute_interval_info() {
    if (states_change) {
        cc->compute(valid_branches);
    }
    states_change = false;
    double t;
    double p;
    for (Interval *i : curr_intervals) {
        p = cc->weight(i->lb, i->ub);
        t = cc->time(i->lb, i->ub);
        i->assign_weight(p);
        i->assign_time(t);
    }
}

void BSP_smc::sanity_check(Recombination &r) {
    for (int i = 0; i < curr_intervals.size(); i++) {
        Interval *interval = curr_intervals[i];
        if (interval->lb == interval->ub and interval->lb == r.inserted_node->time and interval->branch != r.target_branch) {
            forward_probs[curr_index][i] = 0;
        }
    }
}

void BSP_smc::generate_intervals(Recombination &r) {
    Branch b;
    double lb;
    double ub;
    double p;
    vector<Interval *> intervals;
    vector<double> weights;
    Interval_info interval;
    Interval *new_interval = nullptr;
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
            new_interval = new Interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                new_interval->source_weights = move(weights);
                new_interval->source_intervals = move(intervals);
            }
        } else if (p >= cutoff) { // partial intervals
            new_interval = new Interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                new_interval->source_weights = move(weights);
                new_interval->source_intervals = move(intervals);
            }
        }
    }
    forward_probs.push_back(temp);
    curr_intervals = temp_intervals;
    compute_interval_info();
}

/*
void BSP_smc::fast_generate_intervals(Recombination &r) {
    Branch b;
    double lb;
    double ub;
    double p;
    set<Branch> curr_branches = {};
    vector<Interval *> intervals;
    vector<double> weights;
    Interval_info interval;
    Interval *new_interval = nullptr;
    for (auto x : transfer_weights) {
        interval = x.first;
        weights = x.second;
        intervals = transfer_intervals[x.first];
        b = interval.branch;
        lb = interval.lb;
        ub = interval.ub;
        p = accumulate(weights.begin(), weights.end(), 0.0);
        assert(!isnan(p));
        if (valid_branches.count(b) == 0) {
            continue;
        }
        if (lb == max(cut_time, b.lower_node->time) and ub == b.upper_node->time) { // full intervals
            curr_branches.insert(b);
            new_interval = new Interval(b, lb, ub, curr_index);
            curr_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                source_weights[new_interval] = weights;
                source_intervals[new_interval] = intervals;
            }
        } else if (p >= cutoff) { // partial intervals
            new_interval = new Interval(b, lb, ub, curr_index);
            curr_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                source_weights[new_interval] = weights;
                source_intervals[new_interval] = intervals;
            }
        }
    }
    for (Branch b : valid_branches) {
        if (curr_branches.count(b) == 0) {
            lb = max(cut_time, b.lower_node->time);
            ub = b.upper_node->time;
            new_interval = new Interval(b, lb, ub, curr_index);
            curr_intervals.push_back(new_interval);
            temp.push_back(0);
        }
    }
    forward_probs.push_back(temp);
    compute_interval_info();
    temp.clear();
}
 */

double BSP_smc::get_overwrite_prob(Recombination &r, double lb, double ub) {
    if (check_points.count(r.pos) > 0) {
        // cout << "check point" << endl;
        return 0.0;
    }
    double join_time = r.inserted_node->time;
    double p1 = cc->weight(lb, ub);
    double p2 = cc->weight(max(cut_time, r.start_time), join_time);
    if (p1 == 0 and p2 == 0) {
        return 1.0;
    }
    double overwrite_prob = p2/(p1 + p2);
    // assert(!isnan(overwrite_prob));
    return overwrite_prob;
}

void BSP_smc::process_interval(Recombination &r, int i) {
    Branch prev_branch = curr_intervals[i]->branch;
    if (prev_branch == r.source_branch) {
        process_source_interval(r, i);
    } else if (prev_branch == r.target_branch) {
        process_target_interval(r, i);
    } else {
        process_other_interval(r, i);
    }
}

void BSP_smc::process_source_interval(Recombination &r, int i) {
    double w1, w2, lb, ub = 0;
    Interval *prev_interval = curr_intervals[i];
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

void BSP_smc::process_target_interval(Recombination &r, int i) {
    double w0, w1, w2, lb, ub = 0;
    Interval *prev_interval = curr_intervals[i];
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

void BSP_smc::process_other_interval(Recombination &r, int i) {
    double lb, ub = 0;
    Interval *prev_interval = curr_intervals[i];
    double p = forward_probs[curr_index - 1][i];
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

double BSP_smc::random() {
    // return (double) rand()/RAND_MAX;
    double p = uniform_random();
    return p;
}

int BSP_smc::get_prev_breakpoint(int x) {
    auto state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->first;
}

vector<Interval *> &BSP_smc::get_state_space(int x) {
    auto state_it = state_spaces.upper_bound(x);
    state_it--;
    return state_it->second;
}

int BSP_smc::get_interval_index(Interval *interval, vector<Interval *> &intervals) {
    auto it = find(intervals.begin(), intervals.end(), interval);
    // assert(it != intervals.end());
    int index = (int) distance(intervals.begin(), it);
    // assert(intervals[index] == interval);
    return index;
}
 
void BSP_smc::simplify(map<double, Branch> &joining_branches) {
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

Interval *BSP_smc::sample_curr_interval(int x) {
    vector<Interval *> &intervals = get_state_space(x);
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

Interval *BSP_smc::sample_prev_interval(int x) {
    vector<Interval *> &intervals = get_state_space(x);
    double rho = rhos[x];
    double ws = recomb_sums[x];
    double q = random();
    double w = ws*q;
    // assert(intervals.size() == forward_probs[x].size());
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

Interval *BSP_smc::sample_source_interval(Interval *interval, int x) {
    vector<Interval *> &prev_intervals = get_state_space(x);
    // assert(source_intervals.count(interval) > 0);
    // vector<Interval *> intervals = source_intervals[interval];
    // vector<double> weights = source_weights[interval];
    vector<Interval *> &intervals = interval->source_intervals;
    vector<double> &weights = interval->source_weights;
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
}

int BSP_smc::trace_back_helper(Interval *interval, int x) {
    if (!interval->full(cut_time)) {
        return interval->start_pos;
    }
    // int y = get_prev_breakpoint(x);
    double p = random();
    double q = 1;
    double shrinkage = 0;
    double recomb_prob = 0;
    double non_recomb_prob = 0;
    double all_prob = 0;
    while (x > interval->start_pos) {
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
    return interval->start_pos;
}

void BSP_smc::check_intervals() {
    set<Branch> full_branches = {};
    for (Interval *i : curr_intervals) {
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
}
