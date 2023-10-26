//
//  Trace_pruner.cpp
//  SINGER
//
//  Created by Yun Deng on 4/25/23.
//

#include "Trace_pruner.hpp"

#include "Trace_pruner.hpp"

Trace_pruner::Trace_pruner() {}

void Trace_pruner::prune_arg(ARG &a) {
    cut_time = a.cut_time;
    queries = a.removed_branches;
    start = a.removed_branches.begin()->first;
    end = a.removed_branches.rbegin()->first;
    segments.insert({start, end});
    used_seeds = {start, end};
    seed_trees[start] = a.start_tree;
    seed_trees[start].internal_cut(cut_time);
    deletions[end] = {};
    insertions[end] = {};
    build_match_map(a);
    float x = 0;
    while (potential_seeds.size() > 2 or segments.size() > 0) {
        x = find_minimum_match();
        extend(a, x);
    }
    assert(segments.size() == 0);
}

void Trace_pruner::set_check_points(set<float> &p) {
    check_points = p;
}

void Trace_pruner::start_search(ARG &a, float m) {
    seed_scores.clear();
    Node_ptr n = get_node_at(m);
    float mismatch = 0;
    float lb, ub;
    Interval_info interval;
    float x0 = find_closest_reference(m);
    seed_trees[m] = a.internal_modify_tree_to(m, seed_trees[x0], x0);
    length += abs(m - x0);
    float min_mismatch = INT_MAX;
    for (auto &x : seed_trees[m].parents) {
        if (x.second->time > cut_time) {
            mismatch = count_mismatch(Branch(x.first, x.second), n, m);
            min_mismatch = min(mismatch, min_mismatch);
        }
    }
    for (auto &x : seed_trees[m].parents) {
        if (x.second->time > cut_time) {
            Branch b = Branch(x.first, x.second);
            mismatch = count_mismatch(b, n, m);
            if (mismatch == min_mismatch) {
                lb = max(cut_time, b.lower_node->time);
                ub = b.upper_node->time;
                interval = Interval_info(b, lb, ub);
                interval.seed_pos = m;
                seed_scores[interval] = 1;
            }
        }
    }
    restrict_search();
    potential_seeds.erase(m);
    assert(seed_scores.size() > 0 and seed_scores.size() <= band_width);
}

void Trace_pruner::write_reduction_distance(ARG &a, string filename) {
    float start = a.removed_branches.begin()->first;
    float end = a.removed_branches.rbegin()->first;
    auto join_it = a.joining_branches.begin();
    auto reduced_it = reductions.begin();
    auto recomb_it = a.recombinations.lower_bound(start);
    Branch joining_branch = Branch();
    Tree tree = Tree();
    float x = start;
    set<Branch> reduced_set = reduced_it->second;
    ofstream file;
    file.open(filename);
    while (join_it->first < end) {
        x = join_it->first;
        joining_branch = join_it->second;
        while (reduced_it->first <= x) {
            reduced_set = reduced_it->second;
            reduced_it++;
        }
        while (recomb_it->first <= x) {
            Recombination &r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
        int min_distance = INT_MAX;
        int distance = INT_MAX;
        for (Branch b : reduced_set) {
            distance = tree.distance(joining_branch.lower_node, b.lower_node);
            min_distance = min(min_distance, distance);
        }
        file << join_it->first << " " << min_distance << "\n";
        join_it++;
    }
}

void Trace_pruner::write_reduction_size(string filename) {
    ofstream file;
    file.open(filename);
    for (auto x : reductions) {
        file << x.first << " " << x.second.size() << "\n";
    }
}

float Trace_pruner::min_reduction_error() {
    float error = 0;
    set<Branch> &reduced_set = reductions.begin()->second;
    Node_ptr node;
    auto reduced_it = reductions.begin();
    for (auto x : match_map) {
        float m = x.first;
        while (reduced_it->first <= m) {
            reduced_it++;
            reduced_set = reduced_it->second;
        }
        node = get_node_at(m);
        float min_mismatch = 1;
        for (const Branch &b : reduced_set) {
            float mismatch = count_mismatch(b, node, m);
            min_mismatch = min(min_mismatch, mismatch);
        }
        if (min_mismatch > 0) {
            error += 1;
        }
    }
    return error;
}

void Trace_pruner::write_reductions(ARG &a) {
    assert(segments.size() == 0);
    auto delete_it = deletions.begin();
    auto insert_it = insertions.begin();
    int start_index = a.get_index(start);
    int end_index = a.get_index(end);
    set<Interval_info> reduced_set = {};
    bool state_change = false;
    for (int i = start_index; i < end_index; i++) {
        while (delete_it->first <= a.coordinates[i]) {
            for (auto x : delete_it->second) {
                assert(reduced_set.count(x) > 0);
                reduced_set.erase(x);
            }
            for (auto x : insert_it->second) {
                reduced_set.insert(x);
            }
            delete_it++;
            insert_it++;
            state_change = true;
        }
        if (state_change) {
            for (auto x : reduced_set) {
                reductions[a.coordinates[i]].insert(x.branch);
            }
            state_change = false;
        }
    }
    reductions[a.sequence_length] = {};
}

Node_ptr Trace_pruner::get_node_at(float x) {
    auto query_it = queries.upper_bound(x);
    query_it--;
    return query_it->second.lower_node;
}

float Trace_pruner::count_mismatch(Branch branch, Node_ptr n, float m) {
    float s0 = n->get_state(m);
    float sl = branch.lower_node->get_state(m);
    float su = branch.upper_node->get_state(m);
    if (branch.upper_node->index != -1) {
        if (abs(sl - s0) > 0.5 and abs(su - s0) > 0.5) {
            return 1;
        } else {
            return 0;
        }
    } else { // root branch mutation case
        if (sl == 1 or sl == 0.5) { // you can join the root as 1-0
            return 0;
        } else if (s0 == sl) { // else must be a zero
            return 0;
        } else {
            return 1;
        }
    }
}

/*
float Trace_pruner::count_mismatch(Branch branch, Node_ptr n, float m) {
    if (private_mutations.count(m) > 0) {
        return 0.3;
    }
    float s0 = n->get_state(m);
    float sl = branch.lower_node->get_state(m);
    float su = branch.upper_node->get_state(m);
    if (branch.upper_node->index != -1) {
        if (abs(sl - s0) > 0.5 and abs(su - s0) > 0.5) {
            return 1;
        } else {
            return 0;
        }
    } else { // root branch mutation case
        if (sl == 1 or sl == 0.5) { // you can join the root as 1-0
            return 0;
        } else if (s0 == sl) { // else must be a zero
            return 0;
        } else {
            return 1;
        }
    }
}
 */

void Trace_pruner::build_match_map(ARG &a) {
    float m = 0;
    float inf = INT_MAX;
    float lb = 0;
    auto mb_it = a.mutation_branches.lower_bound(start);
    Node_ptr n = nullptr;
    while (mb_it->first < end) {
        m = mb_it->first;
        n = get_node_at(m);
        lb = get_match_time(mb_it->second, m, n);
        if (lb >= 0) { // when negative, this site is uninformative and not suitable for checking either
            match_map[m] = lb;
        }
        mb_it++;
    }
    match_map[end] = inf;
    match_map[start] = inf;
    potential_seeds = match_map;
}

/*
float Trace_pruner::find_closest_reference(float x) {
    auto tree_it = seed_trees.upper_bound(x);
    tree_it--;
    return tree_it->first;
}
 */

float Trace_pruner::find_closest_reference(float x) {
    auto r_it = seed_trees.lower_bound(x);
    auto l_it = seed_trees.upper_bound(x);
    l_it--;
    float ld = abs(l_it->first - x);
    float rd = abs(r_it->first - x);
    if (ld <= rd) {
        return l_it->first;
    } else {
        return r_it->first;
    }
}

float Trace_pruner::find_minimum_match() {
    float x = 0;
    if (potential_seeds.size() == 2) {
        auto it = segments.begin();
        x = 0.5*(it->first + it->second);
        if (floor(x) == x) {
            x += 0.5; // x better not be a integer, just for programming convenience
        }
        match_map[x] = 0.0;
    } else {
        auto it = min_element(potential_seeds.begin(), potential_seeds.end(),
                              [](const auto& l, const auto& r) { return l.second < r.second;});
        x = it->first;
    }
    return x;
}

void Trace_pruner::extend_forward(ARG &a, float x) {
    int index = a.get_index(x);
    float m = x;
    auto recomb_it = a.recombinations.upper_bound(x);
    auto match_it = match_map.lower_bound(x);
    auto used_it = used_seeds.upper_bound(x);
    match_it++;
    float ub = *used_it;
    Node_ptr n;
    float bin_start = 0;
    float bin_end = 0;
    set<float> mutations = {};
    while (a.coordinates[index] < ub and curr_scores.size() > 0) {
        for (float y : mutations) {
            potential_seeds.erase(y);
        }
        bin_start = a.coordinates[index];
        bin_end = a.coordinates[index + 1];
        while (recomb_it->first <= bin_start) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            recombination_forward(r);
        }
        mutations.clear();
        while (match_it->first < bin_end) {
            m = match_it->first;
            n = get_node_at(m);
            mutation_update(n, m);
            mutations.insert(m);
            ++match_it;
        }
        if (a.recombinations.count(bin_end) == 0) {
            forward_prune_states(bin_end);
        }
        ++index;
    }
    if (curr_scores.size() > 0) {
        if (ub == a.sequence_length) {
            bin_end = ub;
        } else {
            index = a.get_index(ub);
            bin_end = a.coordinates[index + 1];
        }
        delete_all(bin_end);
    }
    remove_segment(x, bin_end);
    // cout << "Extend forward from " << x << " to " << m << endl;
}

void Trace_pruner::extend_backward(ARG &a, float x) {
    int index = a.get_index(x);
    float m = x;
    auto recomb_it = a.recombinations.upper_bound(x);
    auto match_it = match_map.lower_bound(x);
    auto used_it = used_seeds.upper_bound(x);
    --recomb_it;
    --used_it;
    float lb = *used_it;
    Node_ptr n;
    float bin_start = 0;
    float bin_end = 0;
    set<float> mutations = {};
    while (a.coordinates[index + 1] > lb and curr_scores.size() > 0) {
        for (float y : mutations) {
            potential_seeds.erase(y);
        }
        bin_start = a.coordinates[index];
        bin_end = a.coordinates[index + 1];
        while (recomb_it->first > bin_start) {
            Recombination &r = recomb_it->second;
            recombination_backward(r);
            --recomb_it;
        }
        mutations.clear();
        while (match_it->first >= bin_start) {
            m = match_it->first;
            n = get_node_at(m);
            mutation_update(n, m);
            mutations.insert(m);
            if (match_it == match_map.begin()) {
                break;
            }
            --match_it;
        }
        if (a.recombinations.count(bin_start) == 0) {
            backward_prune_states(bin_start);
        }
        --index;
    }
    if (curr_scores.size() > 0) {
        index = a.get_index(lb);
        bin_start = a.coordinates[index];
        insert_all(bin_start);
    }
    remove_segment(bin_start, x);
    // cout << "Extend backward from " << x << " to " << m << endl;
}

void Trace_pruner::extend(ARG &a, float x) {
    start_search(a, x);
    curr_scores = seed_scores;
    extend_forward(a, x);
    curr_scores = seed_scores;
    extend_backward(a, x);
    used_seeds.insert(x); // when extending later seeds, don't go beyond previous seeds (to save computation)
}

void Trace_pruner::mutation_update(Node_ptr n, float m) {
    if (n == nullptr) {
        return;
    }
    float mismatch = 0;
    float penalty = 0;
    for (auto &[i, s] : curr_scores) {
        mismatch = count_mismatch(i.branch, n, m);
        penalty = pow(mut_prob, mismatch);
        s *= penalty;
    }
}

void Trace_pruner::recombination_forward(Recombination &r) {
    transition_scores.clear();
    for (auto &[i, s] : curr_scores) {
        forward_transition(r, i);
        assert(!r.affect(i.branch) or deletions[r.pos].count(i) > 0);
    }
    curr_scores = transition_scores;
}

void Trace_pruner::recombination_backward(Recombination &r) {
    transition_scores.clear();
    for (auto &[i, s] : curr_scores) {
        backward_transition(r, i);
        assert(!r.affect(i.branch) or insertions[r.pos].count(i) > 0);
    }
    curr_scores = transition_scores;
}

void Trace_pruner::forward_transition(Recombination &r, const Interval_info &interval) {
    float lb, ub;
    float w0, w1, w2;
    Interval_info new_interval;
    Branch b;
    float p = curr_scores[interval];
    p = max(p, cutoff*0.1f);
    float l = interval.lb;
    float u = interval.ub;
    if (!r.affect(interval.branch)) {
        transition_scores[interval] = curr_scores[interval];
    } else if (interval.branch == r.source_branch) {
        if (l <= r.start_time) {
            lb = min(l, r.start_time);
            ub = min(u, r.start_time);
            w1 = exp_prop(interval.lb, interval.ub, lb, ub);
            new_interval = Interval_info(r.recombined_branch, lb, ub);
            forward_transition_helper(interval, new_interval, r.pos, w1*p);
        }
        if (u > r.start_time) {
            lb = max(r.start_time, l);
            ub = max(r.start_time, u);
            w2 = exp_prop(l, u, lb, ub);
            new_interval = Interval_info(r.merging_branch, r.deleted_node->time, r.deleted_node->time);
            forward_transition_helper(interval, new_interval, r.pos, w2*p);
        }
    } else if (interval.branch == r.target_branch) {
        w0 = forward_overwrite_prob(r, l, u);
        if (l < r.inserted_node->time and w0 < 1) {
            lb = min(l, r.inserted_node->time);
            ub = min(u, r.inserted_node->time);
            w1 = exp_prop(l, u, lb, ub);
            new_interval = Interval_info(r.lower_transfer_branch, lb, ub);
            forward_transition_helper(interval, new_interval, r.pos, (1 - w0)*w1*p);
        }
        if (u > r.inserted_node->time and w0 < 1) {
            lb = max(r.inserted_node->time, l);
            ub = max(r.inserted_node->time, u);
            w2 = exp_prop(l, u, lb, ub);
            new_interval = Interval_info(r.upper_transfer_branch, lb, ub);
            forward_transition_helper(interval, new_interval, r.pos, (1 - w0)*w2*p);
        }
        if (w0*p > 0) {
            lb = max(cut_time, r.start_time);
            ub = r.inserted_node->time;
            new_interval = Interval_info(r.recombined_branch, lb, ub);
            forward_transition_helper(interval, new_interval, r.pos, w0*p);
        }
    } else {
        new_interval = Interval_info(r.merging_branch, l, u);
        forward_transition_helper(interval, new_interval, r.pos, p);
    }
    assert(!r.affect(interval.branch) or deletions[r.pos].count(interval) > 0);
}

void Trace_pruner::backward_transition(Recombination &r, const Interval_info &interval) {
    float lb, ub;
    float w0, w1, w2;
    Interval_info new_interval;
    Branch b;
    float p = curr_scores[interval];
    p = max(p, cutoff*0.1f);
    float l = interval.lb;
    float u = interval.ub;
    float x = r.pos;
    if (!r.create(interval.branch)) {
        transition_scores[interval] = curr_scores[interval];
    } else if (interval.branch == r.recombined_branch) {
        if (l <= r.start_time) {
            lb = min(l, r.start_time);
            ub = min(u, r.start_time);
            w1 = exp_prop(interval.lb, interval.ub, lb, ub);
            new_interval = Interval_info(r.source_branch, lb, ub);
            backward_transition_helper(interval, new_interval, x, w1*p);
        }
        if (u > r.start_time) {
            lb = max(r.start_time, l);
            ub = max(r.start_time, u);
            if (check_points.count(r.pos) > 0) { // forbid transition to point mass when check points
                w2 = 0.01*cutoff;
            } else {
                w2 = exp_prop(l, u, lb, ub);
            }
            new_interval = Interval_info(r.target_branch, r.inserted_node->time, r.inserted_node->time);
            backward_transition_helper(interval, new_interval, x, w2*p);
        }
    } else if (interval.branch == r.merging_branch) {
        w0 = backward_overwrite_prob(r, l, u);
        if (l < r.deleted_node->time and w0 < 1) {
            lb = min(l, r.deleted_node->time);
            ub = min(u, r.deleted_node->time);
            w1 = exp_prop(l, u, lb, ub);
            new_interval = Interval_info(r.source_sister_branch, lb, ub);
            backward_transition_helper(interval, new_interval, r.pos, (1 - w0)*w1*p);
        }
        if (u > r.deleted_node->time and w0 < 1) {
            lb = max(r.deleted_node->time, l);
            ub = max(r.deleted_node->time, u);
            w2 = exp_prop(l, u, lb, ub);
            new_interval = Interval_info(r.source_parent_branch, lb, ub);
            backward_transition_helper(interval, new_interval, r.pos, (1 - w0)*w2*p);
        }
        if (w0*p > 0) {
            lb = max(cut_time, r.start_time);
            ub = r.deleted_node->time;
            new_interval = Interval_info(r.source_branch, lb, ub);
            backward_transition_helper(interval, new_interval, r.pos, w0*p);
        }
    } else {
        new_interval = Interval_info(r.target_branch, l, u);
        backward_transition_helper(interval, new_interval, r.pos, p);
    }
    assert(!r.affect(interval.branch) or insertions[r.pos].count(interval) > 0);
}

void Trace_pruner::forward_transition_helper(Interval_info prev_interval, Interval_info next_interval, float x, float p) {
    if (p == 0) {
        return;
    }
    assert(prev_interval.lb >= cut_time and next_interval.lb >= cut_time);
    next_interval.seed_pos = prev_interval.seed_pos;
    transition_scores[next_interval] += p;
    deletions[x].insert(prev_interval);
    insertions[x].insert(next_interval);
}

void Trace_pruner::backward_transition_helper(Interval_info next_interval, Interval_info prev_interval, float x, float p) {
    if (p == 0) {
        return;
    }
    assert(prev_interval.lb >= cut_time and next_interval.lb >= cut_time);
    prev_interval.seed_pos = next_interval.seed_pos;
    transition_scores[prev_interval] += p;
    deletions[x].insert(prev_interval);
    insertions[x].insert(next_interval);
}

void Trace_pruner::forward_prune_states(float x) {
    auto it = curr_scores.begin();
    while (it != curr_scores.end()) {
        if (it->second < cutoff) {
            deletions[x].insert(it->first);
            it = curr_scores.erase(it);
            insertions[x];
        } else {
            ++it;
        }
    }
}

void Trace_pruner::backward_prune_states(float x) {
    auto it = curr_scores.begin();
    while (it != curr_scores.end()) {
        if (it->second < cutoff) {
            insertions[x].insert(it->first);
            it = curr_scores.erase(it);
            deletions[x];
        } else {
            ++it;
        }
    }
}

void Trace_pruner::delete_all(float x) {
    for (auto &[i, s] : curr_scores) {
        deletions[x].insert(i);
        insertions[x];
    }
}

void Trace_pruner::insert_all(float x) {
    for (auto &[i, s] : curr_scores) {
        insertions[x].insert(i);
        deletions[x];
    }
}

float Trace_pruner::exp_prob(float l, float u) {
    float prob = exp(-l) - exp(-u);
    return prob;
}

float Trace_pruner::exp_prop(float l, float u, float x, float y) {
    if (l == u) {
        if (x == y and x == l) {
            return 1;
        } else {
            return 0;
        }
    }
    float sub_prob = exp_prob(x, y);
    float prob = exp_prob(l, u);
    float prop;
    if (sub_prob == 0) {
        prop = (y - x)/(u - l);
    } else {
        prop = sub_prob/prob;
    }
    assert(!isnan(prop));
    return prop;
}

float Trace_pruner::exp_median(float l, float u) {
    float lcdf = 1 - exp(-l);
    float ucdf = 1 - exp(-u);
    float cdf = 0.5*(lcdf + ucdf);
    float m = -log(1 - cdf);
    return m;
}

float Trace_pruner::forward_overwrite_prob(Recombination &r, float lb, float ub) {
    if (lb > r.inserted_node->time or ub < r.inserted_node->time) { // this dominates
        return 0;
    }
    if (check_points.count(r.pos) > 0) { // reject illegal transitions
        return 0.01*cutoff;
    }
    if (lb == ub and lb == r.inserted_node->time) {
        return 1;
    }
    float p1 = exp_prob(lb, ub);
    float p2 = exp_prob(r.start_time, r.inserted_node->time);
    float p = p2/(p1 + p2);
    assert(!isnan(p) and p <= 1);
    return p;
}

float Trace_pruner::backward_overwrite_prob(Recombination &r, float lb, float ub) {
    if (lb > r.deleted_node->time or ub < r.deleted_node->time) {
        return 0;
    }
    if (lb == ub and lb == r.inserted_node->time) {
        return 1;
    }
    float p1 = exp_prob(lb, ub);
    float p2 = exp_prob(r.start_time, r.deleted_node->time);
    float p = p2/(p1 + p2);
    if (p1 + p2 == 0) {
        p = 0.01*cutoff;
    }
    assert(!isnan(p));
    return p;
}

void Trace_pruner::remove_segment(float x, float y) {
    assert(x < y);
    auto it = segments.begin();
    set<pair<float, float>> overlaps = {};
    set<pair<float, float>> new_segments = {};
    while (it != segments.end() and it->first < y) {
        float start = it->first;
        float end = it->second;
        if (start >= x and end <= y) {
            overlaps.insert(*it);
        }
        if (start <= x and end > x) {
            new_segments.insert({start, x});
            overlaps.insert(*it);
        }
        if (end >= y) {
            new_segments.insert({y, end});
            overlaps.insert(*it);
        }
        ++it;
    }
    for (auto &x : overlaps) {
        segments.erase(x);
    }
    for (auto &x : new_segments) {
        if (x.second > x.first) {
            segments.insert(x);
        }
    }
}

float Trace_pruner::get_match_time(set<Branch> &branches, float m, Node_ptr n) {
    float state = n->get_state(m);
    assert(state == 0 or state == 1);
    int valid_count = 0;
    for (const Branch &b : branches) {
        if (b.upper_node->time > cut_time) {
            valid_count += 1; // needs at least one mutation above cut_time (excluding root mutations)
        }
    }
    if (valid_count == 0) {
        return -1; // essentially a private mutation, not that informative
    }
    if (valid_count > 1) {
        return max_time; // multiple mutations but at least one above cut_time, still for checking
    }
    assert(valid_count == 1);
    const Branch &b = *branches.rbegin(); // the single valid mutation branch, which must be the highest branch
    if (b.lower_node->get_state(m) == state) { // a matching branch
        return max(b.lower_node->time, cut_time) - cut_time;
    } else if (b.upper_node->index == -1) {
        return 0.0;
    } else { // non-matching branch
        return max_time - max(cut_time, b.lower_node->time) - cut_time;
    }
}

/*
float Trace_pruner::get_match_time(set<Branch> &branches, float m, Node_ptr n) {
    float state = n->get_state(m);
    assert(state == 0 or state == 1);
    int valid_count = 0;
    for (const Branch &b : branches) {
        if (b.upper_node->time > cut_time) {
            valid_count += 1; // needs at least one mutation above cut_time (excluding root mutations)
        }
    }
    if (valid_count == 0) {
        if (state == 1) {
            private_mutations.insert(m);
            return max_time; // essentially a private mutation, not that informative
        } else {
            return -1; // this mutation is not in the picture at all
        }
    }
    if (valid_count > 1) {
        return max_time; // multiple mutations but at least one above cut_time, still for checking
    }
    assert(valid_count == 1);
    const Branch &b = *branches.rbegin(); // the single valid mutation branch, which must be the highest branch
    if (b.lower_node->get_state(m) == state) { // a matching branch
        return max(b.lower_node->time, cut_time) - cut_time;
    } else if (b.upper_node->index == -1) {
        return 0.0;
    } else { // non-matching branch
        return max_time - max(cut_time, b.lower_node->time) - cut_time;
    }
}
 */

void Trace_pruner::restrict_search() {
    if (seed_scores.size() <= band_width) {
        return; // If size is less than or equal to band_with, do nothing
    }
    vector<Interval_info> seeds;
    transform(seed_scores.begin(), seed_scores.end(), back_inserter(seeds),
                       [](const auto& pair) { return pair.first; });
    // Shuffle using existing random_engine
    shuffle(seeds.begin(), seeds.end(), random_engine);

    for (int i = band_width; i < seeds.size(); ++i) {
        seed_scores.erase(seeds[i]);
    }
    assert(seed_scores.size() <= band_width);
}
