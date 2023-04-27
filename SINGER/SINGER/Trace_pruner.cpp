//
//  Trace_pruner.cpp
//  SINGER
//
//  Created by Yun Deng on 4/25/23.
//

#include "Trace_pruner.hpp"

Trace_pruner::Trace_pruner() {}

void Trace_pruner::prune_arg(ARG &a) {}

void Trace_pruner::start_search(ARG &a, float m) {
    seed_scores.clear();
    Node *n = get_node_at(m);
    float mismatch = 0;
    float lb, ub;
    Interval_info interval;
    float x0 = find_closest_reference(m);
    seed_trees[m] = seed_trees[x0];
    a.get_tree_at(m, seed_trees[m], x0);
    for (Branch b : seed_trees[m].branches) {
        mismatch = count_mismatch(b, n, m);
        if (mismatch == 0) {
            lb = lb = max(cut_time, b.lower_node->time);
            ub = b.upper_node->time;
            interval = Interval_info(b, lb, ub);
            seed_scores[interval] = 1;
        }
    }
    potential_seeds.erase(m);
}

Node *Trace_pruner::get_node_at(float x) {
    auto query_it = queries.upper_bound(x);
    query_it--;
    return query_it->second.lower_node;
}

float Trace_pruner::count_mismatch(Branch branch, Node *n, float m) {
    float s0 = n->get_state(m);
    float sl = branch.lower_node->get_state(m);
    float su = branch.upper_node->get_state(m);
    if (abs(sl - s0) > 0.5 and abs(su - s0) > 0.5) {
        return 1;
    } else {
        return 0;
    }
}

void Trace_pruner::build_match_map(ARG &a) {
    float m = 0;
    float inf = INT_MAX;
    float lb = 0;
    float state = 0;
    auto mb_it = a.mutation_branches.lower_bound(start);
    Node *n = nullptr;
    while (mb_it->first < end) {
        m = mb_it->first;
        n = get_node_at(m);
        state = n->get_state(m);
        lb = -1;
        for (Branch b : mb_it->second) {
            if (b.upper_node->time < n->time) {
                continue; // don't search if the mutation is absolutely below the query node
            }
            if (b.lower_node->get_state(m) == state) {
                lb = max(lb, max(b.lower_node->time, n->time));
            } else if (b.upper_node->index == -1) {
                lb = max(lb, 0.0f);
            } else {
                lb = max(lb, max_time - b.lower_node->time);
            }
        }
        if (lb >= 0) { // when negative, no mutation is legal, thus no match was found
            match_map[m] = lb - n->time;
        }
        mb_it++;
    }
    match_map[end] = inf;
    match_map[start] = inf;
    potential_seeds = match_map;
}

float Trace_pruner::find_closest_reference(float x) {
    auto tree_it = seed_trees.upper_bound(x);
    tree_it--;
    return tree_it->first;
}

float Trace_pruner::find_minimum_match() {
    auto it = min_element(potential_seeds.begin(), potential_seeds.end(),
                          [](const auto& l, const auto& r) { return l.second < r.second;});
    float x = it->first;
    return x;
}

void Trace_pruner::extend_forward(ARG &a, float x) {
    auto recomb_it = a.recombinations.upper_bound(x);
    auto match_it = match_map.lower_bound(x);
    auto used_it = used_seeds.upper_bound(x);
    float m = x;
    Node *n = nullptr;
    float ub = *used_it;
    int index = a.get_index(x);
    while (curr_scores.size() > 0 and match_it->first < ub) {
        potential_seeds.erase(m);
        match_it++;
        m = match_it->first;
        while (recomb_it->first <= match_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            recombination_forward(r);
        }
        if (m >= a.coordinates[index + 1]) {
            forward_prune_states(a.coordinates[index + 1]);
            index = a.get_index(m);
        }
        n = get_node_at(m);
        mutation_update(n, m);
    }
    // cout << "Extend forward from " << x << " to " << m << endl;
}

void Trace_pruner::extend_backward(ARG &a, float x) {
    auto recomb_it = a.recombinations.upper_bound(x);
    auto match_it = match_map.find(x);
    auto used_it = used_seeds.upper_bound(x);
    recomb_it--;
    used_it--;
    float m = x;
    Node *n = nullptr;
    float lb = *used_it;
    int index = a.get_index(x);
    while (curr_scores.size() > 0 and match_it->first > lb) {
        potential_seeds.erase(m);
        match_it--;
        m = match_it->first;
        while (recomb_it->first > m) {
            Recombination &r = recomb_it->second;
            recomb_it--;
            recombination_backward(r);
        }
        if (m < a.coordinates[index]) {
            backward_prune_states(a.coordinates[index]);
            index = a.get_index(m);
        }
        n = get_node_at(m);
        mutation_update(n, m);
    }
    // cout << "Extend backward from " << x << " to " << m << endl;
}

void Trace_pruner::mutation_update(Node *n, float m) {
    float mismatch = 0;
    for (auto &[i, s] : curr_scores) {
        mismatch = count_mismatch(i.branch, n, m);
        s *= pow(mut_prob, mismatch);
    }
}

void Trace_pruner::recombination_forward(Recombination &r) {
    transition_scores.clear();
    for (auto &[i, s] : curr_scores) {
        forward_transition(r, i);
    }
}

void Trace_pruner::recombination_backward(Recombination &r) {
    
}

void Trace_pruner::forward_transition(Recombination &r, const Interval_info &interval) {
    float lb, ub;
    float w1, w2;
    Interval_info new_interval;
    Branch b;
    float p = curr_scores[interval];
    float l = interval.lb;
    float u = interval.ub;
    if (!r.affect(interval.branch)) {
        transition_scores[interval] = curr_scores[interval];
    } else if (interval.branch == r.source_branch) {
        lb = min(l, r.start_time);
        ub = min(u, r.start_time);
        w1 = exp_prop(interval.lb, interval.ub, lb, ub);
        new_interval = Interval_info(r.recombined_branch, lb, ub);
        transition_helper(new_interval, w1*p);
        lb = max(r.start_time, l);
        ub = max(r.start_time, u);
        w2 = exp_prop(l, u, lb, ub);
        new_interval = Interval_info(r.merging_branch, r.deleted_node->time, r.deleted_node->time);
        transition_helper(new_interval, w2*p);
    } else if (interval.branch == r.target_branch) {
        lb = min(l, r.inserted_node->time);
        ub = min(u, r.inserted_node->time);
        w1 = exp_prop(l, u, lb, ub);
        new_interval = Interval_info(r.lower_transfer_branch, lb, ub);
        transition_helper(new_interval, w1*p);
        lb = max(r.inserted_node->time, l);
        ub = max(r.inserted_node->time, u);
        w2 = exp_prop(l, u, lb, ub);
        new_interval = Interval_info(r.upper_transfer_branch, lb, ub);
        transition_helper(new_interval, w2*p);
    } else {
        new_interval = Interval_info(r.merging_branch, l, u);
        transition_helper(new_interval, p);
    }
}

void Trace_pruner::backward_transition(Recombination &r, const Interval_info &interval) {
    float lb, ub;
    float w1, w2;
    Interval_info new_interval;
    Branch b;
    float p = curr_scores[interval];
    float l = interval.lb;
    float u = interval.ub;
    if (!r.affect(interval.branch)) {
        transition_scores[interval] = curr_scores[interval];
    } else if (interval.branch == r.recombined_branch) {
        lb = min(l, r.start_time);
        ub = min(u, r.start_time);
        w1 = exp_prop(interval.lb, interval.ub, lb, ub);
        new_interval = Interval_info(r.source_branch, lb, ub);
        transition_helper(new_interval, w1*p);
        lb = max(r.start_time, l);
        ub = max(r.start_time, u);
        w2 = exp_prop(l, u, lb, ub);
        new_interval = Interval_info(r.target_branch, r.inserted_node->time, r.inserted_node->time);
        transition_helper(new_interval, w2*p);
    } else if (interval.branch == r.merging_branch) {
        lb = min(l, r.deleted_node->time);
        ub = min(u, r.deleted_node->time);
        w1 = exp_prop(l, u, lb, ub);
        new_interval = Interval_info(r.source_sister_branch, lb, ub);
        transition_helper(new_interval, w1*p);
        lb = max(r.deleted_node->time, l);
        ub = max(r.deleted_node->time, u);
        w2 = exp_prop(l, u, lb, ub);
        new_interval = Interval_info(r.source_parent_branch, lb, ub);
        transition_helper(new_interval, w2*p);
    } else {
        new_interval = Interval_info(r.target_branch, l, u);
        transition_helper(new_interval, p);
    }
}

void Trace_pruner::transition_helper(Interval_info &new_interval, float p) {
    if (p == 0) {
        return;
    }
    transition_scores[new_interval] += p;
}

void Trace_pruner::forward_prune_states(float x) {
    auto it = curr_scores.begin();
    while (it != curr_scores.begin()) {
        if (it->second < cutoff) {
            deleted_branches[x].insert(it->first.branch);
            it = curr_scores.erase(it);
            inserted_branches[x];
        } else {
            ++it;
        }
    }
}

void Trace_pruner::backward_prune_states(float x) {
    auto it = curr_scores.begin();
    while (it != curr_scores.begin()) {
        if (it->second < cutoff) {
            inserted_branches[x].insert(it->first.branch);
            it = curr_scores.erase(it);
            deleted_branches[x];
        } else {
            ++it;
        }
    }
}

float Trace_pruner::exp_prob(float l, float u) {
    float prob = exp(-l) - exp(-u);
    return prob;
}

float Trace_pruner::exp_prop(float l, float u, float x, float y) {
    float sub_prob = exp_prob(x, y);
    float prob = exp_prob(l, u);
    return sub_prob/prob;
}

float Trace_pruner::exp_median(float l, float u) {
    float lcdf = 1 - exp(-l);
    float ucdf = 1 - exp(-u);
    float cdf = 0.5*(lcdf + ucdf);
    float m = -log(1 - cdf);
    return m;
}
