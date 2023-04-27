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
