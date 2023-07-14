//
//  Normalizer.cpp
//  SINGER
//
//  Created by Yun Deng on 7/13/23.
//

#include "Normalizer.hpp"

Normalizer::Normalizer() {}

Normalizer::~Normalizer() {}

void Normalizer::get_root_span(ARG &a) {
    map<Node_ptr, float> root_start = {};
    Branch prev_branch;
    Branch next_branch;
    Recombination &r = a.recombinations.begin()->second;
    prev_branch = *(r.inserted_branches.rbegin());
    root_start[prev_branch.lower_node] = 0;
    auto r_it = a.recombinations.upper_bound(0);
    while (next(r_it)->first < a.sequence_length) {
        Recombination &r = r_it->second;
        prev_branch = *r.deleted_branches.rbegin();
        next_branch = *r.inserted_branches.rbegin();
        if (prev_branch.upper_node == a.root) {
            Node_ptr dn = prev_branch.lower_node;
            Node_ptr in = next_branch.lower_node;
            assert(dn != in);
            assert(next_branch.upper_node == a.root);
            assert(root_start.count(dn) > 0);
            root_span[dn] += r.pos - root_start[dn];
            root_start.erase(dn);
            root_start[in] = r.pos;
        } else {
            assert(next_branch.upper_node != a.root);
        }
        r_it++;
    }
    for (auto &x : root_start) {
        Node_ptr n = x.first;
        root_span[n] += a.sequence_length - x.second;
    }
}

void Normalizer::get_node_span(ARG &a) {
    map<Node_ptr, float> node_start = {};
    for (const Branch &b : a.recombinations.begin()->second.inserted_branches) {
        node_start[b.upper_node] = 0;
    }
    auto r_it = next(a.recombinations.begin());
    while (next(r_it) != a.recombinations.end()) {
        Node_ptr dn = r_it->second.deleted_node;
        Node_ptr in = r_it->second.inserted_node;
        assert(node_start.count(dn) > 0);
        node_span[dn] += r_it->first - node_start[dn];
        node_start.erase(dn);
        node_start[in] = r_it->first;
        r_it++;
    }
    for (auto &x : node_start) {
        Node_ptr n = x.first;
        node_span[n] += a.sequence_length - node_start[n];
    }
    for (auto &x : root_span) {
        Node_ptr n = x.first;
        node_span[n] += x.second;
    }
    node_span.erase(a.root);
}

void Normalizer::count_mutations(ARG &a) {
    for (auto &x : node_span) {
        mutation_counts[x.first->time] = 0;
    }
    mutation_counts[0] = 0;
    mutation_counts[INT_MAX] = 0;
    float lb, ub;
    for (auto &x : a.mutation_branches) {
        for (auto &y : x.second) {
            if (y.upper_node != a.root) {
                lb = y.lower_node->time;
                ub = y.upper_node->time;
                add_mutation(lb, ub);
            }
        }
    }
}

void Normalizer::randomize_mutation_ages(ARG &a) {
    float lb, ub, m;
    for (auto &x : a.mutation_branches) {
        for (auto &y : x.second) {
            if (y.upper_node != a.root) {
                lb = y.lower_node->time;
                ub = y.upper_node->time;
                m = uniform_random()*(ub - lb) + lb;
                mutation_ages.insert(m);
            }
        }
    }
    mutation_ages.insert(INT_MAX);
}

void Normalizer::randomized_normalize(ARG &a) {
    get_root_span(a);
    get_node_span(a);
    randomize_mutation_ages(a);
    float num_lineages = a.sample_nodes.size()*a.sequence_length;
    float theta = 0;
    float count = 0;
    float correction = 0;
    float tau = 0;
    auto n_it = node_span.begin();
    auto m_it = mutation_ages.begin();
    while (n_it != node_span.end()) {
        Node_ptr n = n_it->first;
        count = 0;
        while (*m_it < n->time) {
            count += 1;
            m_it++;
        }
        tau = max(count/num_lineages/theta, 1e-6f);
        /*
        if (correction + tau > 8) {
            correction += 1e-6;
        } else {
            correction += max(count/num_lineages/theta, 1e-6f);
        }
         */
        theta = num_lineages*4e-4;
        correction += max(count/theta, 1e-6f);
        n->time = correction;
        num_lineages -= n_it->second;
        n_it++;
    }
    a.smc_sample_recombinations();
}

void Normalizer::normalize(ARG &a) {
    get_root_span(a);
    get_node_span(a);
    count_mutations(a);
    float num_lineages = a.sample_nodes.size()*a.sequence_length;
    float theta = 4e-4;
    float count = 0;
    float correction = 0;
    auto n_it = node_span.begin();
    auto c_it = mutation_counts.begin();
    while (n_it != node_span.end()) {
        Node_ptr n = n_it->first;
        count = c_it->second;
        correction += max(count/num_lineages/theta, 1e-6f);
        n->time = correction;
        num_lineages -= n_it->second;
        n_it++;
        c_it++;
    }
    a.smc_sample_recombinations();
}

void Normalizer::add_mutation(float lb, float ub) {
    float x, y, l;
    auto it = mutation_counts.lower_bound(lb);
    while (it->first < ub) {
        x = it->first;
        y = next(it)->first;
        l = min(ub, y) - max(lb, x);
        mutation_counts[x] += l/(ub - lb);
        it++;
    }
}
