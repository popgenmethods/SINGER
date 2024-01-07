//
//  Scaler.cpp
//  SINGER
//
//  Created by Yun Deng on 1/3/24.
//

#include "Scaler.hpp"

Scaler::Scaler() {}

void Scaler::compute_deltas(ARG &a) {
    // the base case
    map<Node_ptr, double, compare_node> node_start = {};
    map<Node_ptr, double, compare_node> node_span = {};
    for (const Branch &b : a.recombinations.begin()->second.inserted_branches) {
        node_start[b.upper_node] = 0;
    }
    node_start.erase(a.root);
    auto r_it = next(a.recombinations.begin());
    while (next(r_it) != a.recombinations.end()) {
        Recombination &r = r_it->second;
        Node_ptr dn = r.deleted_node;
        Node_ptr in = r.inserted_node;
        assert(node_start.count(dn) > 0);
        node_span[dn] += r.pos - node_start[dn];
        node_start.erase(dn);
        node_start[in] = r_it->first;
        r_it++;
    }
    for (auto &x : node_start) {
        Node_ptr n = x.first;
        node_span[n] += a.sequence_length - node_start[n];
    }
    // the root case
    map<Node_ptr, double, compare_node> root_start = {};
    Branch prev_branch;
    Branch next_branch;
    Recombination &r = a.recombinations.begin()->second;
    prev_branch = *(r.inserted_branches.rbegin());
    root_start[prev_branch.lower_node] = 0;
    r_it = next(a.recombinations.begin());
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
            node_span[dn] += r.pos - root_start[dn];
            root_start.erase(dn);
            root_start[in] = r.pos;
        } else {
            assert(next_branch.upper_node != a.root);
        }
        r_it++;
    }
    for (auto &x : root_start) {
        Node_ptr n = x.first;
        node_span[n] += a.sequence_length - x.second;
    }
    int num_samples = (int) a.sample_nodes.size();
    sorted_nodes.resize(node_span.size() + num_samples);
    node_deltas.resize(node_span.size() + num_samples);
    int index = num_samples;
    for (auto &x : node_span) {
        sorted_nodes[index] = x.first;
        node_deltas[index] = -x.second;
        index++;
    }
    index = 0;
    for (Node_ptr n : a.sample_nodes) {
        sorted_nodes[index] = n;
        index++;
    }
    double total_rate = accumulate(node_deltas.begin(), node_deltas.end(), 0.0);
    node_deltas[0] = -total_rate;
}

void Scaler::compute_old_grid() {
    expected_arg_length.resize(num_windows);
    rates.resize(sorted_nodes.size());
    accumulated_arg_length.resize(sorted_nodes.size());
    partial_sum(node_deltas.begin(), node_deltas.end(), rates.begin());
    assert(rates.back() > -1); // assures there is no negativitiy problem
    for (int i = 1; i < sorted_nodes.size(); i++) {
        accumulated_arg_length[i] = accumulated_arg_length[i-1] + rates[i-1]*(sorted_nodes[i]->time - sorted_nodes[i-1]->time);
    }
    double unit_arg_length = accumulated_arg_length.back()/num_windows;
    double partial_arg_length = 0;
    int new_index = 0;
    double rate = 0;
    double residue = 0;
    for (int i = 1; i <= num_windows; i++) {
        partial_arg_length = accumulated_arg_length.back()*i/num_windows;
        auto it = upper_bound(accumulated_arg_length.begin(), accumulated_arg_length.end(), partial_arg_length);
        new_index = (int) distance(accumulated_arg_length.begin(), it);
        new_index = min((int) sorted_nodes.size() - 1, new_index);
        rate = rates[new_index - 1];
        residue = accumulated_arg_length[new_index] - partial_arg_length;
        residue = max(0.0, residue);
        expected_arg_length[i-1] = unit_arg_length;
        old_grid.push_back(sorted_nodes[new_index]->time - residue/rate);
    }
    assert(old_grid.size() == num_windows + 1);
}

void Scaler::compute_new_grid(double theta) {
    for (auto &x : observed_arg_length) {
        x /= theta; // convert mutation counts back to arg length
    }
    double base_time = 0;
    double old_window_width = 0, scaling_factor = 0;
    for (int i = 1; i < old_grid.size(); i++) {
        old_window_width = old_grid[i] - old_grid[i-1];
        scaling_factor = observed_arg_length[i-1]/expected_arg_length[i-1];
        scaling_factors.push_back(scaling_factor);
        base_time += old_window_width*scaling_factor;
        new_grid.push_back(base_time);
    }
}

void Scaler::map_mutations(ARG &a) {
    observed_arg_length.resize(num_windows);
    for (auto &x : a.mutation_branches) {
        for (auto &y : x.second) {
            add_mutation(1.0, y.lower_node->time, y.upper_node->time);
        }
    }
}

void Scaler::add_mutation(double w, double lb, double ub) {
    double x, y, l;
    int index;
    auto it = upper_bound(old_grid.begin(), old_grid.end(), lb);
    it--;
    index = (int) distance(old_grid.begin(), it);
    while (old_grid[index] < ub) {
        x = old_grid[index];
        y = old_grid[index + 1];
        l = min(ub, y) - max(lb, x);
        observed_arg_length[index] += w*l/(ub - lb);
        index++;
    }
}

void Scaler::rescale(ARG &a, double theta) {
    compute_deltas(a);
    compute_old_grid();
    map_mutations(a);
    compute_new_grid(theta);
    int k = 0;
    double t;
    for (int i = 0; i < sorted_nodes.size(); i++) {
        while (sorted_nodes[i]->time > old_grid[k+1]) {
            k++;
        }
        t = scaling_factors[k]*(sorted_nodes[i]->time - old_grid[k]) + new_grid[k];
        sorted_nodes[i]->time = t;
    }
    for (int i = 0; i < sorted_nodes.size() - 1; i++) {
        assert(sorted_nodes[i]->time <= sorted_nodes[i+1]->time);
    }
}
