//
//  Normalizer.cpp
//  SINGER
//
//  Created by Yun Deng on 7/13/23.
//

#include "Normalizer.hpp"

Normalizer::Normalizer() {}

Normalizer::~Normalizer() {}

void Normalizer::compute_max_time(ARG &a) {
    double num_muts = a.mutation_sites.size();
    double theta = accumulate(a.thetas.begin(), a.thetas.end(), 0.0f);
    double last_mut_pos = *prev(a.mutation_sites.end(), 2);
    double first_mut_pos = *a.mutation_sites.begin();
    double q = (last_mut_pos - first_mut_pos)/a.sequence_length;
    theta *= q;
    double num_samples = a.sample_nodes.size();
    double scaling = num_muts/2/theta/log(num_samples);
    max_time = 10*scaling;
}

void Normalizer::get_root_span(ARG &a) {
    map<Node_ptr, double, compare_node> root_start = {};
    map<Node_ptr, double, compare_node> root_span = {};
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
    all_root_nodes.resize(root_span.size());
    all_root_spans.resize(root_span.size());
    int index = 0;
    for (auto &x : root_span) {
        all_root_nodes[index] = x.first;
        all_root_spans[index] = x.second;
        index++;
    }
}

void Normalizer::get_node_span(ARG &a) {
    map<Node_ptr, double, compare_node> node_start = {};
    map<Node_ptr, double, compare_node> node_span = {};
    for (const Branch &b : a.recombinations.begin()->second.inserted_branches) {
        node_start[b.upper_node] = 0;
    }
    node_start.erase(a.root);
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
    all_nodes.resize(node_span.size());
    all_spans.resize(node_span.size());
    int index = 0;
    for (auto &x : node_span) {
        all_nodes[index] = x.first;
        all_spans[index] = x.second;
        index++;
    }
}

void Normalizer::randomize_mutation_ages(ARG &a) {
    double lb, ub, m;
    for (auto &x : a.mutation_branches) {
        for (auto &y : x.second) {
            if (y.upper_node != a.root) {
                lb = y.lower_node->time;
                ub = y.upper_node->time;
                m = uniform_random()*(ub - lb) + lb;
                // m = 0.5*(lb + ub);
                mutation_ages.insert(m);
            }
        }
    }
    mutation_ages.insert(INT_MAX);
}

void Normalizer::randomized_normalize(ARG &a) {}

void Normalizer::normalize(ARG &a, double theta) {
    compute_max_time(a);
    get_root_span(a);
    get_node_span(a);
    partition_arg(a);
    count_mutations(a);
    count_recombinations(a);
    double t = 1.0/a.Ne;
    int k = 0;
    new_grid.resize(old_grid.size());
    new_grid[0] = t;
    for (int i = 0; i < new_grid.size() - 1; i++) {
        double e = expected_mutation_counts[i];
        double o = observed_mutation_counts[i]/theta;
        double scale = o/e;
        new_grid[i+1] = new_grid[i] + scale*(old_grid[i+1] - old_grid[i]);
        new_grid[i+1] = min(max_time, new_grid[i+1]); // process grid values
        new_grid[i+1] = max(new_grid[i+1], nextafter(new_grid[i], numeric_limits<double>::infinity())); // process grid values
    }
    int i = 0;
    double p;
    k = 0;
    while (k < new_grid.size() - 1) {
        while (i < all_nodes.size() and all_nodes[i]->time <= old_grid[k+1]) {
            Node_ptr node = all_nodes[i];
            p = (node->time - old_grid[k])/(old_grid[k+1] - old_grid[k]);
            t = (1 - p)*new_grid[k] + p*new_grid[k+1];
            if (i > 0 and t <= all_nodes[i-1]->time) {
                t = nextafter(all_nodes[i-1]->time, numeric_limits<double>::infinity());
                assert(t > all_nodes[i-1]->time); // don't want duplicated coalescence time
            }
            node->time = t;
            i++;
        }
        k++;
    }
    for (auto &x : a.recombinations) {
        x.second.start_time = -1;
    }
    a.adjust_recombinations();
}

void Normalizer::partition_arg(ARG &a) {
    double active_lineages = accumulate(all_spans.begin(), all_spans.end(), 0);
    active_lineages += accumulate(all_root_spans.begin(), all_root_spans.end(), 0);
    double t = 0;
    int j = 0;
    for (int i = 0; i < all_spans.size(); i++) {
        Node_ptr n = all_nodes[i];
        ls += (n->time - t)*active_lineages;
        active_lineages -= all_spans[i];
        if (all_root_nodes[j] == n) {
            active_lineages -= all_root_spans[j];
            j++;
        }
        t = n->time;
    }
    double l = 0;
    double x = 0;
    t = 0;
    j = 0;
    int k = 1;
    active_lineages = accumulate(all_spans.begin(), all_spans.end(), 0);
    active_lineages += accumulate(all_root_spans.begin(), all_root_spans.end(), 0);
    old_grid.push_back(0);
    for (int i = 0; i < all_spans.size(); i++) {
        Node_ptr n = all_nodes[i];
        l += (n->time - t)*active_lineages;
        while (l > k*ls/num_windows and k <= num_windows - 1) {
            x = n->time - (l - k*ls/num_windows)/active_lineages;
            old_grid.push_back(x);
            expected_mutation_counts.push_back(ls/num_windows);
            k++;
        }
        active_lineages -= all_spans[i];
        if (all_root_nodes[j] == n) {
            active_lineages -= all_root_spans[j];
            j++;
        }
        t = n->time;
    }
    old_grid.push_back(all_nodes.back()->time);
    expected_mutation_counts.push_back(l - (1 - 1.0/num_windows)*ls);
}

void Normalizer::normalize_recombinations(ARG &a) {
    count_recombinations(a);
    sample_recombinations(a);
}

void Normalizer::count_mutations(ARG &a) {
    observed_mutation_counts.resize(expected_mutation_counts.size());
    double lb, ub;
    for (auto &x : a.mutation_branches) {
        for (auto &y : x.second) {
            if (y.upper_node != a.root) {
                lb = y.lower_node->time;
                ub = y.upper_node->time;
                add_mutation(lb, ub);
            }
        }
    }
    double num_muts = a.mutation_sites.size() - 1;
    double num_mapped_muts = accumulate(observed_mutation_counts.begin(), observed_mutation_counts.end(), 0.0);
    double r = num_muts/num_mapped_muts;
    for (auto &x : observed_mutation_counts) {
        x *= r;
    }
}

void Normalizer::count_recombinations(ARG &a) {
    observed_recombination_counts.resize(observed_mutation_counts.size());
    double lb, ub;
    for (auto &x : a.recombinations) {
        if (x.first > 0 and x.first < a.sequence_length) {
            Recombination &r = x.second;
            lb = r.source_branch.lower_node->time;
            ub = min(r.inserted_node->time, r.deleted_node->time);
            add_recombination(lb, ub);
        }
    }
    double num_recs = accumulate(observed_recombination_counts.begin(), observed_recombination_counts.end(), 0.0f);
    num_recs /= observed_recombination_counts.size();
    recombination_density.resize(observed_recombination_counts.size());
    for (int i = 0; i < observed_recombination_counts.size(); i++) {
        recombination_density[i] = num_recs/observed_recombination_counts[i];
    }
}

void Normalizer::add_mutation(double lb, double ub) {
    double x, y, l;
    int index;
    auto it = upper_bound(old_grid.begin(), old_grid.end(), lb);
    it--;
    index = (int) distance(old_grid.begin(), it);
    while (old_grid[index] < ub) {
        x = old_grid[index];
        y = old_grid[index + 1];
        l = min(ub, y) - max(lb, x);
        observed_mutation_counts[index] += l/(ub - lb);
        index++;
    }
}

void Normalizer::add_recombination(double lb, double ub) {
    double x, y, l;
    int index;
    auto it = upper_bound(old_grid.begin(), old_grid.end(), lb);
    it--;
    index = (int) distance(old_grid.begin(), it);
    while (old_grid[index] < ub) {
        x = old_grid[index];
        y = old_grid[index + 1];
        l = min(ub, y) - max(lb, x);
        observed_recombination_counts[index] += l/(ub - lb);
        index++;
    }
}

/*
void Normalizer::calculate_branch_length(ARG &a, double Ne) {
    observed_branch_length.resize(expected_mutation_counts.size());
    expected_branch_length.resize(expected_mutation_counts.size());
    double all_lineages = accumulate(all_spans.begin(), all_spans.end(), 0);
    double t = 0;
    int k = 0;
    double lb = 0, ub = 0;
    for (int i = 0; i < all_spans.size(); i++) {
        Node_ptr n = all_nodes[i];
        lb = max(t, old_grid[k]);
        ub = min(n->time, old_grid[k+1]);
        observed_branch_length[k] += (ub - lb)*all_lineages;
        while (n->time > old_grid[k+1] and k <= num_windows - 2) {
            k++;
            lb = max(t, old_grid[k]);
            ub = min(n->time, old_grid[k+1]);
            observed_branch_length[k] += all_lineages*(ub - lb);
        }
        all_lineages -= all_spans[i];
        t = n->time;
    }
    for (k = 0; k < old_grid.size() - 1; k++) {
        lb = a.sample_nodes.size()*(exp(0.5*old_grid[k]) - 1) + 1;
        ub = a.sample_nodes.size()*(exp(0.5*old_grid[k+1]) - 1) + 1;
        expected_branch_length[k] = 2*(log(ub) - log(lb))*a.sequence_length;
    }
}
*/

double Normalizer::sample_recombination_time(double lb, double ub) {
    lb = max(lb, new_grid.front());
    double tau;
    double x, y, l;
    double q = 0;
    int base_index, index;
    auto it = upper_bound(new_grid.begin(), new_grid.end(), lb);
    it--;
    base_index = (int) distance(new_grid.begin(), it);
    index = base_index;
    while (new_grid[index] < ub) {
        x = new_grid[index];
        y = new_grid[index + 1];
        l = min(ub, y) - max(lb, x);
        q += l/(y - x)*recombination_density[index];
        index++;
    }
    double rq = q*uniform_random();
    index = base_index;
    while (new_grid[index] < ub) {
        x = new_grid[index];
        y = new_grid[index + 1];
        l = min(ub, y) - max(lb, x);
        rq -= l/(y - x)*recombination_density[index];
        if (rq < 0) {
            tau = min(ub, y) + rq*l/recombination_density[index];
            return tau;
        }
        index++;
    }
    cerr << "resample recombination time failed" << endl;
    exit(1);
}

void Normalizer::sample_recombinations(ARG &a) {
    double lb, ub;
    for (auto &x : a.recombinations) {
        if (x.first > 0 and x.first < a.sequence_length) {
            Recombination &r = x.second;
            lb = r.source_branch.lower_node->time;
            ub = min(r.inserted_node->time, r.deleted_node->time);
            r.start_time = nextafter(ub, lb);
            assert(r.start_time > lb and r.start_time < ub);
        }
    }
}
