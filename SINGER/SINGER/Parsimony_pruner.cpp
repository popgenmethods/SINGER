//
//  Parsimony_pruner.cpp
//  SINGER
//
//  Created by Yun Deng on 3/14/23.
//

#include "Parsimony_pruner.hpp"

Parsimony_pruner::Parsimony_pruner() {
}

void Parsimony_pruner::prune_arg(ARG &a, map<float, Node *> base_nodes) {
    nodes = base_nodes;
    start = base_nodes.begin()->first;
    end = base_nodes.rbegin()->first;
    used_seeds = {start, end};
    build_match_map(a, base_nodes);
    float x = 0;
    while (potential_seeds.size() > 2) {
        x = find_minimum_match();
        extend(a, x);
    }
    deleted_branches.insert({a.sequence_length, {}});
    inserted_branches.insert({a.sequence_length, {}});
    write_reductions(a);
    cout << endl;
}

void Parsimony_pruner::start_search(ARG &a, float m) {
    curr_mismatch.clear();
    Node *n = get_node_at(m);
    float mismatch = 0;
    Tree tree = a.get_tree_at(m);
    for (Branch branch : tree.branches) {
        mismatch = count_mismatch(branch, n, m);
        if (mismatch == 0 and branch.lower_node != n) {
            curr_mismatch[branch] = mismatch;
        }
    }
    potential_seeds.erase(m);
}

void Parsimony_pruner::extend(ARG &a, float x) {
    start_search(a, x);
    cout << "Seed search size: " << curr_mismatch.size() << endl;
    extend_forward(a, x);
    start_search(a, x);
    extend_backward(a, x);
    used_seeds.insert(x); // when extending later seeds, don't go beyond previous seeds (to save computation)
}

void Parsimony_pruner::mutation_forward(Node *n, float m) {
    float mismatch = 0;
    Branch branch = Branch();
    set<Branch> bad_branches = {};
    for (auto x : curr_mismatch) {
        branch = x.first;
        mismatch = x.second;
        mismatch += count_mismatch(branch, n, m);
        curr_mismatch[branch] = mismatch;
        if (mismatch > max_mismatch) {
            bad_branches.insert(branch);
        }
    }
    for (Branch b : bad_branches) {
        curr_mismatch.erase(b);
    }
    if (bad_branches.size() > 0) {
        write_reduction_change(m, bad_branches, {});
    }
}

void Parsimony_pruner::mutation_backward(Node *n, float m) {
    float mismatch = 0;
    Branch branch = Branch();
    set<Branch> bad_branches = {};
    for (auto x : curr_mismatch) {
        branch = x.first;
        mismatch = x.second;
        mismatch += count_mismatch(branch, n, m);
        curr_mismatch[branch] = mismatch;
        if (mismatch > max_mismatch) {
            bad_branches.insert(branch);
        }
    }
    for (Branch b : bad_branches) {
        curr_mismatch.erase(b);
    }
    if (bad_branches.size() > 0) {
        write_reduction_change(m, {}, bad_branches);
    }
    if (m == start) {
        write_init_set();
    }
}

Node *Parsimony_pruner::get_node_at(float x) {
    map<float, Node *>::iterator node_it = nodes.upper_bound(x);
    node_it--;
    return node_it->second;
}

void Parsimony_pruner::recombination_forward(Recombination &r) {
    transitions.clear();
    Branch branch = Branch();
    set<Branch> db = {};
    set<Branch> ib = {};
    for (auto x : curr_mismatch) {
        branch = x.first;
        if (not r.affect(branch)) {
            transition_helper(branch, branch);
        } else if (branch == r.source_branch) {
            transition_helper(branch, r.recombined_branch);
            db.insert(branch);
            ib.insert(r.recombined_branch);
        } else if (branch == r.target_branch) {
            transition_helper(branch, r.lower_transfer_branch);
            transition_helper(branch, r.upper_transfer_branch);
            db.insert(branch);
            ib.insert(r.lower_transfer_branch);
            ib.insert(r.upper_transfer_branch);
        } else {
            transition_helper(branch, r.merging_branch);
            db.insert(branch);
            ib.insert(r.merging_branch);
        }
    }
    write_reduction_change(r.pos, db, ib);
    update_mismatch();
}

void Parsimony_pruner::recombination_backward(Recombination &r) {
    transitions.clear();
    Branch branch = Branch();
    set<Branch> db = {};
    set<Branch> ib = {};
    for (auto x : curr_mismatch) {
        branch = x.first;
        if (not r.create(branch)) {
            transition_helper(branch, branch);
        }
        if (branch == r.recombined_branch) {
            transition_helper(branch, r.source_branch);
            ib.insert(branch);
            db.insert(r.source_branch);
        }
        if (branch == r.lower_transfer_branch) {
            transition_helper(branch, r.target_branch);
            ib.insert(branch);
            db.insert(r.target_branch);
        }
        if (branch == r.upper_transfer_branch) {
            transition_helper(branch, r.target_branch);
            ib.insert(branch);
            db.insert(r.target_branch);
        }
        if (branch == r.merging_branch) {
            transition_helper(branch, r.source_sister_branch);
            transition_helper(branch, r.source_parent_branch);
            ib.insert(branch);
            db.insert(r.source_sister_branch);
            db.insert(r.source_parent_branch);
        }
    }
    write_reduction_change(r.pos, db, ib);
    update_mismatch();
}

void Parsimony_pruner::check_reduction(map<float, pair<Branch, Node *>> joining_points) {
      float start = joining_points.begin()->first;
    float end = joining_points.rbegin()->first;
    map<float, pair<Branch, Node *>>::iterator join_it = joining_points.begin();
    map<float, set<Branch>>::iterator reduced_it = reductions.begin();
    Branch joining_branch = Branch();
    float x = start;
    set<Branch> reduced_set = reduced_it->second;
    while (join_it->first < end) {
        x = join_it->first;
        joining_branch = join_it->second.first;
        while (reduced_it->first < x) {
            reduced_set = reduced_it->second;
            reduced_it++;
        }
        cout << "Position " << join_it->first << " includes truth: " << reduced_set.count(joining_branch) << endl;
        join_it++;
    }
}

void Parsimony_pruner::print_reduction_size() {
    for (auto x : reductions) {
        cout << x.first << " : " << x.second.size() << endl;
    }
}

void Parsimony_pruner::build_match_map(ARG &a, map<float, Node *> base_nodes) {
    float m = 0;
    float inf = INT_MAX;
    float lb = 0;
    float state = 0;
    map<float, set<Branch>>::iterator mb_it = a.mutation_branches.lower_bound(start);
    Node *n = nullptr;
    while (mb_it->first < end) {
        m = mb_it->first;
        n = get_node_at(m);
        state = n->get_state(m);
        lb = -1;
        for (Branch b : mb_it->second) {
            if (b.lower_node->get_state(m) == state) {
                lb = max(lb, b.lower_node->time);
            } else if (b.upper_node->index == -1) {
                lb = max(lb, 0.0f);
            } else {
                lb = max(lb, max_time - b.lower_node->time);
            }
        }
        if (lb >= 0) {
            match_map[m] = lb;
        }
        mb_it++;
    }
    match_map[end] = inf;
    match_map[start] = inf;
    potential_seeds = match_map;
}

float Parsimony_pruner::find_minimum_match() {
    auto it = min_element(potential_seeds.begin(), potential_seeds.end(),
                          [](const auto& l, const auto& r) { return l.second < r.second;});
    float x = it->first;
    cout << "Seed mutation age: " << it->second << endl;
    return x;
}

float Parsimony_pruner::count_mismatch(Branch branch, Node *n, float m) {
    if (private_mutations.count(m) > 0) { // special penalty for singleton mutations
        return 0.5;
    }
    float s0 = n->get_state(m);
    float sl = branch.lower_node->get_state(m);
    float su = branch.upper_node->get_state(m);
    if (abs(sl - s0) > 0.5 and abs(su - s0) > 0.5) {
        return 1;
    } else {
        return 0;
    }
}

void Parsimony_pruner::transition_helper(Branch sb, Branch tb) {
    if (transitions.count(tb) > 0) {
        transitions[tb].insert(sb);
    } else {
        transitions[tb] = {sb};
    }
}

void Parsimony_pruner::update_mismatch() {
    map<Branch, float> prev_mismatch = curr_mismatch;
    map<Branch, float> next_mismatch = {};
    Branch branch = Branch();
    float mismatch_sum = 0;
    for (auto x : transitions) {
        branch = x.first;
        mismatch_sum = 0;
        for (Branch b : x.second) {
            mismatch_sum += prev_mismatch[b];
        }
        next_mismatch[branch] = mismatch_sum/x.second.size();
    }
    curr_mismatch = next_mismatch;
}

void Parsimony_pruner::write_init_set() {
    set<Branch> init_branches = {};
    for (auto x : curr_mismatch) {
        init_branches.insert(x.first);
    }
    deleted_branches[0] = {};
    inserted_branches[0] = init_branches;
}

void Parsimony_pruner::write_reduction_change(float x, set<Branch> db, set<Branch> ib) {
    if (db.size() == 0 and ib.size() == 0) { // skip the no information change
        return;
    }
    if (deleted_branches.count(x) > 0) {
        deleted_branches[x].insert(db.begin(), db.end());
    } else {
        deleted_branches[x] = db;
    }
    if (inserted_branches.count(x) > 0) {
        inserted_branches[x].insert(ib.begin(), ib.end());
    } else {
        inserted_branches[x] = ib;
    }
}

void Parsimony_pruner::write_reductions(ARG &a) {
    map<float, set<Branch>>::iterator deleted_it = deleted_branches.begin();
    map<float, set<Branch>>::iterator inserted_it = inserted_branches.begin();
    vector<float>::iterator start_it = std::lower_bound(a.coordinates.begin(), a.coordinates.end(), start);
    int start_index = (int) (start_it - a.coordinates.begin());
    vector<float>::iterator end_it = std::lower_bound(a.coordinates.begin(), a.coordinates.end(), end);
    int end_index = (int) (end_it - a.coordinates.begin());
    set<Branch> reduced_set = {};
    for (int i = start_index; i < end_index; i++) {
        while (deleted_it->first < a.coordinates[i+1]) {
            for (auto x : deleted_it->second) {
                reduced_set.erase(x);
            }
            for (auto x : inserted_it->second) {
                reduced_set.insert(x);
            }
            deleted_it++;
            inserted_it++;
            reductions[a.coordinates[i]] = reduced_set;
        }
    }
}

void Parsimony_pruner::extend_forward(ARG &a, float x) {
    map<float, Recombination>::iterator recomb_it = a.recombinations.upper_bound(x);
    map<float, float>::iterator match_it = match_map.lower_bound(x);
    set<float>::iterator used_it = used_seeds.upper_bound(x);
    float m = x;
    Node *n = nullptr;
    float ub = *used_it;
    while (curr_mismatch.size() > 0 and match_it->first < ub) {
        potential_seeds.erase(m);
        match_it++;
        m = match_it->first;
        while (recomb_it->first <= m) {
            Recombination r = recomb_it->second;
            recomb_it++;
            recombination_forward(r);
        }
        n = get_node_at(m);
        mutation_forward(n, m);
    }
    cout << "Extend forward from " << x << " to " << m << endl;
}

void Parsimony_pruner::extend_backward(ARG &a, float x) {
    map<float, Recombination>::iterator recomb_it = a.recombinations.upper_bound(x);
    map<float, float>::iterator match_it = match_map.find(x);
    set<float>::iterator used_it = used_seeds.upper_bound(x);
    recomb_it--;
    used_it--;
    float m = x;
    Node *n = nullptr;
    float lb = *used_it;
    while (curr_mismatch.size() > 0 and match_it->first > lb) {
        potential_seeds.erase(m);
        match_it--;
        m = match_it->first;
        while (recomb_it->first > m) {
            Recombination r = recomb_it->second;
            recomb_it--;
            recombination_backward(r);
        }
        n = get_node_at(m);
        mutation_backward(n, m);
    }
    cout << "Extend backward from " << x << " to " << m << endl;
}
