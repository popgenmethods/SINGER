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
    used_seeds = {0.0, a.sequence_length};
    get_match_map(a, base_nodes);
    float x = 0;
    while (potential_seeds.size() > 2) {
        x = find_minimum_match();
        extend(a, x);
    }
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
    cout << "Search size: " << curr_mismatch.size() << endl;
    extend_forward(a, x);
    start_search(a, x);
    cout << "Search size: " << curr_mismatch.size() << endl;
    extend_backward(a, x);
    used_seeds.insert(x);
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
}

Node *Parsimony_pruner::get_node_at(float x) {
    map<float, Node *>::iterator node_it = nodes.upper_bound(x);
    node_it--;
    return node_it->second;
}

void Parsimony_pruner::recombination_forward(Recombination &r) {
    transitions.clear();
    Branch branch = Branch();
    for (auto x : curr_mismatch) {
        branch = x.first;
        if (not r.affect(branch)) {
            transition_helper(branch, branch);
        } else if (branch == r.source_branch) {
            transition_helper(branch, r.recombined_branch);
            transition_helper(branch, r.merging_branch);
        } else if (branch == r.target_branch) {
            transition_helper(branch, r.lower_transfer_branch);
            transition_helper(branch, r.upper_transfer_branch);
            transition_helper(branch, r.recombined_branch);
        } else {
            transition_helper(branch, r.merging_branch);
        }
    }
    update_mismatch();
}

void Parsimony_pruner::recombination_backward(Recombination &r) {
    transitions.clear();
    Branch branch = Branch();
    for (auto x : curr_mismatch) {
        branch = x.first;
        if (not r.create(branch)) {
            transition_helper(branch, branch);
        } else if (branch == r.recombined_branch) {
            transition_helper(branch, r.source_branch);
            transition_helper(branch, r.target_branch);
        } else if (branch == r.lower_transfer_branch) {
            transition_helper(branch, r.target_branch);
        } else if (branch == r.upper_transfer_branch) {
            transition_helper(branch, r.target_branch);
        } else {
            for (Branch b : r.deleted_branches) {
                if (b != r.source_branch and b != r.target_branch) {
                    transition_helper(branch, b);
                }
            }
        }
    }
    update_mismatch();
}

void Parsimony_pruner::get_match_map(ARG &a, map<float, Node *> base_nodes) {
    float start = base_nodes.begin()->first;
    float end = base_nodes.rbegin()->first;
    map<float, Recombination>::iterator recomb_it = a.recombinations.upper_bound(start);
    set<float>::iterator mut_it = a.mutation_sites.lower_bound(start);
    map<float, Node *>::iterator base_it = base_nodes.begin();
    float m = *mut_it;
    Tree tree = a.get_tree_at(start);
    Recombination r;
    Node *base_node = nullptr;
    float match = 0;
    while (*mut_it < end) {
        while (base_it->first < m) {
            base_node = base_it->second;
            base_it++;
        }
        while (recomb_it->first < m) {
            r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
        m = *mut_it;
        match = 0;
        mut_it++;
        for (Branch b : tree.branches) {
            if (b.lower_node != base_node and b.upper_node->time > base_node->time) {
                match += 1 - count_mismatch(b, base_node, m);
            }
        }
        if (match > 0) {
            match_map[m] = match;
        }
    }
    match_map[INT_MAX] = INT_MAX;
    match_map[-1.0] = INT_MAX;
    potential_seeds = match_map;
}

float Parsimony_pruner::find_minimum_match() {
    auto it = min_element(potential_seeds.begin(), potential_seeds.end(),
                          [](const auto& l, const auto& r) { return l.second < r.second;});
    float x = it->first;
    return x;
}

float Parsimony_pruner::count_mismatch(Branch branch, Node *n, float m) {
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

void Parsimony_pruner::extend_forward(ARG &a, float x) {
    map<float, Recombination>::iterator recomb_it = a.recombinations.upper_bound(x);
    map<float, float>::iterator match_it = match_map.upper_bound(x);
    set<float>::iterator used_it = used_seeds.upper_bound(x);
    float m = x;
    Node *n = nullptr;
    float ub = *used_it;
    while (curr_mismatch.size() > 0 and match_it->first < ub) {
        m = match_it->first;
        while (recomb_it->first < m) {
            Recombination r = recomb_it->second;
            recomb_it++;
            recombination_forward(r);
        }
        n = get_node_at(m);
        mutation_forward(n, m);
        potential_seeds.erase(m);
        match_it++;
    }
    cout << "Forward search from " << x << " to " << m << endl;
}

void Parsimony_pruner::extend_backward(ARG &a, float x) {
    map<float, Recombination>::iterator recomb_it = a.recombinations.upper_bound(x);
    map<float, float>::iterator match_it = match_map.find(x);
    set<float>::iterator used_it = used_seeds.upper_bound(x);
    recomb_it--;
    match_it--;
    used_it--;
    float m = x;
    Node *n = nullptr;
    float lb = *used_it;
    while (curr_mismatch.size() > 0 and match_it->first > lb) {
        m = match_it->first;
        while (recomb_it->first >= m) {
            Recombination r = recomb_it->second;
            recomb_it--;
            recombination_backward(r);
        }
        n = get_node_at(m);
        mutation_backward(n, m);
        potential_seeds.erase(m);
        match_it--;
    }
    cout << "Backward search from " << x << " to " << m << endl;
}
