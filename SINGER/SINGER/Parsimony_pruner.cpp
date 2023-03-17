//
//  Parsimony_pruner.cpp
//  SINGER
//
//  Created by Yun Deng on 3/14/23.
//

#include "Parsimony_pruner.hpp"

Parsimony_pruner::Parsimony_pruner() {}

void Parsimony_pruner::start_search(Node *n, float m, set<Branch> branches) {
    int mismatch = 0;
    for (Branch branch : branches) {
        mismatch = count_mismatch(branch, n, m);
        curr_mismatch[branch] = mismatch;
    }
}

void Parsimony_pruner::mutation_forward(Node *n, float m) {
    int mismatch = 0;
    Branch branch = Branch();
    for (auto x : curr_mismatch) {
        branch = x.first;
        mismatch = x.second;
        if (mismatch > max_mismatch) {
            curr_mismatch.erase(branch);
        } else {
            mismatch += count_mismatch(branch, n, m);
            curr_mismatch[branch] = curr_mismatch[branch] + mismatch;
        }
    }
}

void Parsimony_pruner::recombination_forward(Recombination &r) {
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
