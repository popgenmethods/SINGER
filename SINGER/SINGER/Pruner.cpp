//
//  Pruner.cpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#include "Pruner.hpp"

Pruner::Pruner() {}

void Pruner::prune(ARG a, map<int, Node *> lower_nodes) {
    
}

void Pruner::mutation_update(float x, Node *n) {
    for (Branch_node *bn : curr_nodes) {
        bn->mutation_update(x, n);
    }
}

void Pruner::mutation_update(set<float> mutations, Node *n) {
    for (float x : mutations) {
        mutation_update(x, n);
    }
}

void Pruner::recombination_update(Recombination &r) {
    int pos = r.pos;
    set<Branch_node *, compare_branch_node> new_nodes = {};
    Branch_node *nbn = nullptr;
    for (Branch_node *bn : curr_nodes) {
        if (!r.affect(bn->branch)) { // when unaffected by recombination, go to the same branch
            nbn = new Branch_node(bn->branch, pos);
            nbn->add_prev_node(bn);
            new_nodes.insert(nbn);
        } else if (bn->branch == r.source_branch) { // when source branch, go to the recombined branch, and the merging branch
            nbn = new Branch_node(r.recombined_branch, pos);
            new_nodes.insert(nbn);
            nbn = new Branch_node(r.merging_branch, pos);
            new_nodes.insert(nbn);
        } else if (bn->branch == r.source_branch) { // when target branch, go to lower/upper transfer and the recombined branch;
            nbn = new Branch_node(r.lower_transfer_branch, pos);
            new_nodes.insert(nbn);
            nbn = new Branch_node(r.upper_transfer_branch, pos);
            new_nodes.insert(nbn);
            nbn = new Branch_node(r.recombined_branch, pos);
            new_nodes.insert(nbn);
        } else { // when else branch, go to merging branch
            nbn = new Branch_node(r.merging_branch, pos);
            new_nodes.insert(nbn);
        }
    }
}

void Pruner::update_helper(Branch_node *bn, Branch b, int x) {
    
}
