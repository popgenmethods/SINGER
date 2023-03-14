//
//  Branch_node.cpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#include "Branch_node.hpp"

Branch_node::Branch_node(Branch b, int x) {
    branch = b;
    start_pos = x;
}

void Branch_node::add_prev_node(Branch_node *bn) {
    prev_nodes.insert(bn);
}

void Branch_node::mutation_update(float x, Node *n) {
    float s0 = n->get_state(x);
    float su = branch.upper_node->get_state(x);
    float sl = branch.lower_node->get_state(x);
    if (branch.upper_node->index == -1) { // the root branch can always tolerate a mutation
        ;
    } else if (abs(s0 - sl) > 0.5 and abs(s0 - su) > 0.5) {
        curr_mismatch += 1;
    }
}

bool operator<(const Branch_node& b, const Branch_node& c) {
    return b.branch < c.branch;
}
