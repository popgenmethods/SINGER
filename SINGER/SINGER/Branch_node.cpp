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

bool operator<(const Branch_node& b, const Branch_node& c) {
    return b.branch < c.branch;
}
