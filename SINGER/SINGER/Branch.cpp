//
//  Branch.cpp
//  SINGER
//
//  Created by Yun Deng on 3/31/22.
//

#include "Branch.hpp"

Branch::Branch() {
    lower_node = nullptr;
    upper_node = nullptr;
}

Branch::Branch(Node *l, Node *u) {
    assert(l == nullptr or u == nullptr or l->time < u->time);
    lower_node = l;
    upper_node = u;
}

float Branch::length() {
    return upper_node->time - lower_node->time;
}

bool operator==(const Branch& b, const Branch& c) {
    return b.upper_node == c.upper_node and b.lower_node == c.lower_node;
}

bool operator!=(const Branch& b, const Branch& c) {
    return b.upper_node != c.upper_node or b.lower_node != c.lower_node;
}

bool operator<(const Branch& b, const Branch& c) {
    if (b.upper_node->time != c.upper_node->time) {
        return b.upper_node->time < c.upper_node->time;
    } else if (b.lower_node->time != c.lower_node->time) {
        return b.lower_node->time < c.lower_node->time;
    } else if (b.upper_node->index != c.upper_node->index) {
        return b.upper_node->index < c.upper_node->index;
    } else if (b.lower_node->index != c.lower_node->index) {
        return b.lower_node->index < c.lower_node->index;
    } else {
        return b.lower_node < c.lower_node;
    }
}
