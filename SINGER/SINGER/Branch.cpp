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

Branch::Branch(Node_ptr l, Node_ptr u) {
    assert((l == nullptr and u == nullptr) or l->time < u->time);
    lower_node = l;
    upper_node = u;
}

float Branch::length() {
    return upper_node->time - lower_node->time;
}

/*
bool Branch::operator<(const Branch &other) const {
    if (upper_node->time != other.upper_node->time) {
        return upper_node->time < other.upper_node->time;
    } else if (lower_node->time != other.lower_node->time) {
        return lower_node->time < other.lower_node->time;
    } else if (upper_node->index != other.upper_node->index) {
        return upper_node->index < other.upper_node->index;
    } else if (lower_node->index != other.lower_node->index) {
        return lower_node->index < other.lower_node->index;
    } else {
        return lower_node < other.lower_node;
    }
}
 */

bool Branch::operator<(const Branch &other) const {
    compare_node cn;
    if (cn(upper_node, other.upper_node)) {
        return true;
    } else if (cn(other.upper_node, upper_node)) {
        return false;
    } else if (cn(lower_node, other.lower_node)) {
        return true;
    } else if (cn(other.lower_node, lower_node)) {
        return false;
    } else {
        return false;
    }
}


bool Branch::operator==(const Branch &other) const {
    if (lower_node == other.lower_node and upper_node == other.upper_node) {
        return true;
    }
    return false;
}

bool Branch::operator!=(const Branch &other) const {
    if (lower_node != other.lower_node or upper_node != other.upper_node) {
        return true;
    }
    return false;
}
