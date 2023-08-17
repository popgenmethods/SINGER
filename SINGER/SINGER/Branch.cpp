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

bool Branch::operator<(const Branch &other) const {
    static compare_node cn;
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
