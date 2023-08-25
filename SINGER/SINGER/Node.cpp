//
//  Node.cpp
//  SINGER
//
//  Created by Yun Deng on 3/31/22.
//

#include "Node.hpp"

Node::Node(float t) {
    time = t;
}

void Node::set_index(int index) {
    this->index = index;
}

/*
void Node::add_mutation(float pos) {
    mutation_sites.insert(pos);
}
 */

void Node::add_mutation(float pos) {
    mutation_sites[pos] = 1;
}


/*
float Node::get_state(float pos) {
    if (mutation_sites.count(pos) > 0) {
        return 1;
    } else if (ambiguous_sites.count(pos) > 0) {
        return 0.5;
    } else {
        return 0;
    }
}

void Node::write_state(float pos, float s) {
    if (s == 0) {
        mutation_sites.erase(pos);
        ambiguous_sites.erase(pos);
        return;
    } else if (s == 0.5) {
        mutation_sites.erase(pos);
        ambiguous_sites.insert(pos);
    } else if (s == 1) {
        ambiguous_sites.erase(pos);
        mutation_sites.insert(pos);
    }
    return;
}
 */

float Node::get_state(float pos) {
    auto it = mutation_sites.find(pos);
    if (it == mutation_sites.end()) {
        return 0;
    } else {
        return it->second;
    }
}

void Node::write_state(float pos, float s) {
    auto it = mutation_sites.find(pos);
    if (s == 0) {
        if (it != mutation_sites.end()) {
            mutation_sites.erase(pos);
        }
    } else {
        if (it == mutation_sites.end()) {
            mutation_sites[pos] = s;
        } else {
            it->second = s;
        }
    }
    return;
}

/*
float Node::get_state(float pos) {
    if (it == mutation_sites.begin()) {
        it = mutation_sites.upper_bound(pos);
        it--;
    } else if (*prev(it) > pos or *next(it) < pos) {
        it = mutation_sites.upper_bound(pos);
        it--;
    }
    if (*it == pos) {
        return 1;
    } else if (pos > *it) {
        if (*next(it) == pos) {
            it++;
            return 1;
        }
    } else if (pos < *it) {
        if (*prev(it) == pos) {
            it--;
            return 1;
        }
    }
    return 0;
}

void Node::write_state(float pos, float s) {
    if (s == 0) {
        mutation_sites.erase(pos);
        return;
    } else if (s == 1) {
        mutation_sites.insert(pos);
    }
    return;
}
 */

void Node::read_mutation(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    float x;
    while (fin >> x) {
        add_mutation(x);
    }
}

shared_ptr<Node> new_node(float t) {
    return make_shared<Node>(t);
}
