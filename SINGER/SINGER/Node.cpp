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

void Node::add_mutation(float pos) {
    mutation_sites.insert(pos);
}

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
