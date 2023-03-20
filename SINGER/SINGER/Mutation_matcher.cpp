//
//  Mutation_matcher.cpp
//  SINGER
//
//  Created by Yun Deng on 3/16/23.
//

#include "Mutation_matcher.hpp"

Mutation_matcher::Mutation_matcher(ARG &a) {
    float start = 0;
    float end = a.sequence_length;
    map<float, Recombination>::iterator recomb_it = a.recombinations.lower_bound(0);
    set<float>::iterator mut_it = a.mutation_sites.lower_bound(0);
    Tree tree = Tree();
    float curr_pos = start;
    float x = 0;
    while (curr_pos < end) {
        Recombination r = recomb_it->second;
        tree.forward_update(r);
        recomb_it++;
        curr_pos = recomb_it->first;
        while (*mut_it < curr_pos) {
            x = *mut_it;
            classify_branches(tree, x);
            mut_it++;
        }
    }
}

set<Branch> Mutation_matcher::get_matches(float x, float state) {
    assert(state == 0 or state == 1);
    if (state == 0) {
        return ancestral_branches.at(x);
    } else {
        return derived_branches.at(x);
    }
}

float Mutation_matcher::count_match(float x, float state) {
    float match_count = 0;
    if (state == 0) {
        match_count = ancestral_branches.at(x).size();
    } else if (state == 1) {
        match_count = derived_branches.at(x).size();
    } else {
        match_count = ancestral_branches.at(x).size() + derived_branches.at(x).size();
    }
    return match_count;
}


void Mutation_matcher::classify_branches(Tree tree, float x) {
    set<Branch> derived = {};
    set<Branch> ancestral = {};
    float sl = 0;
    float su = 0;
    for (Branch b : tree.branches) {
        sl = b.lower_node->get_state(x);
        su = b.upper_node->get_state(x);
        if (sl != su or sl == 0.5 or su == 0.5) {
            derived.insert(b);
            ancestral.insert(b);
        } else if (sl == 0) {
            ancestral.insert(b);
        } else {
            derived.insert(b);
        }
    }
    ancestral_branches[x] = ancestral;
    derived_branches[x] = derived;
}

map<float, float> Mutation_matcher::build_match_map(map<float, Node *> base_nodes) {
    map<float, float> match_map = {};
    float state = 0;
    float match = 0;
    float start = base_nodes.begin()->first;
    float end = base_nodes.rbegin()->first;
    map<float, set<Branch>>::iterator anc_it = ancestral_branches.lower_bound(start);
    map<float, Node *>::iterator base_it = base_nodes.begin();
    Node *base_node = nullptr;
    float m = anc_it->first;
    while (m < end) {
        while (base_it->first < m) {
            base_node = base_it->second;
            base_it++;
        }
        state = base_node->get_state(m);
        match = count_match(m, state);
        match_map[m] = match;
        anc_it++;
        m = anc_it->first;
    }
    return match_map;
}
