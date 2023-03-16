//
//  Mutation_matcher.cpp
//  SINGER
//
//  Created by Yun Deng on 3/16/23.
//

#include "Mutation_matcher.hpp"

Mutation_matcher::Mutation_matcher(ARG &a) {
    int start = 0;
    int end = a.bin_num;
    map<int, Recombination>::iterator recomb_it = a.recombination_info.upper_bound(0);
    map<int, set<float>>::iterator mut_it = a.mutation_info.lower_bound(0);
    Tree tree = Tree();
    for (int i = start; i < end; i++) {
        if (i == recomb_it->first) {
            Recombination r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
        if (i == mut_it->first) {
            set<float> mutations = mut_it->second;
            for (float x : mutations) {
                classify_branches(tree, x);
            }
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
        match_count = ancestral_branches.at(x).size();
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
}
