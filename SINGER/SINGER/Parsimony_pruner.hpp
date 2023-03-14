//
//  Parsimony_pruner.hpp
//  SINGER
//
//  Created by Yun Deng on 3/14/23.
//

#ifndef Parsimony_pruner_hpp
#define Parsimony_pruner_hpp

#include <stdio.h>
#include "ARG.hpp"
#include "Branch_node.hpp"

class Parsimony_pruner {
    
    int max_mismatch = 0;
    map<int, set<Branch_node*, compare_branch_node>> prune_graph = {};
    set<Branch_node *, compare_branch_node> curr_nodes = {};

    Parsimony_pruner();

    void prune(ARG a, map<int, Node *> lower_nodes);

    void mutation_update(float x, Node *n);

    void mutation_update(set<float> mutations, Node *n);

    void recombination_update(Recombination &r);

    // private:

    void update_helper(Branch_node *bn, Branch b, int x);
    
};

#endif /* Parsimony_pruner_hpp */
