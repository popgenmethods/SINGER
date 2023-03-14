//
//  Branch_node.hpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#ifndef Branch_node_hpp
#define Branch_node_hpp

#include <stdio.h>
#include "ARG.hpp"

class Branch_node {
    
public:
    
    Branch branch = Branch();
    set<Branch_node *> prev_nodes = {};
    int start_pos = 0;
    int curr_mismatch = 0;
    
    Branch_node(Branch b, int x);
    
    void add_prev_node(Branch_node *bn);
    
    void mutation_update(float x, Node *n);
};

bool operator<(const Branch_node& b, const Branch_node& c);

struct compare_branch_node {
    
    bool operator() (const Branch_node *n1, const Branch_node *n2) const {
        return n1->branch < n2->branch;
    }
};

#endif /* Branch_node_hpp */
