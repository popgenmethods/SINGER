//
//  Parsimony_pruner.hpp
//  SINGER
//
//  Created by Yun Deng on 3/14/23.
//

#ifndef Parsimony_pruner_hpp
#define Parsimony_pruner_hpp

#include <stdio.h>
#include "Pruner.hpp"
#include "ARG.hpp"

class Parsimony_pruner {
    
public:
    
    int max_mismatch = 0;
    map<float, Node *> nodes = {};
    map<Branch, float> curr_mismatch = {};
    map<float, float> match_map = {};
    
    map<Branch, set<Branch>> transitions = {};

    Parsimony_pruner(ARG &a, map<float, Node *> base_nodes);
    
    void start_search(ARG &a, Node *n, float m);
    
    void extend(ARG &a);

    void mutation_forward(Node *n, float m);

    void recombination_forward(Recombination &r);

    // private:
    
    Node *get_node_at(float x);
    
    void get_match_map(ARG &a, map<float, Node *> base_nodes);
    
    float find_minimum_match();
    
    float count_mismatch(Branch branch, Node *n, float m);
    
    void transition_helper(Branch sb, Branch tb);
    
    void update_mismatch();
    
    void extend_forward(ARG &a, Node *n, float x);
    
    void extend_backward(ARG &a, Node *n, float x);
    
};

#endif /* Parsimony_pruner_hpp */