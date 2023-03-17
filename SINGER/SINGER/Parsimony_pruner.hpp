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

class Parsimony_pruner : public Pruner {
    
    int max_mismatch = 1;
    map<Branch, float> curr_mismatch = {};
    
    map<Branch, set<Branch>> transitions = {};

    Parsimony_pruner();
    
    void start_search(Node *n, float m, set<Branch> branches);

    void mutation_forward(Node *n, float m);

    void recombination_forward(Recombination &r);

    // private:
    
    float count_mismatch(Branch branch, Node *n, float m);
    
    void transition_helper(Branch sb, Branch tb);
    
    void update_mismatch();
    
};

#endif /* Parsimony_pruner_hpp */
