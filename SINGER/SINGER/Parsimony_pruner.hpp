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
#include "Pruner.hpp"

class Parsimony_pruner : public Pruner {
    
    int max_mismatch = 1;
    int curr_pos = 0;
    map<Branch, float> curr_mismatch = {};
    map<int, set<Branch>> reduced_set = {};

    Parsimony_pruner();
    
    void start_search(Node *n, set<float> mutations, set<Branch> branches);

    void mutation_update(set<float> mutations, Node *n);

    void recombination_update(Recombination &r);

    // private:
    
    float count_mismatch(Branch branch, set<float> mutations, Node *n);
    
};

#endif /* Parsimony_pruner_hpp */
