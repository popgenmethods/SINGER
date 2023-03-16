//
//  Mutation_matcher.hpp
//  SINGER
//
//  Created by Yun Deng on 3/16/23.
//

#ifndef Mutation_matcher_hpp
#define Mutation_matcher_hpp

#include <stdio.h>
#include "ARG.hpp"

class Mutation_matcher {
    
    map<float, set<Branch>> derived_branches = {};
    map<float, set<Branch>> ancestral_branches = {};
    
    Mutation_matcher(ARG &a);
    
    set<Branch> get_matches(float x, float state);
    
    float count_match(float x, float state);
    
    // private:
    
    void classify_branches(Tree tree, float x);
    
};

#endif /* Mutation_matcher_hpp */
