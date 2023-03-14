//
//  Pruner.hpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#ifndef Pruner_hpp
#define Pruner_hpp

#include <stdio.h>
#include "ARG.hpp"
#include "Branch_node.hpp"

class Pruner {
    
public:
    
    int max_mismatch = 0;
    set<Branch_node> prune_graph = {};
    
    Pruner();
    
    void prune();
    
};

#endif /* Pruner_hpp */
