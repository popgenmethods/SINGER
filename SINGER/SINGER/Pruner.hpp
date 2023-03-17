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

class Pruner {
    
public:
    
    virtual map<int, set<Branch>> prune_arg(ARG a, map<int, Node *> lower_nodes) = 0;
    
};

#endif /* Pruner_hpp */
