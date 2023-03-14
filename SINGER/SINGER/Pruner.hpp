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
    
    virtual void prune_arg(ARG a) = 0;
    
};

#endif /* Pruner_hpp */
