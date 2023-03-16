//
//  BSP.hpp
//  SINGER
//
//  Created by Yun Deng on 3/15/23.
//

#ifndef BSP_hpp
#define BSP_hpp

#include <stdio.h>
#include "ARG.hpp"

class BSP {
    
    virtual map<int, pair<Branch, Node *>> sample_joining_branches() = 0;
    
};

#endif /* BSP_hpp */
