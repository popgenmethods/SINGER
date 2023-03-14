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
    
    Branch branch = Branch();
    Branch prev_branch = Branch();
    Branch next_branch = Branch();
    int mismatch_sum = 0;
    
};


#endif /* Branch_node_hpp */
