//
//  Threader.hpp
//  SINGER
//
//  Created by Yun Deng on 1/2/23.
//

#ifndef Threader_hpp
#define Threader_hpp

#include <stdio.h>
#include "ARG.hpp"

class Threader {
    
public:
    
    virtual void thread(ARG &a, Node *n) = 0;
    
    virtual void internal_rethread(ARG &a, tuple<int, Branch, float> cut_point) = 0;
    
    virtual void terminal_rethread(ARG &old_arg, tuple<int, Branch, float> cut_point) = 0;
    
};

#endif /* Threader_hpp */
