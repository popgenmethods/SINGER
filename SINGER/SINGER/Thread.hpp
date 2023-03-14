//
//  Thread.hpp
//  SINGER
//
//  Created by Yun Deng on 1/2/23.
//

#ifndef Thread_hpp
#define Thread_hpp

#include <stdio.h>
#include "ARG.hpp"

class Thread {
    
public:
    
    virtual void thread(ARG &a, Node *n);
    
    virtual void internal_rethread(ARG &a, tuple<int, Branch, float> cut_point);
    
    virtual void terminal_rethread(ARG &old_arg, tuple<int, Branch, float> cut_point);
    
};

#endif /* Thread_hpp */
