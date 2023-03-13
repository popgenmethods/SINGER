//
//  Branch.hpp
//  SINGER
//
//  Created by Yun Deng on 3/31/22.
//

#ifndef Branch_hpp
#define Branch_hpp

#include <stdio.h>
#include "Node.hpp"

class Branch {
    
public:
    
    Node *lower_node;
    Node *upper_node;
    
    Branch();
    
    Branch(Node *l, Node *u);
    
    float length();
};

bool operator==(const Branch& b, const Branch& c);

bool operator!=(const Branch& b, const Branch& c);

bool operator<(const Branch& b, const Branch& c);

#endif /* Branch_hpp */
