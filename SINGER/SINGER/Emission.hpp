//
//  Emission.hpp
//  SINGER
//
//  Created by Yun Deng on 4/3/23.
//

#ifndef Emission_hpp
#define Emission_hpp

#include <stdio.h>
#include "Branch.hpp"

class Emission {
    
public:
    
    virtual ~Emission() = default;
    virtual float null_emit(Branch branch, float time, float theta, Node *node) const = 0;
    virtual float mut_emit(Branch branch, float time, float theta, float mut_pos, Node *node) const = 0;
    
};

#endif /* Emission_hpp */
