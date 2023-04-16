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
    
    virtual float null_emit(Branch branch, float time, float theta, Node *node) = 0;
    virtual float mut_emit(Branch branch, float time, float theta, float bin_size, set<float> &mut_set, Node *node) = 0;
    virtual float emit(Branch branch, float time, float theta, float bin_size, vector<float> &emissions, Node *node) = 0;
    
};

#endif /* Emission_hpp */
