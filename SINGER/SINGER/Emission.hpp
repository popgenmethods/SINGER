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
#include "ARG.hpp"

class Emission {
    
public:
    
    virtual double null_emit(Branch &branch, double time, double theta, Node_ptr node) = 0;
    virtual double mut_emit(Branch &branch, double time, double theta, double bin_size, set<double> &mut_set, Node_ptr node) = 0;
    virtual double emit(Branch &branch, double time, double theta, double bin_size, vector<double> &emissions, Node_ptr node) = 0;
};

#endif /* Emission_hpp */
