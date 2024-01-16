//
//  Reconstruction.hpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#ifndef Reconstruction_hpp
#define Reconstruction_hpp

#include <stdio.h>
#include "Tree.hpp"

class Reconstruction {
    
    virtual void update(Recombination &r) = 0;
    
    virtual void reconstruct(double pos) = 0;
    
};

#endif /* Reconstruction_hpp */
