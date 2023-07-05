//
//  Distribution.hpp
//  SINGER
//
//  Created by Yun Deng on 7/2/23.
//

#ifndef Distribution_hpp
#define Distribution_hpp

#include <stdio.h>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

class Distribution {
    
public:
    
    map<float, float> quantile = {};
    map<float, float> cdf = {};
    
    Distribution();
    
    void read(string filename);
    
    float get_cdf(float x);
    
    float get_quantile(float q);
    
};

#endif /* Distribution_hpp */
