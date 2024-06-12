//
//  Rate_map.hpp
//  SINGER
//
//  Created by Yun Deng on 6/4/24.
//

#ifndef Rate_map_hpp
#define Rate_map_hpp

#include <stdio.h>
#include "Node.hpp"

class Rate_map {
    
public:
    
    double sequence_length = INT_MAX;
    vector<double> coordinates = {};
    vector<double> rate_distances = {};
    
    Rate_map();
    
    void load_map(string mut_map_file);
    
    int find_index(double x);
    
    double cumulative_distance(double x);
    
    double segment_distance(double x, double y);
    
    double mean_rate();
    
};

#endif /* Rate_map_hpp */
