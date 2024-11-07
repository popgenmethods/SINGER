//
//  Rate_map.hpp
//  SINGER
//
//  Created by Yun Deng on 6/4/24.
//

#ifndef Rate_map_hpp
#define Rate_map_hpp

#include <stdio.h>
#include <limits>
#include "Node.hpp"

class Rate_map {
    
public:
    
    // The rate map is for the entire chromosome, typically, but you may be running SINGER on only
    // a subset of that chromosome. This offset is where SINGER's concept of 0 starts.
    double offset = std::numeric_limits<double>::max();
    double sequence_length = std::numeric_limits<double>::max();
    vector<double> coordinates = {};
    vector<double> rate_distances = {};
    
    Rate_map(double poffset = 0.0);
    
    void load_map(string mut_map_file);
    
    double cumulative_distance(double x);
    
    double segment_distance(double x, double y);
    
    double mean_rate();
    
private:
    int find_index(double x);
};

#endif /* Rate_map_hpp */
