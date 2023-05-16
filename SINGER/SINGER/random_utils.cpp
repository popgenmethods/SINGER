//
//  random_utils.cpp
//  SINGER
//
//  Created by Yun Deng on 5/10/23.
//

#include "random_utils.hpp"

std::mt19937 random_engine;
std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);

float uniform_random() {
    return uniform_distribution(random_engine);
    // return (float) std::rand()/RAND_MAX;
}

void set_seed(unsigned seed) {  // Implement the set_seed function
    random_engine.seed(seed);
}
