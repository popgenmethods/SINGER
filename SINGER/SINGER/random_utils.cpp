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

std::string get_time() {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
    auto timer = system_clock::to_time_t(now);
    std::tm bt = *std::localtime(&timer);
    std::ostringstream oss;
    oss << "[" << std::put_time(&bt, "%H:%M:%S"); // HH:MM:SS
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count() << "]";
    return oss.str();
}
