//
//  Distribution.cpp
//  SINGER
//
//  Created by Yun Deng on 7/2/23.
//

#include "Distribution.hpp"

Distribution::Distribution() {}

void Distribution::read(string filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }
    string line;
    while (getline(file, line)) {
        istringstream ss(line);
        float x, q;
        if (!(ss >> x >> q)) {
            cerr << "Failed to parse line: " << line << endl;
            continue;
        }
        cdf[x] = q;
        quantile[q] = x;
    }
    file.close();
}

float Distribution::get_cdf(float x) {
    auto u_it = cdf.upper_bound(x);
    auto l_it = prev(u_it);
    float q = l_it->second + (x - l_it->first)/(u_it->first - l_it->first)*(u_it->second - l_it->second);
    return q;
}

float Distribution::get_quantile(float q) {
    auto u_it = quantile.upper_bound(q);
    auto l_it = prev(u_it);
    float x = l_it->second + (q - l_it->first)/(u_it->first - l_it->first)*(u_it->second - l_it->second);
    return x;
}
