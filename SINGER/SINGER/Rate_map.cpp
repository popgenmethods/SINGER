//
//  Rate_map.cpp
//  SINGER
//
//  Created by Yun Deng on 6/4/24.
//

#include "Rate_map.hpp"

// Even if there are regions of truly no recombination, we have to pretend there is
// a really low rate to make the math work out currently.
static constexpr double MIN_RECOMB_RATE = 1e-16;

Rate_map::Rate_map(double poffset)
    : offset(poffset) {}

void Rate_map::load_map(string mut_map_file) {
    ifstream fin(mut_map_file);
    if (!fin.good()) {
        cerr << "input rate map file not found" << endl;
        exit(1);
    }
    rate_distances.push_back(0);
    double left;
    double right;
    double rate;
    double mut_dist;
    while (fin >> left >> right >> rate) {
        if (rate < MIN_RECOMB_RATE) {
            rate = MIN_RECOMB_RATE;
        }
        coordinates.push_back(left);
        mut_dist = rate_distances.back() + rate*(right - left);
        assert(mut_dist != rate_distances.back());
        rate_distances.push_back(mut_dist);
    }
    sequence_length = right;
    coordinates.push_back(sequence_length);
}

int Rate_map::find_index(double x) {
    auto it = upper_bound(coordinates.begin(), coordinates.end(), x);
    it--;
    int index = (int) distance(coordinates.begin(), it);
    assert(index >= 0 and index <= coordinates.size() - 1);
    return index;
}

double Rate_map::cumulative_distance(double x) {
    const double actualCoordinate = x + offset;
    int index = find_index(actualCoordinate);
    double prev_dist = rate_distances[index];
    double next_dist = rate_distances[index+1];
    double p = (actualCoordinate - coordinates[index])/(coordinates[index+1] - coordinates[index]);
    double dist = (1-p)*prev_dist + p*next_dist;
    return dist;
}

double Rate_map::segment_distance(double x, double y) {
    double cd = cumulative_distance(y) - cumulative_distance(x);
    assert(cd > 0.0);
    return cd;
}

double Rate_map::mean_rate() {
    double mr = rate_distances.back()/sequence_length;
    return mr;
}
