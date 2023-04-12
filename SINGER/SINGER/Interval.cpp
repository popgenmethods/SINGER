//
//  Interval.cpp
//  SINGER
//
//  Created by Yun Deng on 4/9/22.
//

#include "Interval.hpp"

Interval::Interval(Branch b, float tl, float tu, float init_pos) {
    branch = b;
    lb = tl;
    ub = tu;
    assert(lb <= ub);
    start_pos = init_pos;
}

void Interval::assign_weight(float w) {
    weight = w;
}

void Interval::assign_time(float t) {
    assert(t >= lb and t <= ub);
    assert(t != numeric_limits<float>::infinity());
    time = t;
}

void Interval::fill_time() {
    if (ub == numeric_limits<float>::infinity()) {
        time = lb + log(2);
    } else if (abs(lb - ub) < 1e-3) {
        time = 0.5*(lb + ub);
    } else {
        float lq = 1 - exp(-lb);
        float uq = 1 - exp(-ub);
        if (uq - lq < 1e-3) {
            time = 0.5*(lb + ub);
        } else {
            float q = 0.5*(lq + uq);
            time = -log(1 - q);
        }
    }
    assert(time != numeric_limits<float>::infinity());
    assert(time >= lb and time <= ub);
}

bool Interval::full(float t) {
    assert(lb >= t);
    if (lb == t) {
        return ub == branch.upper_node->time;
    } else {
        return lb == branch.lower_node->time and ub == branch.upper_node->time;
    }
}

void Interval::set_source(vector<Interval *> intervals, vector<float> weights) {
    source_intervals = intervals;
    source_weights = weights;
}

void Interval::set_node(Node *n) {
    assert(n != nullptr);
    node = n;
}

Interval* Interval::sample_source() {
    if (source_intervals.size() == 0) {
        return nullptr;
    }
    if (source_intervals.size() == 1) {
        return source_intervals.front();
    }
    float weight_sum = accumulate(source_weights.begin(), source_weights.end(), 0.0);
    float p = (float) rand()/RAND_MAX;
    weight_sum *= p;
    for (int i = 0; i < source_intervals.size(); i++) {
        weight_sum -= source_weights[i];
        if (weight_sum < 0) {
            return source_intervals[i];
        }
    }
    cerr << "sample source interval failed!" << endl;
    exit(1);
}

Interval_info::Interval_info() {
}

Interval_info::Interval_info(Branch b, float tl, float tu) {
    assert(tl <= tu);
    branch = b;
    lb = tl;
    ub = tu;
}

bool operator==(const Interval_info& i1, const Interval_info& i2) {
    if (i1.branch != i2.branch) {
        return false;
    }
    if (i1.ub != i2.ub) {
        return false;
    }
    if (i1.lb != i2.lb) {
        return false;
    }
    return true;
}

bool operator!=(const Interval_info& i1, const Interval_info& i2) {
    if (i1.branch != i2.branch) {
        return true;
    }
    if (i1.ub != i2.ub) {
        return true;
    }
    if (i1.lb != i2.lb) {
        return true;
    }
    return false;
}

bool operator<(const Interval_info& i1, const Interval_info& i2) {
    if (i1.branch != i2.branch) {
        return i1.branch < i2.branch;
    }
    if (i1.ub != i2.ub) {
        return i1.ub < i2.ub;
    }
    return i1.lb < i2.lb;
}

// private methods:
