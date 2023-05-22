//
//  Interval.cpp
//  SINGER
//
//  Created by Yun Deng on 4/9/22.
//

#include "Interval.hpp"

Interval::Interval() {}

Interval::Interval(Branch b, float tl, float tu, int init_pos) {
    branch = b;
    lb = tl;
    ub = tu;
    assert(lb <= ub);
    // assert(ub >= 0.001);
    assert(b.lower_node->time <= tl and b.upper_node->time >= tu);
    start_pos = init_pos;
}

void Interval::assign_weight(float w) {
    weight = w;
}

void Interval::assign_time(float t) {
    assert(t >= lb and t <= ub);
    assert(!isinf(time));
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
    assert(!isinf(time));
    assert(time >= lb and time <= ub);
}

bool Interval::full(float t) {
    assert(lb >= t);
    return lb == max(t, branch.lower_node->time) and ub == branch.upper_node->time;
}

bool Interval::operator<(const Interval &other) const {
    if (start_pos != other.start_pos) {
        return start_pos < other.start_pos;
    }
    if (branch != other.branch) {
        return branch < other.branch;
    }
    if (ub != other.ub) {
        return ub < other.ub;
    }
    return lb < other.lb;
}

bool Interval::operator==(const Interval &other) const {
    if (start_pos != other.start_pos) {
        return false;
    }
    if (branch != other.branch) {
        return false;
    }
    if (ub != other.ub) {
        return false;
    }
    if (lb != other.lb) {
        return false;
    }
    return true;
}

bool Interval::operator!=(const Interval &other) const {
    if (start_pos != other.start_pos) {
        return true;
    }
    if (branch != other.branch) {
        return true;
    }
    if (ub != other.ub) {
        return true;
    }
    if (lb != other.lb) {
        return true;
    }
    return false;
}

shared_ptr<Interval> create_interval(Branch b, float tl, float tu, int init_pos) {
    return make_shared<Interval>(b, tl, tu, init_pos);
}

Interval_info::Interval_info() {
}

Interval_info::Interval_info(Branch b, float tl, float tu) {
    assert(tl <= tu);
    assert(tl >= b.lower_node->time and tu <= b.upper_node->time);
    branch = b;
    lb = tl;
    ub = tu;
}

bool Interval_info::operator==(const Interval_info& other) const {
    if (seed_pos != other.seed_pos) {
        return false;
    }
    if (branch != other.branch) {
        return false;
    }
    if (ub != other.ub) {
        return false;
    }
    if (lb != other.lb) {
        return false;
    }
    return true;
}

bool Interval_info::operator!=(const Interval_info& other) const {
    if (seed_pos != other.seed_pos) {
        return true;
    }
    if (branch != other.branch) {
        return true;
    }
    if (ub != other.ub) {
        return true;
    }
    if (lb != other.lb) {
        return true;
    }
    return false;
}

bool Interval_info::operator<(const Interval_info& other) const {
    if (seed_pos != other.seed_pos) {
        return seed_pos < other.seed_pos;
    }
    if (branch != other.branch) {
        return branch < other.branch;
    }
    if (ub != other.ub) {
        return ub < other.ub;
    }
    return lb < other.lb;
}
