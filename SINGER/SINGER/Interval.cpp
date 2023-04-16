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
