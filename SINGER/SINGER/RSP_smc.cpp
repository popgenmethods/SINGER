//
//  RSP_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#include "RSP_smc.hpp"

RSP_smc::RSP_smc() {
}

double RSP_smc::sample_start_time(Branch b, int density, double join_time, double cut_time) {
    double lb = b.lower_node->time;
    double ub = b.upper_node->time;
    lb = max(cut_time, lb);
    ub = min(join_time, ub);
    assert(lb < ub);
    double p;
    double t;
    double w;
    double start_time = 0;
    vector<double> start_times = {};
    vector<double> weights = {};
    double weight_sum = 0;
    for (int i = 0; i < density; i++) {
        t = random_time(lb, ub);
        w = recomb_pdf(t, ub);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
    }
    p = uniform_random();
    weight_sum = weight_sum*p;
    for (int i = 0; i < density; i++) {
        weight_sum -= weights[i];
        if (weight_sum <= 0) {
            start_time = start_times[i];
            break;
        }
    }
    assert(start_time >= cut_time);
    assert(start_time <= join_time);
    assert(start_time <= ub);
    return start_time;
}

pair<Branch, double> RSP_smc::sample_start_time(Branch b1, Branch b2, int density, double join_time, double cut_time) {
    double lb1 = max(b1.lower_node->time, cut_time);
    double ub1 = min(b1.upper_node->time, join_time);
    double lb2 = max(b2.lower_node->time, cut_time);
    double ub2 = min(b2.upper_node->time, join_time);
    assert(lb1 < ub1);
    assert(lb2 < ub2);
    double q = (ub1 - lb1)/(ub1 + ub2 - lb1 - lb2);
    int n1 = round(density*q);
    int n2 = density - n1;
    double p;
    double t;
    double w;
    Branch source_branch;
    double start_time = 0;
    vector<double> start_times(0);
    vector<double> weights(0);
    vector<int> branch_indices(0);
    double weight_sum = 0;
    for (int i = 0; i < n1; i++) {
        t = random_time(lb1, ub1);
        w = recomb_pdf(t, ub1);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
        branch_indices.push_back(1);
    }
    for (int i = 0; i < n2; i++) {
        t = random_time(lb2, ub2);
        w = recomb_pdf(t, ub2);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
        branch_indices.push_back(2);
    }
    p = uniform_random();
    weight_sum = weight_sum*p;
    for (int i = 0; i < density; i++) {
        weight_sum -= weights[i];
        if (weight_sum <= 0) {
            start_time = start_times[i];
            if (branch_indices[i] == 1) {
                source_branch = b1;
            } else {
                source_branch = b2;
            }
            break;
        }
    }
    assert(weight_sum <= 0);
    assert(start_time >= cut_time);
    assert(start_time <= join_time);
    assert(start_time <= min(ub1, ub2));
    return {source_branch, start_time};
}

void RSP_smc::sample_recombination(Recombination &r, double cut_time, Tree &tree) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    get_coalescence_rate(tree, r, cut_time);
    vector<Branch> source_candidates;
    for (Branch b : r.deleted_branches) {
        if (b.upper_node == r.deleted_node and b.lower_node->time < r.inserted_node->time) {
            Branch candidate_recombined_branch = Branch(b.lower_node, r.inserted_node);
            if (r.create(candidate_recombined_branch)) {
                source_candidates.push_back(b);
            }
        }
    }
    if (source_candidates.size() == 1) {
        r.source_branch = source_candidates[0];
        r.start_time = sample_start_time(r.source_branch, 20, r.inserted_node->time, cut_time);
    } else if (source_candidates.size() == 2) {
        pair<Branch, double> breakpoint = sample_start_time(source_candidates[0], source_candidates[1], 40, r.inserted_node->time, cut_time);
        r.source_branch = breakpoint.first;
        r.start_time = breakpoint.second;
    } else {
        cout << r.pos << " " << source_candidates.size() << endl;
        cerr << "no candidates in smc sampling" << endl;
        exit(1);
    }
    r.find_target_branch();
    r.find_recomb_info();
    assert(r.target_branch != Branch());
    assert(r.merging_branch != Branch());
    assert(r.start_time <= r.inserted_node->time);
    assert(r.start_time >= cut_time);
}

void RSP_smc::approx_sample_recombination(Recombination &r, double cut_time) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    vector<Branch> source_candidates;
    for (Branch b : r.deleted_branches) {
        if (b.upper_node == r.deleted_node and b.lower_node->time < r.inserted_node->time) {
            Branch candidate_recombined_branch = Branch(b.lower_node, r.inserted_node);
            if (r.create(candidate_recombined_branch)) {
                source_candidates.push_back(b);
            }
        }
    }
    double lb1, lb2, ub1, ub2;
    if (source_candidates.size() == 1) {
        r.source_branch = source_candidates[0];
        lb1 = max(cut_time, r.source_branch.lower_node->time);
        ub1 = min(r.source_branch.upper_node->time, r.inserted_node->time);
        // r.start_time = random_time(lb1, ub1, 0.5);
        r.start_time = choose_time(lb1, ub1);
    } else if (source_candidates.size() == 2) {
        lb1 = max(cut_time, source_candidates[0].lower_node->time);
        lb2 = max(cut_time, source_candidates[1].lower_node->time);
        ub1 = min(source_candidates[0].upper_node->time, r.inserted_node->time);
        ub2 = min(source_candidates[1].upper_node->time, r.inserted_node->time);
        double q = (ub1 - lb1)/(ub1 + ub2 - lb1 - lb2);
        double p = 0.5;
        if (p <= q) {
            r.source_branch = source_candidates[0];
            // r.start_time = random_time(lb1, ub1, 0.5);
            r.start_time = choose_time(lb1, ub1);
        } else {
            r.source_branch = source_candidates[1];
            // r.start_time = random_time(lb2, ub2, 0.5);
            r.start_time = choose_time(lb2, ub2);
        }
    } else {
        cout << r.pos << " " << source_candidates.size() << endl;
        cerr << "no candidates in smc sampling" << endl;
        exit(1);
    }
    if (r.deleted_node->time == r.inserted_node->time) {
        r.inserted_node->time = nextafter(r.inserted_node->time, numeric_limits<double>::infinity());
    }
    r.find_target_branch();
    r.find_recomb_info();
    assert(r.target_branch != Branch());
    assert(r.merging_branch != Branch());
    assert(r.start_time >= cut_time);
    double ub = min(r.deleted_node->time, r.inserted_node->time);
    if (r.start_time >= ub) {
        r.start_time = nextafter(ub, -numeric_limits<double>::infinity());
    }
    assert(r.start_time <= ub);
}

void RSP_smc::adjust(Recombination &r, double cut_time) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    double lb, ub;
    lb = max(cut_time, r.source_branch.lower_node->time);
    ub = min(r.deleted_node->time, r.inserted_node->time);
    r.start_time = choose_time(lb, ub);
    if (r.start_time >= ub or r.start_time <= lb) {
        r.start_time = 0.5*(lb + ub);
    }
    assert(r.start_time <= ub);
}

void RSP_smc::approx_sample_recombination(Recombination &r, double cut_time, double n) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    vector<Branch> source_candidates;
    for (Branch b : r.deleted_branches) {
        if (b.upper_node == r.deleted_node and b.lower_node->time < r.inserted_node->time) {
            Branch candidate_recombined_branch = Branch(b.lower_node, r.inserted_node);
            if (r.create(candidate_recombined_branch)) {
                source_candidates.push_back(b);
            }
        }
    }
    double lb1, lb2, ub1, ub2;
    if (source_candidates.size() == 1) {
        r.source_branch = source_candidates[0];
        lb1 = max(cut_time, r.source_branch.lower_node->time);
        ub1 = min(r.source_branch.upper_node->time, r.inserted_node->time);
        r.start_time = choose_time(lb1, ub1, n);
    } else if (source_candidates.size() == 2) {
        lb1 = max(cut_time, source_candidates[0].lower_node->time);
        lb2 = max(cut_time, source_candidates[1].lower_node->time);
        ub1 = min(source_candidates[0].upper_node->time, r.inserted_node->time);
        ub2 = min(source_candidates[1].upper_node->time, r.inserted_node->time);
        double q = (ub1 - lb1)/(ub1 + ub2 - lb1 - lb2);
        double p = 0.5;
        if (p <= q) {
            r.source_branch = source_candidates[0];
            r.start_time = choose_time(lb1, ub1, n);
        } else {
            r.source_branch = source_candidates[1];
            r.start_time = choose_time(lb2, ub2, n);
        }
    } else {
        cout << r.pos << " " << source_candidates.size() << endl;
        cerr << "no candidates in smc sampling" << endl;
        exit(1);
    }
    if (r.deleted_node->time == r.inserted_node->time) {
        r.inserted_node->time = nextafter(r.inserted_node->time, numeric_limits<double>::infinity());
    }
    r.find_target_branch();
    r.find_recomb_info();
    assert(r.target_branch != Branch());
    assert(r.merging_branch != Branch());
    assert(r.start_time >= cut_time);
    double ub = min(r.deleted_node->time, r.inserted_node->time);
    if (r.start_time >= ub) {
        r.start_time = nextafter(ub, -numeric_limits<double>::infinity());
    }
    assert(r.start_time <= ub);
}

void RSP_smc::adjust(Recombination &r, double cut_time, double n) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    double lb, ub;
    lb = max(cut_time, r.source_branch.lower_node->time);
    ub = min(r.deleted_node->time, r.inserted_node->time);
    r.start_time = choose_time(lb, ub, n);
    if (r.start_time >= ub or r.start_time <= lb) {
        r.start_time = 0.5*(lb + ub);
    }
    assert(r.start_time <= ub);
}

// private methods:

double RSP_smc::recomb_pdf(double s, double t) {
    double pdf = 1.0;
    double curr_time = s;
    double next_coalescence_time = s;
    map<double, int>::iterator rate_it = coalescence_rates.upper_bound(s);
    rate_it--;
    int rate;
    while (next_coalescence_time < t) {
        rate = rate_it->second;
        rate_it++;
        next_coalescence_time = rate_it->first;
        pdf *= exp(-rate*(-curr_time + min(t, next_coalescence_time)));
        curr_time = next_coalescence_time;
    }
    return pdf;
}

void RSP_smc::get_coalescence_rate(Tree &tree, Recombination &r, double cut_time) {
    coalescence_rates.clear();
    vector<double> coalescence_times = {cut_time};
    for (auto &x : tree.parents) {
        if (x.first->time > cut_time and x.second != r.deleted_node) {
            coalescence_times.push_back(x.first->time);
        }
    }
    coalescence_times.push_back(numeric_limits<double>::infinity());
    sort(coalescence_times.begin(), coalescence_times.end());
    int n = (int) coalescence_times.size();
    for (int i = 0; i < n; i++) {
        coalescence_rates[coalescence_times[i]] = n-i-1;
    }
}

double RSP_smc::random_time(double lb, double ub) {
    double t = (uniform_random()*0.01 + 0.99)*(ub - lb) + lb;
    return t;
}

double RSP_smc::random_time(double lb, double ub, double q) {
    double t = (q*uniform_random() + (1 - q))*(ub - lb) + lb;
    return t;
}

double RSP_smc::choose_time(double lb, double ub) {
    assert(!isinf(ub));
    double l = exp(lb);
    double u = exp(ub);
    double m = 0.5*(l + u);
    double mt = log(m);
    if (ub - lb < 0.01) {
        mt = 0.5*(lb + ub);
    }
    assert(mt >= lb and mt <= ub);
    return mt;
}

double RSP_smc::choose_time(double lb, double ub, double n) {
    double lambda = 0;
    double lambda_l = 0;
    double lambda_u = 0;
    double mt = 0;
    if (n > 1) {
        lambda_l = n/(n + (1 - n)*exp(-0.5*lb));
        lambda_u = n/(n + (1 - n)*exp(-0.5*ub));
        lambda = lambda_l;
    } else {
        lambda = 1;
    }
    double l = exp(lambda*lb - lambda*ub);
    double u = 1;
    double m = 0.5*(l + u);
    mt = ub + log(m)/lambda;
    if (ub - lb < 0.01) {
        mt = 0.5*(lb + ub);
    }
    assert(mt >= lb and mt <= ub);
    return mt;
}
