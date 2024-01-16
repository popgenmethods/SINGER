//
//  Tree.cpp
//  SINGER
//
//  Created by Yun Deng on 4/12/22.
//

#include "Tree.hpp"

Tree::Tree() {
}

double Tree::length() {
    double l = 0;
    for (auto &x : parents) {
        if (x.second->index != -1) {
            l += x.second->time - x.first->time;
        }
    }
    return l;
}

void Tree::delete_branch(const Branch &b) {
    assert(b.upper_node != nullptr and b.lower_node != nullptr);
    parents.erase(b.lower_node);
    unordered_set<Node_ptr> &children_nodes = children[b.upper_node];
    if (children_nodes.size() == 1) {
        children.erase(b.upper_node);
    } else {
        children_nodes.erase(b.lower_node);  
    }
}

void Tree::insert_branch(const Branch &b) {
    assert(b.upper_node != nullptr and b.lower_node != nullptr);
    parents[b.lower_node] = b.upper_node;
    children[b.upper_node].insert(b.lower_node);
}

void Tree::internal_insert_branch(const Branch &b, double cut_time) {
    if (b.upper_node->time <= cut_time) {
        return;
    }
    parents[b.lower_node] = b.upper_node;
    children[b.upper_node].insert(b.lower_node);
}

void Tree::internal_delete_branch(const Branch &b, double cut_time) {
    if (b.upper_node->time <= cut_time) {
        return;
    }
    parents.erase(b.lower_node);
    unordered_set<Node_ptr> &children_nodes = children[b.upper_node];
    if (children_nodes.size() == 1) {
        children.erase(b.upper_node);
    } else {
        children_nodes.erase(b.lower_node);
    }
}

void Tree::forward_update(Recombination &r) {
    int prev_size = (int) parents.size();
    for (const Branch &b : r.deleted_branches) {
        delete_branch(b);
    }
    for (const Branch &b : r.inserted_branches) {
        insert_branch(b);
    }
    int after_size = (int) parents.size();
    assert(prev_size == after_size or r.pos == 0);
}

void Tree::backward_update(Recombination &r) {
    int prev_size = (int) parents.size();
    for (const Branch &b : r.inserted_branches) {
        delete_branch(b);
    }
    for (const Branch &b : r.deleted_branches) {
        insert_branch(b);
    }
    int after_size = (int) parents.size();
    assert(prev_size == after_size or r.pos == 0);
}

void Tree::remove(Branch b, Node_ptr n) {
    // assert(branches.count(b) > 0);
    assert(b.upper_node->index >= 0);
    Branch joining_branch = find_joining_branch(b);
    Node_ptr sibling = find_sibling(b.lower_node);
    Node_ptr parent = parents[b.upper_node];
    Branch sibling_branch = Branch(sibling, b.upper_node);
    Branch parent_branch = Branch(b.upper_node, parent);
    Branch cut_branch = Branch(b.lower_node, n);
    delete_branch(b);
    // assert(branches.count(sibling_branch) > 0);
    // assert(branches.count(parent_branch) > 0);
    delete_branch(sibling_branch);
    delete_branch(parent_branch);
    insert_branch(joining_branch);
    insert_branch(cut_branch);
}

void Tree::add(Branch added_branch, Branch joining_branch, Node_ptr n) {
    Branch lower_branch = Branch(joining_branch.lower_node, added_branch.upper_node);
    Branch upper_branch = Branch(added_branch.upper_node, joining_branch.upper_node);
    delete_branch(joining_branch);
    insert_branch(lower_branch);
    insert_branch(upper_branch);
    insert_branch(added_branch);
    if (n != nullptr) {
        Branch cut_branch = Branch(added_branch.lower_node, n);
        delete_branch(cut_branch);
    }
}

/*
Node_ptr Tree::find_sibling(Node_ptr n) {
    Node_ptr p = parents[n];
    Branch b = Branch(n, p);
    set<Branch>::iterator branch_it = branches.find(b);
    branch_it++;
    Branch candidate = *branch_it;
    if (candidate.upper_node != p) {
        branch_it--;
        branch_it--;
    }
    candidate = *branch_it;
    Node_ptr s = (*branch_it).lower_node;
    return s;
}
 */

Node_ptr Tree::find_sibling(Node_ptr n) {
    Node_ptr p = parents[n];
    unordered_set<Node_ptr> &candidates = children[p];
    Node_ptr c;
    auto c_it = candidates.begin();
    if (*c_it != n) {
        c = *c_it;
    } else {
        c = *(next(c_it));
    }
    return c;
}

Branch Tree::find_joining_branch(Branch removed_branch) {
    if (removed_branch == Branch()) {
        return Branch();
    }
    Node_ptr p = parents[removed_branch.upper_node];
    Node_ptr c = find_sibling(removed_branch.lower_node);
    assert(parents[c] == removed_branch.upper_node);
    return Branch(c, p);
}

pair<Branch, double> Tree::sample_cut_point() {
    double root_time = parents.rbegin()->first->time;
    double cut_time = random()*root_time;
    vector<Branch> candidates = {};
    for (auto &x : parents) {
        if (x.second->time > cut_time and x.first->time <= cut_time) {
            candidates.push_back(Branch(x.first, x.second));
        }
    }
    int index = (int) floor(candidates.size()*uniform_random());
    index = min((int) candidates.size() - 1, index);
    return {candidates[index], cut_time};
}

void Tree::internal_cut(double cut_time) {
    for (auto it = parents.begin(); it != parents.end();) {
        if (it->second->time <= cut_time) {
            it = parents.erase(it);
        } else {
            ++it;
        }
    }
}

void Tree::internal_forward_update(Recombination &r, double cut_time) {
    for (const Branch &b : r.deleted_branches) {
        internal_delete_branch(b, cut_time);
    }
    for (const Branch &b : r.inserted_branches) {
        internal_insert_branch(b, cut_time);
    }
}

void Tree::internal_backward_update(Recombination &r, double cut_time) {
    for (const Branch &b : r.inserted_branches) {
        internal_delete_branch(b, cut_time);
    }
    for (const Branch &b : r.deleted_branches) {
        internal_insert_branch(b, cut_time);
    }
}

double Tree::prior_likelihood() {
    double log_likelihood = 0;
    set<double> coalescence_times = {};
    int num_leaves = (int) (parents.size() + 1)/2;
    for (auto &x : parents) {
        coalescence_times.insert(x.first->time);
        coalescence_times.insert(x.second->time);
    }
    vector<double> sorted_coalescence_times = vector(coalescence_times.begin(), coalescence_times.end());
    for (int i = 0; i < num_leaves - 1; i++) {
        double lambda = (num_leaves - i)*(num_leaves - i - 1);
        double t = sorted_coalescence_times[i+1] - sorted_coalescence_times[i];
        log_likelihood += log_exp(lambda, t);
    }
    return log_likelihood;
}

double Tree::data_likelihood(double theta, double pos) {
    double log_likelihood = 0;
    double branch_likelihood = 0;
    for (auto &x : parents) {
        Branch b = Branch(x.first, x.second);
        if (b.length() != numeric_limits<double>::infinity()) {
            double sl = b.lower_node->get_state(pos);
            double su = b.upper_node->get_state(pos);
            if (sl != su) {
                branch_likelihood = log(theta) + log(b.length()) - theta*b.length();
            } else {
                branch_likelihood = -theta*b.length();
            }
            log_likelihood += branch_likelihood;
        }
    }
    return log_likelihood;
}

double Tree::null_likelihood(double theta) {
    return -theta*length();
}

double Tree::data_likelihood(double theta, double bin_size, set<double> mutations) {
    double log_likelihood = 0;
    for (double x : mutations) {
        log_likelihood += data_likelihood(theta/bin_size, x);
    }
    double prop = 1 - mutations.size()/bin_size;
    log_likelihood += null_likelihood(-theta*prop);
    return log_likelihood;
}

double Tree::transition_likelihood(Recombination &r) {
    double log_likelihood = 0;
    log_likelihood -= log(length());
    set<double> coalescence_times = {};
    for (auto &x : parents) {
        if (x.second->time > r.start_time) {
            coalescence_times.insert(x.second->time);
        }
    }
    vector<double> sorted_coalescence_times = vector(coalescence_times.begin(), coalescence_times.end());
    int num_leaves = (int) sorted_coalescence_times.size();
    double base_time = r.start_time;
    double join_time = r.inserted_node->time;
    for (int i = 0; i < sorted_coalescence_times.size(); i++) {
        if (sorted_coalescence_times[i] > join_time) {
            log_likelihood += log_exp(num_leaves, join_time - base_time);
            break;
        } else {
            log_likelihood += log_exp(num_leaves, sorted_coalescence_times[i] - base_time);
            base_time = sorted_coalescence_times[i];
            num_leaves -= 1;
        }
    }
    return log_likelihood;
}

// private methods:

double Tree::log_exp(double lambda, double x) {
    return -lambda*x + log(lambda);
}

int Tree::depth(Node_ptr n) {
    int depth = 0;
    while (!isinf(n->time)) {
        depth += 1;
        n = parents[n];
    }
    return depth;
}

Node_ptr Tree::LCA(Node_ptr n1, Node_ptr n2) {
    set<Node_ptr > ancestors = {};
    while (!isinf(n1->time)) {
        ancestors.insert(n1);
        n1 = parents[n1];
    }
    while (!isinf(n2->time)) {
        if (ancestors.count(n2) > 0) {
            return n2;
        }
        n2 = parents[n2];
    }
    return n1;
}

int Tree::distance(Node_ptr n1, Node_ptr n2) {
    if (n1 == n2) {
        return 0;
    }
    Node_ptr lca = LCA(n1, n2);
    int depth1 = depth(n1) - depth(lca);
    int depth2 = depth(n2) - depth(lca);
    return depth1 + depth2;
}

void Tree::impute_states(double m, set<Branch> &mutation_branches) {
    map<Node_ptr, double> states = {};
    for (const Branch &b : mutation_branches) {
        states[b.lower_node] = b.lower_node->get_state(m);
        states[b.upper_node] = b.upper_node->get_state(m);
    }
    for (auto &x : parents) {
        impute_states_helper(x.first, states);
    }
    for (auto &x : states) {
        x.first->write_state(m, x.second);
    }
}

void Tree::impute_states_helper(Node_ptr n, map<Node_ptr, double> &states) {
    if (n->index == -1) {
        states[n] = 0;
        return;
    }
    if (states.count(n) > 0) {
        return;
    }
    Node_ptr p = parents[n];
    impute_states_helper(p, states);
    states[n] = states[p];
}

double Tree::random() {
    double p = uniform_random();
    return p;
}
