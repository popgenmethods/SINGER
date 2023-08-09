//
//  Tree.cpp
//  SINGER
//
//  Created by Yun Deng on 4/12/22.
//

#include "Tree.hpp"

Tree::Tree() {
}

float Tree::length() {
    float l = 0;
    for (Branch b : branches) {
        if (b.upper_node->index != -1) {
            l += b.length();
        }
    }
    return l;
}

void Tree::delete_branch(const Branch &b) {
    assert(b.upper_node != nullptr and b.lower_node != nullptr);
    assert(branches.count(b) > 0);
    branches.erase(b);
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
    assert(branches.count(b) == 0);
    branches.insert(b);
    parents[b.lower_node] = b.upper_node;
    children[b.upper_node].insert(b.lower_node);
}

void Tree::internal_insert_branch(const Branch &b, float cut_time) {
    if (b.upper_node->time <= cut_time) {
        return;
    }
    assert(branches.count(b) == 0);
    branches.insert(b);
    parents[b.lower_node] = b.upper_node;
    children[b.upper_node].insert(b.lower_node);
}

void Tree::internal_delete_branch(const Branch &b, float cut_time) {
    if (b.upper_node->time <= cut_time) {
        return;
    }
    assert(branches.count(b) > 0);
    branches.erase(b);
    parents.erase(b.lower_node);
    unordered_set<Node_ptr> &children_nodes = children[b.upper_node];
    if (children_nodes.size() == 1) {
        children.erase(b.upper_node);
    } else {
        children_nodes.erase(b.lower_node);
    }
}

void Tree::forward_update(Recombination &r) {
    for (const Branch &b : r.deleted_branches) {
        delete_branch(b);
    }
    for (const Branch &b : r.inserted_branches) {
        insert_branch(b);
    }
}

void Tree::backward_update(Recombination &r) {
    for (const Branch &b : r.inserted_branches) {
        delete_branch(b);
    }
    for (const Branch &b : r.deleted_branches) {
        insert_branch(b);
    }
}

void Tree::remove(Branch b, Node_ptr n) {
    assert(branches.count(b) > 0);
    assert(b.upper_node->index >= 0);
    Branch joining_branch = find_joining_branch(b);
    Node_ptr sibling = find_sibling(b.lower_node);
    Node_ptr parent = parents[b.upper_node];
    Branch sibling_branch = Branch(sibling, b.upper_node);
    Branch parent_branch = Branch(b.upper_node, parent);
    Branch cut_branch = Branch(b.lower_node, n);
    delete_branch(b);
    assert(branches.count(sibling_branch) > 0);
    assert(branches.count(parent_branch) > 0);
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

/*
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
 */

Branch Tree::find_joining_branch(Branch removed_branch) {
    if (removed_branch == Branch()) {
        return Branch();
    }
    Node_ptr p = parents[removed_branch.upper_node];
    Node_ptr c = find_sibling(removed_branch.lower_node);
    assert(parents[c] == removed_branch.upper_node);
    return Branch(c, p);
}

/*
pair<Branch, float> Tree::sample_cut_point() {
    float p = random();
    float q = random();
    float l = length()*q;
    Branch branch;
    float time = 0.0f;
    for (Branch b : branches) {
        if (b.upper_node->index != -1) {
            l -= b.length();
        }
        if (l < 0) {
            branch = b;
            time = b.lower_node->time + p*(b.upper_node->time - b.lower_node->time);
            break;
        }
    }
    assert(branch != Branch());
    return {branch, time};
}
*/

pair<Branch, float> Tree::sample_cut_point() {
    float t_max = branches.rbegin()->lower_node->time;
    float cut_time = random()*t_max;
    vector<Branch> candidates = {};
    for (const Branch &b : branches) {
        if (b.upper_node->time > cut_time and b.lower_node->time <= cut_time) {
            candidates.push_back(b);
        }
    }
    int index = (int) floor(candidates.size()*uniform_random());
    return {candidates[index], cut_time};
}

void Tree::internal_cut(float cut_time) {
    for (auto it = branches.begin(); it != branches.end();) {
        if (it->upper_node->time <= cut_time) {
            parents.erase(it->lower_node);
            it = branches.erase(it);
        } else {
            ++it;
        }
    }
}

void Tree::internal_forward_update(Recombination &r, float cut_time) {
    for (const Branch &b : r.deleted_branches) {
        internal_delete_branch(b, cut_time);
    }
    for (const Branch &b : r.inserted_branches) {
        internal_insert_branch(b, cut_time);
    }
}

void Tree::internal_backward_update(Recombination &r, float cut_time) {
    for (const Branch &b : r.deleted_branches) {
        internal_insert_branch(b, cut_time);
    }
    for (const Branch &b : r.inserted_branches) {
        internal_delete_branch(b, cut_time);
    }
}

float Tree::prior_likelihood() {
    float log_likelihood = 0;
    set<float> coalescence_times = {};
    int num_leaves = (int) (branches.size() + 1)/2;
    for (Branch b : branches) {
        coalescence_times.insert(b.lower_node->time);
        coalescence_times.insert(b.upper_node->time);
    }
    vector<float> sorted_coalescence_times = vector(coalescence_times.begin(), coalescence_times.end());
    for (int i = 0; i < num_leaves - 1; i++) {
        float lambda = (num_leaves - i)*(num_leaves - i - 1);
        float t = sorted_coalescence_times[i+1] - sorted_coalescence_times[i];
        log_likelihood += log_exp(lambda, t);
    }
    return log_likelihood;
}

float Tree::data_likelihood(float theta, float pos) {
    float log_likelihood = 0;
    float branch_likelihood = 0;
    for (Branch b : branches) {
        if (b.length() != numeric_limits<float>::infinity()) {
            float sl = b.lower_node->get_state(pos);
            float su = b.upper_node->get_state(pos);
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

float Tree::null_likelihood(float theta) {
    return -theta*length();
}

float Tree::data_likelihood(float theta, float bin_size, set<float> mutations) {
    float log_likelihood = 0;
    for (float x : mutations) {
        log_likelihood += data_likelihood(theta/bin_size, x);
    }
    float prop = 1 - mutations.size()/bin_size;
    log_likelihood += null_likelihood(-theta*prop);
    return log_likelihood;
}

float Tree::transition_likelihood(Recombination &r) {
    float log_likelihood = 0;
    log_likelihood -= log(length());
    set<float> coalescence_times = {};
    for (Branch b : branches) {
        if (b.upper_node->time > r.start_time) {
            coalescence_times.insert(b.upper_node->time);
        }
    }
    vector<float> sorted_coalescence_times = vector(coalescence_times.begin(), coalescence_times.end());
    int num_leaves = (int) sorted_coalescence_times.size();
    float base_time = r.start_time;
    float join_time = r.inserted_node->time;
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

float Tree::log_exp(float lambda, float x) {
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

void Tree::impute_states(float m, set<Branch> &mutation_branches) {
    map<Node_ptr, float> states = {};
    for (const Branch &b : mutation_branches) {
        states[b.lower_node] = b.lower_node->get_state(m);
        states[b.upper_node] = b.upper_node->get_state(m);
    }
    for (const Branch &b : branches) {
        impute_states_helper(b.lower_node, states);
    }
    for (auto &x : states) {
        x.first->write_state(m, x.second);
    }
}

void Tree::impute_states_helper(Node_ptr n, map<Node_ptr, float> &states) {
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

float Tree::random() {
    float p = uniform_random();
    return p;
}
