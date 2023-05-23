//
//  ARG.cpp
//  SINGER
//
//  Created by Yun Deng on 4/14/22.
//

#include "ARG.hpp"

ARG::ARG() {}

ARG::ARG(float N, float l) {
    root->set_index(-1);
    Ne = N;
    sequence_length = l;
    end = l;
    Recombination r = Recombination({}, {});
    r.set_pos(0.0);
    recombinations[0] = r;
    r = Recombination({}, {});
    r.set_pos(INT_MAX);
    recombinations[INT_MAX] = r;
    mutation_sites.insert(INT_MAX);
    mutation_branches[INT_MAX] = {};
}

ARG::~ARG() {
}

void ARG::discretize(float s) {
    auto recomb_it = recombinations.upper_bound(0);
    float curr_pos = 0;
    while (curr_pos < sequence_length) {
        // coordinates.push_back(max(curr_pos - 0.1f, 0.0f));
        coordinates.push_back(curr_pos);
        if (recomb_it->first < curr_pos + s) {
            curr_pos = recomb_it->first;
            recomb_it++;
        } else {
            curr_pos = min(curr_pos + s, sequence_length);
        }
    }
    coordinates.push_back(sequence_length);
    bin_num = (int) coordinates.size() - 1;
}

int ARG::get_index(float x) {
    /*
    int left = 0;
    int right = (int) coordinates.size() - 1;
    int mid;
    while (left <= right) {
        mid = left + (right - left) / 2;
        if (coordinates[mid] == x) {
            return mid;
        }
        if (coordinates[mid] < x) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return right;
     */
    auto it = upper_bound(coordinates.begin(), coordinates.end(), x);
    --it;
    int index = (int) distance(coordinates.begin(), it);
    return index;
}

void ARG::compute_rhos_thetas(float r, float m) {
    int n = (int) coordinates.size() - 1;
    for (int i = 0; i < n; i++) {
        rhos.push_back(r*(coordinates[i+1] - coordinates[i]));
        thetas.push_back(m*(coordinates[i+1] - coordinates[i]));
    }
}

void ARG::build_singleton_arg(Node_ptr n) {
    add_sample(n);
    Branch branch = Branch(n, root);
    Recombination r = Recombination({}, {branch});
    r.set_pos(0.0);
    recombinations[0] = r;
    for (float x : mutation_sites) {
        mutation_branches[x] = {branch};
    }
}

void ARG::add_sample(Node_ptr n) {
    sample_nodes.insert(n);
    for (float x : n->mutation_sites) {
        mutation_sites.insert(x);
    }
    removed_branches.clear();
    removed_branches[0] = Branch(n, root);
    removed_branches[sequence_length] = Branch();
    start_tree = get_tree_at(0);
    cut_pos = 0;
    start = 0;
    end = sequence_length;
}

void ARG::add_node(Node_ptr n) {
    if (n != root and n != nullptr) {
        node_set.insert(n);
    }
}

void ARG::add_new_node(float t) {
    if (!isinf(t)) {
        Node_ptr n = new_node(t);
        node_set.insert(n);
    }
}
 
Tree ARG::get_tree_at(float x) {
    Tree tree = Tree();
    auto recomb_it = recombinations.begin();
    while (recomb_it->first <= x) {
        Recombination &r = recomb_it->second;
        tree.forward_update(r);
        recomb_it++;
    }
    return tree;
}

Node_ptr ARG::get_query_node_at(float x) {
    auto query_it = removed_branches.upper_bound(x);
    query_it--;
    return query_it->second.lower_node;
}

Tree ARG::modify_tree_to(float x, Tree &reference_tree, float x0) {
    Tree tree = reference_tree;
    if (x == x0) {
        return tree;
    } else if (x > x0) {
        auto recomb_it = recombinations.upper_bound(x0);
        while (recomb_it->first <= x) {
            Recombination &r = recomb_it->second;
            tree.forward_update(r);
            ++recomb_it;
        }
        return tree;
    } else {
        auto recomb_it = recombinations.upper_bound(x0);
        --recomb_it;
        while (recomb_it->first > x) {
            Recombination &r = recomb_it->second;
            tree.backward_update(r);
            --recomb_it;
        }
        return tree;
    }
}

/*
void ARG::remove(tuple<float, Branch, float> cut_point) {
    float pos;
    Branch center_branch;
    float t;
    tie(pos, center_branch, t) = cut_point;
    cut_time = t;
    cut_node = new Node(cut_time);
    cut_node->set_index(-2);
    add_node(cut_node);
    Tree forward_tree = cut_tree;
    Tree backward_tree = cut_tree;
    auto f_it = recombinations.upper_bound(pos);
    auto b_it = recombinations.upper_bound(pos);
    Branch prev_joining_branch;
    Branch next_joining_branch;
    Branch prev_removed_branch = center_branch;
    Branch next_removed_branch = center_branch;
    while (next_removed_branch != Branch()) {
        Recombination &r = f_it->second;
        prev_joining_branch = forward_tree.find_joining_branch(prev_removed_branch);
        forward_tree.forward_update(r);
        next_removed_branch = r.trace_forward(t, prev_removed_branch);
        if (next_removed_branch.upper_node == root) {
            next_removed_branch = Branch();
        }
        next_joining_branch = forward_tree.find_joining_branch(next_removed_branch);
        r.remove(prev_removed_branch, next_removed_branch, prev_joining_branch, next_joining_branch, cut_node);
        removed_branches[min(r.pos, sequence_length)] = next_removed_branch;
        joining_branches[min(r.pos, sequence_length)] = next_joining_branch;
        f_it++;
        prev_removed_branch = next_removed_branch;
    }
    next_removed_branch = center_branch;
    prev_removed_branch = center_branch;
    while (prev_removed_branch != Branch()) {
        b_it--;
        Recombination &r = b_it->second;
        removed_branches[r.pos] = prev_removed_branch;
        next_joining_branch = backward_tree.find_joining_branch(next_removed_branch);
        joining_branches[r.pos] = next_joining_branch;
        backward_tree.backward_update(r);
        prev_removed_branch = r.trace_backward(t, next_removed_branch);
        if (prev_removed_branch.upper_node == root) {
            prev_removed_branch = Branch();
        }
        prev_joining_branch = backward_tree.find_joining_branch(prev_removed_branch);
        r.remove(prev_removed_branch, next_removed_branch, prev_joining_branch, next_joining_branch, cut_node);
        next_removed_branch = prev_removed_branch;
    }
    remove_empty_recombinations();
    remap_mutations();
    start = removed_branches.begin()->first;
    end = removed_branches.rbegin()->first;
    cut_tree.remove(center_branch, cut_node);
    start_tree = modify_tree_to(start, cut_tree, cut_pos);
}
 */

void ARG::remove(tuple<float, Branch, float> cut_point) {
    float pos;
    Branch center_branch;
    float t;
    tie(pos, center_branch, t) = cut_point;
    cut_time = t;
    cut_node = new_node(cut_time);
    cut_node->set_index(-2);
    Tree forward_tree = cut_tree;
    Tree backward_tree = cut_tree;
    auto f_it = recombinations.upper_bound(pos);
    auto b_it = recombinations.upper_bound(pos);
    Branch prev_joining_branch;
    Branch next_joining_branch;
    Branch prev_removed_branch = center_branch;
    Branch next_removed_branch = center_branch;
    while (next_removed_branch != Branch()) {
        Recombination &r = f_it->second;
        prev_joining_branch = forward_tree.find_joining_branch(prev_removed_branch);
        forward_tree.forward_update(r);
        next_removed_branch = r.trace_forward(t, prev_removed_branch);
        if (next_removed_branch.upper_node == root) {
            next_removed_branch = Branch();
        }
        next_joining_branch = forward_tree.find_joining_branch(next_removed_branch);
        r.remove(prev_removed_branch, next_removed_branch, prev_joining_branch, next_joining_branch, cut_node);
        removed_branches[min(r.pos, sequence_length)] = next_removed_branch;
        joining_branches[min(r.pos, sequence_length)] = next_joining_branch;
        f_it++;
        prev_removed_branch = next_removed_branch;
    }
    next_removed_branch = center_branch;
    prev_removed_branch = center_branch;
    while (prev_removed_branch != Branch()) {
        b_it--;
        Recombination &r = b_it->second;
        removed_branches[r.pos] = prev_removed_branch;
        next_joining_branch = backward_tree.find_joining_branch(next_removed_branch);
        joining_branches[r.pos] = next_joining_branch;
        backward_tree.backward_update(r);
        prev_removed_branch = r.trace_backward(t, next_removed_branch);
        if (prev_removed_branch.upper_node == root) {
            prev_removed_branch = Branch();
        }
        prev_joining_branch = backward_tree.find_joining_branch(prev_removed_branch);
        r.remove(prev_removed_branch, next_removed_branch, prev_joining_branch, next_joining_branch, cut_node);
        next_removed_branch = prev_removed_branch;
    }
    remove_empty_recombinations();
    remap_mutations();
    start = removed_branches.begin()->first;
    end = removed_branches.rbegin()->first;
    cut_tree.remove(center_branch, cut_node);
    start_tree = modify_tree_to(start, cut_tree, cut_pos);
}

void ARG::remove(map<float, Branch> seed_branches) {
    // insights here: the coordinates of removed branches is the same as recombinations
    Tree tree = Tree();
    auto recomb_it = recombinations.lower_bound(start);
    auto seed_it = seed_branches.begin();
    Branch prev_removed_branch = Branch();
    Branch next_removed_branch = Branch();
    Branch prev_joining_branch = Branch();
    Branch next_joining_branch = Branch();
    while (recomb_it->first < end) {
        next_removed_branch = seed_it->second;
        seed_it++;
        Recombination &r = recomb_it->second;
        tree.forward_update(r);
        next_joining_branch = tree.find_joining_branch(next_removed_branch);
        recomb_it++;
        r.remove(prev_removed_branch, next_removed_branch, prev_joining_branch, next_joining_branch);
        removed_branches[r.pos] = next_removed_branch;
        joining_branches[r.pos] = next_joining_branch;
        prev_removed_branch = next_removed_branch;
        prev_joining_branch = next_joining_branch;
    }
    removed_branches[end] = next_removed_branch;
    joining_branches[end] = next_joining_branch;
    remove_empty_recombinations();
    remap_mutations();
    start = removed_branches.begin()->first;
    end = removed_branches.rbegin()->first;
}

void ARG::remove_leaf(int index) {
    Node_ptr s = nullptr;
    for (Node_ptr n : sample_nodes) {
        if (n->index == index) {
            s = n;
        }
    }
    Tree tree = Tree();
    auto recomb_it = recombinations.begin();
    Branch removed_branch = Branch();
    Node_ptr joining_node = nullptr;
    map<float, Branch> removed_branches = {};
    while (recomb_it->first < sequence_length) {
        Recombination r = recomb_it->second;
        tree.forward_update(r);
        recomb_it++;
        joining_node = tree.parents[s];
        removed_branch = Branch(s, joining_node);
        removed_branches[r.pos] = removed_branch;
    }
    remove(removed_branches);
}

float ARG::get_updated_length() {
    float length = removed_branches.rbegin()->first - removed_branches.begin()->first;
    return length;
}

void ARG::add(map<float, Branch> &new_joining_branches, map<float, Branch> &added_branches) {
    /*
    for (auto x : added_branches) {
        add_node(x.second.upper_node);
    }
     */
    auto join_it = new_joining_branches.begin();
    auto add_it = added_branches.begin();
    auto recomb_it = recombinations.lower_bound(start);
    Branch prev_joining_branch = Branch();
    Branch next_joining_branch = Branch();
    Branch prev_added_branch = Branch();
    Branch next_added_branch = Branch();
    while (add_it != added_branches.end() and add_it->first < sequence_length) {
        if (recomb_it->first == add_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            if (join_it->first == add_it->first) {
                next_joining_branch = join_it->second;
                join_it++;
            }
            next_added_branch = add_it->second;
            add_it++;
            if (prev_added_branch.upper_node == r.inserted_node) {
                if (joining_branches.begin()->first != r.pos and joining_branches.rbegin()->first != r.pos) {
                    assert(prev_joining_branch == r.target_branch);
                }
            }
            r.add(prev_added_branch, next_added_branch, prev_joining_branch, next_joining_branch, cut_node);
            prev_joining_branch = next_joining_branch;
            prev_added_branch = next_added_branch;
        } else {
            if (join_it->first == add_it->first) {
                next_joining_branch = join_it->second;
                join_it++;
            }
            next_added_branch = add_it->second;
            new_recombination(add_it->first, prev_added_branch, prev_joining_branch, next_added_branch, next_joining_branch);
            add_it++;
            prev_joining_branch = next_joining_branch;
            prev_added_branch = next_added_branch;
        }
    }
    remove_empty_recombinations();
    impute(new_joining_branches, added_branches);
    // impute_nodes(start, end);
    // impute_nodes(0, sequence_length);
    start_tree.add(added_branches.begin()->second, new_joining_branches.begin()->second, cut_node);
}

/*
void ARG::smc_sample_recombinations() {
    RSP_smc rsp = RSP_smc();
    Tree tree = Tree();
    for (auto &x : recombinations) {
        if (x.first != 0 and x.first < sequence_length) {
            rsp.sample_recombination(x.second, cut_time, tree);
            assert(x.second.start_time > 0);
        }
        tree.forward_update(x.second);
    }
}
 */

void ARG::smc_sample_recombinations() {
    RSP_smc rsp = RSP_smc();
    Tree tree = start_tree;
    auto it = recombinations.upper_bound(start);
    while (it->first < end) {
        Recombination &r = it->second;
        if (r.pos != 0 and r.pos < sequence_length) {
            rsp.sample_recombination(r, cut_time, tree);
            assert(r.start_time > 0);
        }
        tree.forward_update(r);
        it++;
    }
}

int ARG::count_incompatibility() {
    Tree tree = Tree();
    auto recomb_it = recombinations.begin();
    auto mut_it = mutation_sites.begin();
    float start = 0;
    float end = 0;
    float x = 0;
    int count = 0;
    for (int i = 0; i < bin_num; i++) {
        start = coordinates[i];
        end = coordinates[i+1];
        if (start  == recomb_it->first) {
            Recombination r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
        while (*mut_it < end) {
            x = *mut_it;
            count += count_incompatibility(tree, x);
            mut_it++;
        }
    }
    return count;
}

void ARG::write(string node_file, string branch_file, string recomb_file) {
    sort_nodes();
    write_nodes(node_file);
    write_branches(branch_file);
    write_recombs(recomb_file);
}

void ARG::read(string node_file, string branch_file) {
    read_nodes(node_file);
    read_branches(branch_file);
    start_tree = get_tree_at(0);
}

void ARG::read(string node_file, string branch_file, string recomb_file) {
    read_nodes(node_file);
    read_branches(branch_file);
    read_recombs(recomb_file);
    start_tree = get_tree_at(0);
}

// private methods:

void ARG::impute_nodes(float x, float y) {
    Tree start_tree = get_tree_at(x);
    Branch null_branch = Branch();
    Fitch_reconstruction rc = Fitch_reconstruction(start_tree);
    auto recomb_it = recombinations.upper_bound(x);
    auto mut_it = mutation_sites.lower_bound(x);
    float curr_pos = x;
    float m = 0;
    while (curr_pos < y) {
        curr_pos = recomb_it->first;
        while (*mut_it < curr_pos) {
            m = *mut_it;
            rc.reconstruct(m);
            mut_it++;
        }
        rc.update(recomb_it->second);
        recomb_it++;
    }
    return;
}

void ARG::impute(map<float, Branch> &new_joining_branches, map<float, Branch> &added_branches) {
    float start = added_branches.begin()->first;
    float end = added_branches.rbegin()->first;
    auto mut_it = mutation_sites.lower_bound(start);
    auto join_it = new_joining_branches.begin();
    auto add_it = added_branches.begin();
    float m = 0;
    Branch joining_branch = Branch();
    Branch added_branch = Branch();
    while (add_it->first < end) {
        added_branch = add_it->second;
        if (join_it->first == add_it->first) {
            joining_branch = join_it->second;
            join_it++;
        }
        add_it++;
        while (*mut_it < add_it->first) {
            m = *mut_it;
            map_mutation(m, joining_branch, added_branch);
            mut_it++;
        }
    }
}

void ARG::map_mutations(float x, float y) {
    Tree tree = get_tree_at(x);
    auto recomb_it = recombinations.upper_bound(x);
    auto mut_it = mutation_sites.lower_bound(x);
    float m = *mut_it;
    while (*mut_it < y) {
        m = *mut_it;
        while (recomb_it->first < m) {
            Recombination r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
        map_mutation(tree, m);
        mut_it++;
    }
}

void ARG::map_mutation(float x, Branch joining_branch, Branch added_branch) {
    float sl, su, s0, sm;
    Branch new_branch;
    sl = joining_branch.lower_node->get_state(x);
    su = joining_branch.upper_node->get_state(x);
    s0 = added_branch.lower_node->get_state(x);
    if (sl + su + s0 > 1) {
        sm = 1;
    } else {
        sm = 0;
    }
    added_branch.upper_node->write_state(x, sm);
    if (sl != su) {
        mutation_branches[x].erase(joining_branch);
    }
    if (sm != sl) {
        new_branch = Branch(joining_branch.lower_node, added_branch.upper_node);
        mutation_branches[x].insert(new_branch);
    }
    if (sm != su) {
        new_branch = Branch(added_branch.upper_node, joining_branch.upper_node);
        mutation_branches[x].insert(new_branch);
    }
    if (sm != s0) {
        mutation_branches[x].insert(added_branch);
    }
}

void ARG::remap_mutations() {
    float x = joining_branches.begin()->first;
    float y = joining_branches.rbegin()->first;
    auto mut_it = mutation_branches.lower_bound(x);
    auto join_it = joining_branches.begin();
    auto remove_it = removed_branches.begin();
    Node_ptr joining_node = nullptr;
    Branch joining_branch = Branch();
    Branch removed_branch = Branch();
    Branch lower_branch = Branch();
    Branch upper_branch = Branch();
    while (mut_it->first < y) {
        while (join_it->first < mut_it->first) {
            joining_branch = join_it->second;
            join_it++;
        }
        while (remove_it->first < mut_it->first) {
            removed_branch = remove_it->second;
            joining_node = removed_branch.upper_node;
            remove_it++;
        }
        lower_branch = Branch(joining_branch.lower_node, joining_node);
        upper_branch = Branch(joining_node, joining_branch.upper_node);
        if (mut_it->second.count(removed_branch) > 0) {
            mut_it->second.erase(removed_branch);
        }
        if (mut_it->second.count(lower_branch) > 0 or mut_it->second.count(upper_branch) > 0) {
            mut_it->second.erase(lower_branch);
            mut_it->second.erase(upper_branch);
            mut_it->second.insert(joining_branch);
        }
        mut_it++;
    }
}

void ARG::map_mutation(Tree tree, float x) {
    set<Branch> branches = {};
    float sl = 0;
    float su = 0;
    for (Branch b : tree.branches) {
        sl = b.lower_node->get_state(x);
        su = b.upper_node->get_state(x);
        if (sl != su) {
            branches.insert({b});
        }
    }
    mutation_branches[x] = branches;
}

/*
void ARG::clear_memory() {
    set<Node_ptr, compare_node> reduced_node_set = {};
    set<Node_ptr > deleted_nodes = {};
    Recombination &r = recombinations.begin()->second;
    for (Branch b : r.inserted_branches) {
        reduced_node_set.insert(b.lower_node);
    }
    for (auto &x : recombinations) {
        if (x.first > 0 and x.first < sequence_length) {
            reduced_node_set.insert(x.second.inserted_node);
        }
    }
    for (Node_ptr n : node_set) {
        if (reduced_node_set.count(n) == 0 and n->index != -1) {
            deleted_nodes.insert(n);
        }
    }
    for (Node_ptr n : deleted_nodes) {
        delete n;
    }
    node_set.clear();
    node_set = reduced_node_set;
}
 */

void ARG::check_mapping() {
    for (Node_ptr n : node_set) {
        assert(n->ambiguous_sites.size() == 0);
    }
}

void ARG::check_incompatibility() {
    Tree tree = Tree();
    auto recomb_it = recombinations.begin();
    auto mut_it = mutation_sites.begin();
    int total_count = 0;
    while (mut_it != prev(mutation_sites.end())) {
        float m = *mut_it;
        while (recomb_it->first < m) {
            Recombination &r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
        total_count += count_incompatibility(tree, m);
        mut_it++;
    }
    cout << "Number of incompatibilities: " << total_count << endl;
}

/*
void ARG::clear_memory(map<float, Branch> added_branches) {
    Node_ptr node = nullptr;
    for (auto x : added_branches) {
        if (x.second.upper_node != node) {
            node = x.second.upper_node;
            if (node != nullptr and node_set.count(node) == 0) {
                delete node;
            }
        }
    }
}
 */

void ARG::clear_remove_info() {
    removed_branches.clear();
    joining_branches.clear();
    start = 0;
    end = 0;
}
 
float ARG::smc_prior_likelihood(float r) {
    Tree tree = get_tree_at(0);
    float rho = 0;
    float log_likelihood = 0;
    log_likelihood += tree.prior_likelihood();
    float tree_length = tree.length();
    auto recomb_it = recombinations.upper_bound(0);
    float bin_start = 0;
    float bin_end = 0;
    for (int i = 0; i < bin_num; i++) {
        bin_start = coordinates[i];
        bin_end = coordinates[i+1];
        rho = (bin_end - bin_start)*r*Ne;
        if (bin_start == recomb_it->first) {
            Recombination r = recomb_it->second;
            recomb_it++;
            log_likelihood -= rho*tree_length;
            log_likelihood += log(rho*tree_length);
            log_likelihood += tree.transition_likelihood(r);
            tree.forward_update(r);
            tree_length = tree.length();
            assert(tree_length > 0);
        } else {
            log_likelihood -= rho*tree_length;
        }
        assert(!isnan(log_likelihood));
    }
    return log_likelihood;
}

float ARG::data_likelihood(float m) {
    impute_nodes(0, bin_num);
    float theta = 0;
    float log_likelihood = 0;
    Tree tree = get_tree_at(0);
    auto recomb_it = recombinations.upper_bound(0);
    auto mut_it = mutation_sites.begin();
    float bin_start = 0;
    float bin_end = 0;
    float prev_mut_pos = 0;
    float next_mut_pos = 0;
    for (int i = 0; i < bin_num; i++) {
        bin_start = coordinates[i];
        bin_end = coordinates[i+1];
        if (bin_start == recomb_it->first) {
            Recombination r = recomb_it->second;
            recomb_it++;
            tree.forward_update(r);
        }
        while (*mut_it < bin_end) {
            next_mut_pos = *mut_it;
            theta = m*Ne;
            log_likelihood += tree.data_likelihood(theta, next_mut_pos);
            theta = (next_mut_pos - prev_mut_pos - 1)*m*Ne;
            log_likelihood += tree.null_likelihood(theta);
            prev_mut_pos = next_mut_pos;
        }
    }
    return log_likelihood;
}

float ARG::smc_likelihood(float r, float m) {
    return smc_prior_likelihood(r) + data_likelihood(m);
}

set<float> ARG::get_check_points() {
    float start_pos = removed_branches.begin()->first;
    float end_pos = removed_branches.rbegin()->first;
    auto recomb_it = recombinations.lower_bound(start_pos);
    map<Node_ptr , float> deleted_nodes = {};
    vector<tuple<Node_ptr , float, float>> node_span = {};
    while (recomb_it->first <= end_pos) {
        Recombination r = recomb_it->second;
        deleted_nodes[r.deleted_node] = r.pos;
        Node_ptr inserted_node = r.inserted_node;
        if (deleted_nodes.count(inserted_node) > 0 and inserted_node != root) {
            node_span.push_back({inserted_node, deleted_nodes[inserted_node], r.pos});
            deleted_nodes.erase(inserted_node);
        }
        recomb_it++;
    }
    Node_ptr n;
    float x;
    float y;
    set<float> check_points = {};
    for (auto ns : node_span) {
        tie(n, x, y) = ns;
        if (!check_disjoint_nodes(x, y)) {
            check_points.insert(y);
        }
    }
    return check_points;
}

bool ARG::check_disjoint_nodes(float x, float y) {
    auto recomb_it = recombinations.lower_bound(x);
    float t = recomb_it->second.deleted_node->time;
    Branch b = recomb_it->second.merging_branch;
    while (recomb_it->first < y) {
        b = recomb_it->second.trace_forward(t, b);
        if (b == Branch()) {
            return false;
        }
        recomb_it++;
    }
    if (b != recomb_it->second.target_branch) {
        return false;
    }
    return true;
}

// private methods:

void ARG::new_recombination(float pos, Branch prev_added_branch, Branch prev_joining_branch, Branch next_added_branch, Branch next_joining_branch) {
    set<Branch> deleted_branches;
    set<Branch> inserted_branches;
    deleted_branches.insert(prev_added_branch);
    deleted_branches.insert(Branch(prev_joining_branch.lower_node, prev_added_branch.upper_node));
    deleted_branches.insert(Branch(prev_added_branch.upper_node, prev_joining_branch.upper_node));
    deleted_branches.insert(next_joining_branch);
    inserted_branches.insert(next_added_branch);
    inserted_branches.insert(Branch(next_joining_branch.lower_node, next_added_branch.upper_node));
    inserted_branches.insert(Branch(next_added_branch.upper_node, next_joining_branch.upper_node));
    inserted_branches.insert(prev_joining_branch);
    Recombination r = Recombination(deleted_branches, inserted_branches);
    r.set_pos(pos);
    recombinations[pos] = r;
    return;
}

float ARG::random() {
    // return (float) rand()/RAND_MAX;
    float p = uniform_random();
    return p;
}

void ARG::remove_empty_recombinations() {
    auto recomb_it = recombinations.begin();
    while (recomb_it != recombinations.end()) {
        Recombination &r = recomb_it->second;
        if (r.deleted_branches.size() == 0 and r.inserted_branches.size() == 0 and r.pos < sequence_length) {
            recomb_it = recombinations.erase(recomb_it);
        } else {
            ++recomb_it;
        }
    }
}

int ARG::count_incompatibility(Tree tree, float x) {
    int count = -1;
    for (Branch b : tree.branches) {
        if (b.upper_node->index >= 0) {
            int i1 = b.upper_node->get_state(x);
            int i2 = b.lower_node->get_state(x);
            if (i1 != i2) {
                count += 1;
            }
        }
    }
    return max(0, count);
}

void ARG::sort_nodes() {
    int index = 0;
    for (Node_ptr n : node_set) {
        n->set_index(index);
        index += 1;
    }
}

void ARG::write_nodes(string filename) {
    node_set.clear();
    for (auto x : recombinations) {
        if (x.second.pos > 0 and x.second.pos < sequence_length) {
            add_node(x.second.inserted_node);
        } else {
            for (Branch b : x.second.inserted_branches) {
                add_node(b.upper_node);
            }
        }
    }
    for (Node_ptr n : sample_nodes) {
        add_node(n);
    }
    ofstream file;
    file.open(filename);
    int index = 0;
    for (Node_ptr n : node_set) {
        n->set_index(index);
        file << setprecision(15) << n->time*Ne << "\n";
        index += 1;
    }
    file.close();
    node_set.clear();
}

void ARG::write_branches(string filename) {
    map<Branch, int> branch_map;
    vector<tuple<int, int, float, float>> branch_info;
    int pos;
    for (auto x : recombinations) {
        if (x.first < sequence_length) {
            pos = x.first;
            Recombination &r = x.second;
            for (Branch b : r.inserted_branches) {
                branch_map[b] = pos;
            }
            for (Branch b : r.deleted_branches) {
                // assert((node_set.count(b.lower_node) > 0 and node_set.count(b.upper_node) > 0) or b.upper_node == root);
                int k1 = b.upper_node->index;
                int k2 = b.lower_node->index;
                assert(k1 < 1e5 and k2 < 1e5);
                branch_info.push_back({k1, k2, branch_map.at(b), pos});
                branch_map.erase(b);
            }
        }
    }
    for (auto x : branch_map) {
        Branch b = x.first;
        // assert((node_set.count(b.lower_node) > 0 and node_set.count(b.upper_node) > 0) or b.upper_node == root);
        int k1 = b.upper_node->index;
        int k2 = b.lower_node->index;
        assert(k1 < 1e5 and k2 < 1e5);
        branch_info.push_back({k1, k2, x.second, sequence_length});
    }
    sort(branch_info.begin(), branch_info.end(), compare_edge);
    ofstream file;
    file.open(filename);
    file << std::setprecision(20) << std::fixed;
    for (int i = 0; i < branch_info.size(); i++) {
        auto [k1, k2, x, l] = branch_info[i];
        assert(k1 < 1e5 and k2 < 1e5);
        file << x << " " << l << " " << k1 << " " << k2 << "\n";
    }
    file.close();
}

void ARG::write_recombs(string filename) {
    ofstream file;
    file.open(filename);
    for (auto x : recombinations) {
        Recombination &r = x.second;
        if (x.first > 0 and x.first < sequence_length) {
            file << r.pos << " " << r.source_branch.lower_node->index << " " << r.source_branch.upper_node->index << " " << Ne*r.start_time << endl;
        }
    }
    file.close();
}

void ARG::read_nodes(string filename) {
    root->set_index(-1);
    node_set.clear();
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    float x;
    while (fin >> x) {
        add_new_node(x/Ne);
    }
}

void ARG::read_branches(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    vector<Node_ptr> nodes = vector<Node_ptr>(node_set.begin(), node_set.end());
    float x;
    float y;
    float p;
    float c;
    float left;
    float right;
    Node_ptr parent_node;
    Node_ptr child_node;
    Branch b;
    map<float, set<Branch>> deleted_branches = {{0, {}}};
    map<float, set<Branch>> inserted_branches = {};
    while (fin >> x >> y >> p >> c) {
        left = x;
        right = y;
        if (p < 0) {
            parent_node = root;
        } else {
            parent_node = nodes[int(p)];
        }
        child_node = nodes[int(c)];
        b = Branch(child_node, parent_node);
        if (deleted_branches.count(right) > 0) {
            deleted_branches[right].insert(b);
        } else {
            deleted_branches[right] = {b};
        }
        if (inserted_branches.count(left) > 0) {
            inserted_branches[left].insert(b);
        } else {
            inserted_branches[left] = {b};
        }
    }
    deleted_branches.erase(sequence_length);
    for (auto x : deleted_branches) {
        float pos = x.first;
        set<Branch> db = deleted_branches.at(pos);
        set<Branch> ib = inserted_branches.at(pos);
        Recombination r = Recombination(db, ib);
        r.set_pos(pos);
        recombinations[pos] = r;
    }
}

void ARG::read_recombs(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    vector<Node_ptr> nodes = vector<Node_ptr>(node_set.begin(), node_set.end());
    map<int, Branch> source_branches = {};
    map<int, float> start_times = {};
    float pos;
    int n1;
    int n2;
    float t;
    Node_ptr ln;
    Node_ptr un;
    Branch b;
    while (fin >> pos >> n1 >> n2 >> t) {
        ln = nodes[n1];
        if (n2 == -1) {
            un = root;
        } else {
            un = nodes[n2];
        }
        b = Branch(nodes[n1], nodes[n2]);
        source_branches[pos] = b;
        start_times[pos] = t;
    }
    for (auto &x : recombinations) {
        pos = x.first;
        if (pos > 0 and pos < sequence_length) {
            t = start_times.at(pos);
            b = source_branches.at(pos);
            if (pos > 0 and pos < sequence_length) {
                x.second.start_time = t/Ne;
                x.second.source_branch = b;
                x.second.find_nodes();
                x.second.find_target_branch();
                x.second.find_recomb_info();
            }
        }
    }
}

float ARG::get_arg_length() {
    Tree tree = get_tree_at(0);
    auto recomb_it = recombinations.upper_bound(0);
    float l = 0, span = 0;
    float tree_length = tree.length();
    float prev_pos = 0;
    float next_pos = 0;
    while (next_pos <= sequence_length) {
        next_pos = recomb_it->first;
        span = min(sequence_length, next_pos) - prev_pos;
        l += tree_length*span;
        Recombination &r = recomb_it->second;
        recomb_it++;
        tree.forward_update(r);
        tree_length = tree.length();
        prev_pos = next_pos;
    }
    return l;
}

float ARG::get_arg_length(float x, float y) {
    Tree tree = start_tree;
    auto recomb_it = recombinations.upper_bound(x);
    float l = 0, span = 0;
    float tree_length = tree.length();
    float prev_pos = x;
    float next_pos = x;
    while (next_pos <= y) {
        next_pos = recomb_it->first;
        span = min(sequence_length, next_pos) - prev_pos;
        l += tree_length*span;
        Recombination &r = recomb_it->second;
        recomb_it++;
        tree.forward_update(r);
        tree_length = tree.length();
        prev_pos = next_pos;
    }
    return l;
}

float ARG::get_arg_length(map<float, Branch> &new_joining_branches, map<float, Branch> &new_added_branches) {
    float x = new_added_branches.begin()->first;
    float y = new_added_branches.rbegin()->first;
    float l = get_arg_length(x, y);
    auto add_it = new_added_branches.begin();
    auto join_it = new_joining_branches.begin();
    Branch joining_branch;
    float span = 0, h = 0, join_time = 0;
    while (add_it->first < y) {
        span = next(add_it)->first - add_it->first;
        joining_branch = join_it->second;
        join_time = add_it->second.upper_node->time;
        if (joining_branch.upper_node == root) {
            h = 2*join_time - joining_branch.lower_node->time - cut_time;
        } else {
            h = join_time - cut_time;
        }
        l += h*span;
        if (next(join_it)->first == next(add_it)->first) {
            join_it++;
        }
        add_it++;
    }
    return l;
}

tuple<float, Branch, float> ARG::sample_internal_cut() {
    // float arg_length = get_arg_length(0, sequence_length);
    float arg_length = get_arg_length();
    cut_tree = Tree();
    float p = random();
    p = 0.01 + 0.98*p; // smooth p away from extreme values
    float l = arg_length*p;
    auto recomb_it = recombinations.begin();
    float tree_length = 0;
    Branch branch;
    float prev_pos = 0;
    float next_pos = 0;
    float pos = 0;
    float time;
    while (next_pos < sequence_length) {
        Recombination &r = recomb_it->second;
        cut_tree.forward_update(r);
        recomb_it++;
        next_pos = recomb_it->first;
        tree_length = cut_tree.length();
        l -= tree_length*(next_pos - prev_pos);
        if (l < 0) {
            pos = 0.5*(prev_pos + next_pos);
            tie(branch, time) = cut_tree.sample_cut_point();
            cut_pos = pos;
            assert(pos < next_pos);
            assert(recombinations.count(cut_pos) == 0);
            return {pos, branch, time};
        }
        prev_pos = next_pos;
    }
    cerr << "sample internal cut failed" << endl;
    exit(1);
}

/*
tuple<float, Branch, float> ARG::sample_internal_cut() {
    auto recomb_it = recombinations.begin();
    float p = uniform_random();
    int dist = (recombinations.size() - 2)*p;
    dist = max(dist, 1);
    advance(recomb_it, dist);
    Recombination &r = recomb_it->second;
    assert(r.pos > 0 and r.pos < sequence_length);
    float x = r.pos + 1;
    float t = (r.start_time + r.recombined_branch.lower_node->time)/2;
    cut_pos = x;
    cut_tree = get_tree_at(x);
    return {x, r.recombined_branch, t};
}
 */

/*
tuple<float, Branch, float> ARG::sample_internal_cut() {
    float p = 0;
    float replace_prob = 0;
    int mapping_size = 0;
    int count = 0;
    auto mb_it = mutation_branches.begin();
    Branch b;
    float x = 0, t = 0;
    while (mb_it->first < sequence_length) {
        mapping_size = (int) mb_it->second.size();
        replace_prob = 0;
        if (mapping_size > 1) {
            replace_prob = (float) mapping_size/(mapping_size + count);
            count += mb_it->second.size();
        }
        p = uniform_random();
        if (p < replace_prob) {
            auto b_it = mb_it->second.begin();
            advance(b_it, (mapping_size - 1)*uniform_random());
            b = *b_it;
            x = mb_it->first;
            t = b.lower_node->time + 1e-3;
        }
        mb_it++;
    }
    cut_pos = x;
    cut_tree = get_tree_at(x);
    // cout << x << " " << b.lower_node->time << " " << b.upper_node->time << " " << t << endl;
    return {x, b, t};
}
 */

tuple<float, Branch, float> ARG::sample_terminal_cut() {
    Branch branch;
    float time = 1e-10;
    vector<Node_ptr > nodes = vector<Node_ptr >(sample_nodes.begin(), sample_nodes.end());
    int index = rand() % nodes.size();
    Node_ptr terminal_node = nodes[index];
    cut_tree = get_tree_at(0);
    for (Branch b : cut_tree.branches) {
        if (b.lower_node == terminal_node) {
            branch = b;
            break;
        }
    }
    return {0, branch, time};
}

bool compare_edge(const tuple<int, int, float, float>& edge1, const tuple<int, int, float, float>& edge2) {
    if (get<0>(edge1) < get<0>(edge2)) {
        return true;
    } else if (get<0>(edge1) > get<0>(edge2)) {
        return false;
    }
    if (get<1>(edge1) < get<1>(edge2)) {
        return true;
    } else if (get<1>(edge1) > get<1>(edge2)) {
        return false;
    }
    return get<2>(edge1) < get<2>(edge2);
}
