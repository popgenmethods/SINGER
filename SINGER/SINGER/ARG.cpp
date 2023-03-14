//
//  ARG.cpp
//  conditional-coalescent
//
//  Created by Yun Deng on 4/14/22.
//

#include "ARG.hpp"

ARG::ARG() {
}

ARG::ARG(Node *n, float s, float l) {
    root->set_index(-1);
    sequence_length = l;
    bin_size = s;
    bin_num = ceil((float) sequence_length/bin_size);
    Recombination r = Recombination();
    r.set_pos(INT_MAX);
    recombination_info.insert({INT_MAX, r});
    mutation_info.insert({INT_MAX, {}});
    add_sample(n);
    add_node(n);
    Recombination init_recombination = Recombination({}, {Branch(n, root)});
    init_recombination.set_pos(0);
    recombination_info.insert({0, init_recombination});
}

ARG::ARG(float s, float l) {
    sequence_length = l;
    bin_size = s;
    bin_num = ceil((float) sequence_length/bin_size);
    Recombination r = Recombination();
    r.set_pos(INT_MAX);
    recombination_info.insert({INT_MAX, r});
    mutation_info.insert({INT_MAX, {}});
}

ARG::~ARG() {
}

void ARG::set_rates(float r, float m) {
    recomb_rate = r;
    mut_rate = m;
    float rho = recomb_rate*bin_size;
    float theta = mut_rate*bin_size;
    rho_unit = recomb_rate*bin_size;
    for (int i = 0; i < bin_num; i++) {
        bin_sizes.push_back(bin_size);
        rhos.push_back(rho);
        thetas.push_back(theta);
    }
}

void ARG::add_sample(Node *n) {
    sample_nodes.insert(n);
    add_node(n);
    set<float> mutations = n->mutation_sites;
    int pos;
    for (float x : mutations) {
        pos = floor(x/bin_size);
        if (mutation_info.count(pos) > 0) {
            mutation_info.at(pos).insert(x);
        } else {
            mutation_info.insert({pos, {x}});
        }
    }
    base_nodes.clear();
    base_nodes.insert({0, n});
    base_nodes.insert({bin_num, n});
}

void ARG::add_node(Node *n) {
    if (n != root) {
        node_set.insert(n);
    }
}

Node *ARG::select_node(int index) {
    for (Node *n : node_set) {
        if (n->index == index) {
            return n;
        }
    }
    cerr << "node not found" << endl;
    exit(1);
}

Node *ARG::add_new_node(float t) {
    if (t == numeric_limits<float>::infinity()) {
        return root;
    }
    Node *new_node = new Node(t);
    new_node->set_index((int) node_set.size());
    add_node(new_node);
    if (t == 0) {
        sample_nodes.insert(new_node);
    }
    return new_node;
}

Tree ARG::get_tree_at(int x) {
    Tree tree = Tree();
    map<int, Recombination>::iterator recomb_it = recombination_info.begin();
    for (int i = 0; i <= x; i++) {
        if (i == recomb_it->first) {
            Recombination r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
    }
    return tree;
}

map<int, pair<Branch, Node *>> ARG::remove(tuple<int, Branch, float> cut_point) {
    map<int, pair<Branch, Node *>> joining_points = {};
    float pos;
    Branch center_branch;
    float t;
    tie(pos, center_branch, t) = cut_point;
    cut_time = t;
    cut_node = new Node(cut_time);
    cut_node->set_index(-2);
    add_node(cut_node);
    Tree forward_tree = get_tree_at(pos);
    Tree backward_tree = forward_tree;
    map<int, Recombination>::iterator f_it = recombination_info.upper_bound(pos);
    map<int, Recombination>::iterator b_it = recombination_info.upper_bound(pos);
    Branch prev_split_branch;
    Branch next_split_branch;
    Branch prev_removed_branch = center_branch;
    Branch next_removed_branch = center_branch;
    pair<Branch, Node *> joining_point = {};
    while (next_removed_branch != Branch()) {
        Recombination &r = f_it->second;
        prev_split_branch = forward_tree.find_split_branch(prev_removed_branch);
        forward_tree.forward_update(r);
        next_removed_branch = r.trace_forward(t, prev_removed_branch);
        next_split_branch = forward_tree.find_split_branch(next_removed_branch);
        joining_point = {next_split_branch, next_removed_branch.upper_node};
        r.remove(prev_removed_branch, next_removed_branch, prev_split_branch, next_split_branch, cut_node);
        base_nodes.insert({min(r.pos, bin_num), next_removed_branch.lower_node});
        upper_nodes.insert({min(r.pos, bin_num), next_removed_branch.upper_node});
        joining_points.insert({min(r.pos, bin_num), joining_point});
        f_it++;
        prev_removed_branch = next_removed_branch;
    }
    next_removed_branch = center_branch;
    prev_removed_branch = center_branch;
    while (prev_removed_branch != Branch()) {
        b_it--;
        Recombination &r = b_it->second;
        base_nodes.insert({b_it->first, prev_removed_branch.lower_node});
        upper_nodes.insert({b_it->first, prev_removed_branch.upper_node});
        next_split_branch = backward_tree.find_split_branch(next_removed_branch);
        joining_point = {next_split_branch, next_removed_branch.upper_node};
        joining_points.insert({b_it->first, joining_point});
        backward_tree.backward_update(r);
        prev_removed_branch = r.trace_backward(t, next_removed_branch);
        prev_split_branch = backward_tree.find_split_branch(prev_removed_branch);
        r.remove(prev_removed_branch, next_removed_branch, prev_split_branch, next_split_branch, cut_node);
        next_removed_branch = prev_removed_branch;
    }
    remove_empty_recombinations();
    int start_pos = base_nodes.begin()->first;
    int end_pos = base_nodes.rbegin()->first;
    impute_nodes(start_pos, end_pos);
    joining_points.erase(end_pos);
    joining_points.insert({end_pos, joining_points.rbegin()->second});
    return joining_points;
}


void ARG::add(map<int, pair<Branch, Node*>> joining_points) {
    for (auto x : joining_points) {
        add_node(x.second.second);
    }
    int start_pos = base_nodes.begin()->first;
    int end_pos = base_nodes.rbegin()->first;
    map<int, Node*>::iterator base_it = base_nodes.begin();
    map<int, pair<Branch, Node*>>::iterator join_it = joining_points.begin();
    Node *base_node = base_nodes.begin()->second;
    Node *init_node = joining_points.begin()->second.second;
    Node *next_joining_node = nullptr;
    Branch prev_joining_branch = Branch();
    Branch next_joining_branch;
    Branch prev_added_branch = Branch();
    Branch next_added_branch = Branch(base_node, init_node);
    for (int i = start_pos; i <= end_pos; i++) {
        if (i == base_it->first) {
            base_node = base_it->second;
            base_it++;
        }
        if (i == join_it->first) {
            tie(next_joining_branch, next_joining_node) = join_it->second;
            join_it++;
            next_added_branch = Branch(base_node, next_joining_node);
            if (recombination_info.count(i) > 0) {
                Recombination &r = recombination_info.at(i);
                if (i == bin_num) {
                    break;
                } else if (i == start_pos) {
                    r.fix_front(next_added_branch, next_joining_branch, cut_node);
                } else if (i == end_pos) {
                    r.fix_end(prev_added_branch, prev_joining_branch, cut_node);
                } else {
                    r.add(prev_added_branch, next_added_branch, prev_joining_branch, next_joining_branch, cut_node);
                }
                assert(r.deleted_branches.size() == r.inserted_branches.size() or r.pos == 0 or r.pos >= bin_num);
                assert(r.pos == 0 or (node_set.count(r.deleted_node) > 0 and node_set.count(r.inserted_node) > 0));
                prev_joining_branch = next_joining_branch;
                prev_added_branch = next_added_branch;
            } else if (i != bin_num) {
                new_recombination(i, prev_added_branch, prev_joining_branch, next_added_branch, next_joining_branch);
                prev_joining_branch = next_joining_branch;
                prev_added_branch = next_added_branch;
            }
        }
    }
    remove_empty_recombinations();
    for (auto x : recombination_info) {
        assert(x.first == 0 or x.first >= bin_num or x.second.deleted_branches.size() == 3 or x.second.inserted_branches.size() == 4);
    }
    impute_nodes(start_pos, end_pos);
}

void ARG::sample_recombinations() {
    RSP rsp = RSP();
    Tree tree = Tree();
    for (auto &x : recombination_info) {
        if (x.second.deleted_branches.size() != 0 and x.second.pos < bin_num) {
            rsp.exact_sample_recombination(x.second, cut_time, tree);
        }
        tree.forward_update(x.second);
    }
}

void ARG::check_incompatibility() {
    Tree tree = Tree();
    map<int, Recombination>::iterator recomb_it = recombination_info.begin();
    map<int, set<float>>::iterator mut_it = mutation_info.begin();
    int total_count = 0;
    for (int i = 0; i < bin_num; i++) {
        if (i == recomb_it->first) {
            Recombination r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
        if (i == mut_it->first) {
            set<float> mutations = mut_it->second;
            for (float x : mutations) {
                int count = count_incompatibility(tree, x);
                total_count += count;
                /*
                if (count > 0) {
                    cout << "incompatibility at " << x << endl;
                }
                 */
            }
            mut_it++;
        }
    }
    cout << "Number of incompatibilities: " << total_count << endl;
}

int ARG::count_incompatibility() {
    Tree tree = Tree();
    map<int, Recombination>::iterator recomb_it = recombination_info.begin();
    map<int, set<float>>::iterator mut_it = mutation_info.begin();
    int count = 0;
    for (int i = 0; i < bin_num; i++) {
        if (i == recomb_it->first) {
            Recombination r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
        }
        if (i == mut_it->first) {
            set<float> mutations = mut_it->second;
            for (float x : mutations) {
                count += count_incompatibility(tree, x);
            }
            mut_it++;
        }
    }
    return count;
}

void ARG::write_mutations(string mut_file) {
    ofstream file;
    file.open(mut_file);
    for (int i = 0; i < bin_num; i++) {
        if (mutation_info.count(i) > 0) {
            file << mutation_info.at(i).size() << endl;
        } else {
            file << 0 << endl;
        }
    }
    file.close();
}

void ARG::write(float Ne, string node_file, string branch_file, string recomb_file) {
    write_nodes(node_file);
    write_branches(branch_file);
    write_recombs(recomb_file);
}

void ARG::read(string node_file, string branch_file) {
    read_nodes(node_file);
    read_branches(branch_file);
}

void ARG::read(string node_file, string branch_file, string recomb_file) {
    read_nodes(node_file);
    read_branches(branch_file);
    read_recombs(recomb_file);
}

// private methods:

void ARG::impute_nodes(int x, int y) {
    Tree start_tree = get_tree_at(x);
    Branch null_branch = Branch();
    Fitch_reconstruction rc = Fitch_reconstruction(start_tree);
    map<int, Recombination>::iterator recomb_it = recombination_info.upper_bound(x);
    map<int, set<float>>::iterator mut_it = mutation_info.lower_bound(x);
    for (int i = x; i < y; i++) {
        if (i == recomb_it->first) {
            rc.update(recomb_it->second);
            recomb_it++;
        }
        if (i == mut_it->first) {
            set<float> mutations = mut_it->second;
            for (float pos : mutations) {
                rc.reconstruct(pos);
            }
            mut_it++;
        }
    }
    return;
}

void ARG::clear_memory() {
    set<Node*, compare_node> reduced_node_set = {};
    set<Node *> deleted_nodes = {};
    Recombination r = recombination_info.begin()->second;
    for (Branch b : r.inserted_branches) {
        reduced_node_set.insert(b.lower_node);
    }
    for (auto x : recombination_info) {
        if (x.first > 0 and x.first < bin_num) {
            reduced_node_set.insert(x.second.inserted_node);
        }
    }
    for (Node *n : node_set) {
        if (reduced_node_set.count(n) == 0 and n->index != -1) {
            deleted_nodes.insert(n);
        }
    }
    for (Node *n : deleted_nodes) {
        delete n;
    }
    node_set.clear();
    node_set = reduced_node_set;
}

void ARG::check_mapping() {
    for (Node *n : node_set) {
        assert(n->ambiguous_sites.size() == 0);
    }
}

void ARG::clear_memory(map<int, pair<Branch, Node *>> joining_points) {
    Node *node = nullptr;
    for (auto x : joining_points) {
        if (x.second.second != node) {
            node = x.second.second;
            if (node != root and node_set.count(node) == 0) {
                delete node;
            }
        }
    }
}

void ARG::clear_upper_nodes() {
    upper_nodes.clear();
}

void ARG::clear_base_nodes() {
    base_nodes.clear();
}
 
float ARG::prior_likelihood() {
    Tree tree = get_tree_at(0);
    float log_likelihood = 0;
    log_likelihood += tree.prior_likelihood();
    float tree_length = tree.length();
    map<int, Recombination>::iterator recomb_it = recombination_info.upper_bound(0);
    log_likelihood -= rhos[0]*tree_length;
    for (int i = 0; i < bin_num; i++) {
        if (i == recomb_it->first) {
            Recombination r = recomb_it->second;
            recomb_it++;
            log_likelihood -= rhos[i]*tree_length;
            log_likelihood += log(rhos[i]*tree_length);
            log_likelihood += tree.transition_likelihood(r);
            tree.forward_update(r);
            tree_length = tree.length();
            assert(tree_length > 0);
        } else {
            log_likelihood -= rhos[i]*tree_length;
        }
        assert(!isnan(log_likelihood));
    }
    return log_likelihood;
}

float ARG::data_likelihood() {
    impute_nodes(0, bin_num);
    float log_likelihood = 0;
    Tree tree = get_tree_at(0);
    map<int, Recombination>::iterator recomb_it = recombination_info.upper_bound(0);
    map<int, set<float>>::iterator mut_it = mutation_info.begin();
    for (int i = 0; i < bin_num; i++) {
        if (i == recomb_it->first) {
            Recombination r = recomb_it->second;
            recomb_it++;
            tree.forward_update(r);
        }
        if (i == mut_it->first) {
            set<float> mutations = mut_it->second;
            mut_it++;
            log_likelihood += tree.data_likelihood(thetas[i], bin_sizes[i], mutations);
        } else {
            log_likelihood += tree.null_likelihood(thetas[i]);
        }
    }
    return log_likelihood;
}

float ARG::likelihood() {
    return prior_likelihood() + data_likelihood();
}

int ARG::get_update_length() {
    if (upper_nodes.size() == 0) {
        return 0;
    }
    int update_length = 0;
    int pos = upper_nodes.begin()->first;
    int end_pos = upper_nodes.rbegin()->first;
    map<int, Node *>::iterator remove_it = upper_nodes.begin();
    while (pos < end_pos) {
        assert(remove_it->second != nullptr);
        if (remove_it->second != root) {
            remove_it++;
            update_length += remove_it->first - pos;
        } else {
            remove_it++;
        }
        pos = remove_it->first;
    }
    return update_length;
}

bool ARG::node_check() {
    Node *n;
    Branch b;
    map<Node *, Branch> node_branch = {};
    set<Node *> deleted_nodes = {};
    // int start_pos = base_nodes.begin()->first;
    // int end_pos = base_nodes.rbegin()->first;
    for (auto x : recombination_info) {
        if (x.first <= 0 or x.first >= bin_num) {
            continue;
        }
        deleted_nodes.clear();
        Recombination r = x.second;
        for (auto y : node_branch) {
            if (y.first == r.inserted_node and y.second != r.target_branch) {
                // return false;
                assert(false);
            }
        }
        if (node_branch.count(r.inserted_node) > 0) {
            deleted_nodes.insert(r.inserted_node);
        }
        for (auto &y : node_branch) {
            b = r.trace_forward(y.first->time, y.second);
            if (b == Branch()) {
                deleted_nodes.insert(y.first);
            }
            y.second = b;
        }
        for (Node *dn : deleted_nodes) {
            node_branch.erase(dn);
        }
        n = x.second.deleted_node;
        b = x.second.merging_branch;
        assert(b != Branch());
        assert(node_branch.count(n) == 0);
        node_branch.insert({n, b});
    }
    return true;
}

set<int> ARG::get_check_points() {
    int start_pos = base_nodes.begin()->first;
    int end_pos = base_nodes.rbegin()->first;
    map<int, Recombination>::iterator recomb_it = recombination_info.lower_bound(start_pos);
    map<Node *, int> deleted_nodes = {};
    vector<tuple<Node *, int, int>> node_span = {};
    while (recomb_it->first <= end_pos) {
        Recombination r = recomb_it->second;
        deleted_nodes.insert({r.deleted_node, r.pos});
        Node *inserted_node = r.inserted_node;
        if (deleted_nodes.count(inserted_node) > 0 and inserted_node != root) {
            node_span.push_back({inserted_node, deleted_nodes.at(inserted_node), r.pos});
            deleted_nodes.erase(inserted_node);
        }
        recomb_it++;
    }
    Node *n;
    int x;
    int y;
    set<int> check_points = {};
    for (auto ns : node_span) {
        tie(n, x, y) = ns;
        if (!check_disjoint_nodes(x, y)) {
            check_points.insert(y);
        }
    }
    return check_points;
}

bool ARG::check_disjoint_nodes(int x, int y) {
    map<int, Recombination>::iterator recomb_it = recombination_info.lower_bound(x);
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
        // cout << "Disjoint nodes warning" << endl;
        return false;
    }
    return true;
}

// private methods:

void ARG::new_recombination(int pos, Branch prev_added_branch, Branch prev_joining_branch, Branch next_added_branch, Branch next_joining_branch) {
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
    // assert(r.inserted_branches.count(r.merging_branch) > 0);
    recombination_info.insert({pos, r});
    return;
}

float ARG::random() {
    return (float) rand()/RAND_MAX;
}

void ARG::remove_empty_recombinations() {
    map<int, Recombination>::iterator recomb_it = recombination_info.begin();
    while (recomb_it != recombination_info.end()) {
        Recombination r = recomb_it->second;
        if (r.deleted_branches.size() == 0 and r.inserted_branches.size() == 0 and r.pos < bin_num) {
            recomb_it = recombination_info.erase(recomb_it);
        } else {
            ++recomb_it;
        }
    }
}

int ARG::count_incompatibility(Tree tree, float x) {
    int count = -1;
    for (Branch b : tree.branches) {
        if (b.upper_node->index > 0) {
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
    for (Node *n : node_set) {
        n->set_index(index);
        index += 1;
    }
}

void ARG::write_nodes(string filename) {
    ofstream file;
    file.open(filename);
    int index = 0;
    for (Node *n : node_set) {
        n->set_index(index);
        file << setprecision(10) << n->time << "\n";
        index += 1;
    }
    file.close();
}

void ARG::write_branches(string filename) {
    map<Branch, int> branch_map;
    vector<tuple<int, int, float, float>> branch_info;
    int pos;
    for (auto x : recombination_info) {
        if (x.first != bin_num) {
            pos = x.first;
            Recombination r = x.second;
            for (Branch b : r.inserted_branches) {
                branch_map.insert({b, pos*bin_size});
            }
            for (Branch b : r.deleted_branches) {
                assert((node_set.count(b.lower_node) > 0 and node_set.count(b.upper_node) > 0) or b.upper_node == root);
                int k1 = b.upper_node->index;
                int k2 = b.lower_node->index;
                assert(k1 < 1e4 and k2 < 1e4);
                branch_info.push_back({k1, k2, branch_map.at(b), pos*bin_size});
                branch_map.erase(b);
            }
        }
    }
    for (auto x : branch_map) {
        Branch b = x.first;
        assert((node_set.count(b.lower_node) > 0 and node_set.count(b.upper_node) > 0) or b.upper_node == root);
        int k1 = b.upper_node->index;
        int k2 = b.lower_node->index;
        assert(k1 < 1e4 and k2 < 1e4);
        branch_info.push_back({k1, k2, x.second, sequence_length});
    }
    sort(branch_info.begin(), branch_info.end(), compare_edge);
    ofstream file;
    file.open(filename);
    for (int i = 0; i < branch_info.size(); i++) {
        auto [k1, k2, x, l] = branch_info[i];
        assert(k1 < 1e4 and k2 < 1e4);
        file << x << " " << l << " " << k1 << " " << k2 << "\n";
    }
    file.close();
}

void ARG::write_recombs(string filename) {
    ofstream file;
    file.open(filename);
    for (auto x : recombination_info) {
        Recombination r = x.second;
        if (x.first > 0 and x.first < bin_num) {
            file << r.pos << " " << r.source_branch.lower_node->index << " " << r.source_branch.upper_node->index << " " << r.start_time << endl;
        }
    }
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
        add_new_node(x);
    }
}

void ARG::read_branches(string filename) {
    recombination_info.clear();
    Recombination r = Recombination();
    r.set_pos(INT_MAX);
    recombination_info.insert({INT_MAX, r});
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    vector<Node*> nodes = vector<Node*>(node_set.begin(), node_set.end());
    float x;
    float y;
    float p;
    float c;
    int left;
    int right;
    Node *parent_node;
    Node *child_node;
    Branch b;
    map<int, set<Branch>> deleted_branches = {{0, {}}};
    map<int, set<Branch>> inserted_branches = {};
    while (fin >> x >> y >> p >> c) {
        left = round(x/bin_size);
        right = round(y/bin_size);
        if (p < 0) {
            parent_node = root;
        } else {
            parent_node = nodes[int(p)];
        }
        child_node = nodes[int(c)];
        b = Branch(child_node, parent_node);
        if (deleted_branches.count(right) > 0) {
            deleted_branches.at(right).insert(b);
        } else {
            deleted_branches.insert({right, {b}});
        }
        if (inserted_branches.count(left) > 0) {
            inserted_branches.at(left).insert(b);
        } else {
            inserted_branches.insert({left, {b}});
        }
    }
    deleted_branches.erase(bin_num);
    for (auto x : deleted_branches) {
        int pos = x.first;
        set<Branch> db = deleted_branches.at(pos);
        set<Branch> ib = inserted_branches.at(pos);
        Recombination r = Recombination(db, ib);
        r.set_pos(pos);
        recombination_info.insert({pos, r});
    }
}

void ARG::read_recombs(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    vector<Node*> nodes = vector<Node*>(node_set.begin(), node_set.end());
    map<int, Branch> source_branches = {};
    map<int, float> start_times = {};
    int pos;
    int n1;
    int n2;
    float t;
    Node *ln;
    Node *un;
    Branch b;
    while (fin >> pos >> n1 >> n2 >> t) {
        ln = nodes[n1];
        if (n2 == -1) {
            un = root;
        } else {
            un = nodes[n2];
        }
        b = Branch(nodes[n1], nodes[n2]);
        source_branches.insert({pos, b});
        start_times.insert({pos, t});
    }
    for (auto &x : recombination_info) {
        pos = x.first;
        if (pos > 0 and pos < bin_num) {
            t = start_times.at(pos);
            b = source_branches.at(pos);
            if (pos > 0 and pos < bin_num) {
                x.second.start_time = t;
                x.second.source_branch = b;
                x.second.find_nodes();
                x.second.find_target_branch();
                x.second.find_recomb_info();
            }
        }
    }
}

float ARG::get_arg_length() {
    Tree tree = Tree();
    map<int, Recombination>::iterator recomb_it = recombination_info.begin();
    Recombination r;
    float l = 0;
    float tree_length = 0;
    for (int i = 0; i < bin_num; i++) {
        if (recomb_it->first == i) {
            r = recomb_it->second;
            recomb_it++;
            tree.forward_update(r);
            tree_length = tree.length();
        }
        l += tree_length;
    }
    return l;
}

float ARG::get_arg_length(int x, int y) {
    Tree tree = get_tree_at(x);
    map<int, Recombination>::iterator recomb_it = recombination_info.upper_bound(x);
    Recombination r;
    float l = 0;
    float tree_length = tree.length();
    for (int i = x; i < y; i++) {
        if (recomb_it->first == i) {
            r = recomb_it->second;
            recomb_it++;
            tree.forward_update(r);
            tree_length = tree.length();
        }
        l += tree_length;
    }
    return l;
}

tuple<int, Branch, float> ARG::sample_internal_cut() {
    float arg_length = get_arg_length();
    float p = (float) rand()/RAND_MAX;
    p = 0.01 + 0.98*p;
    float l = arg_length*p;
    Tree tree = Tree();
    map<int, Recombination>::iterator recomb_it = recombination_info.begin();
    float tree_length = 0;
    Recombination r;
    Branch branch;
    float time;
    for (int i = 0; i < bin_num; i++) {
        if (recomb_it->first == i) {
            r = recomb_it->second;
            recomb_it++;
            tree.forward_update(r);
            tree_length = tree.length();
        }
        l -= tree_length;
        if (l < 0) {
            tie(branch, time) = tree.sample_cut_point();
            return {i, branch, time};
        }
    }
    // time = branch.lower_node->time + 1e-6;
    cerr << "sample internal cut failed" << endl;
    exit(1);
    //return {0, Branch(), 0};
}

/*
tuple<int, Branch, float> ARG::sample_internal_cut() {
    float p = random();
    int pos = floor(bin_num*p);
    Branch branch;
    float time;
    Tree tree = get_tree_at(pos);
    float q = random();
    q *= (tree.branches.size() - 1);
    for (Branch b : tree.branches) {
        if (b.upper_node->time < numeric_limits<float>::infinity()) {
            q -= 1;
        }
        if (q < 0) {
            branch = b;
            break;
        }
    }
    // time = branch.lower_node->time + 1e-6;
    tie(branch, time) = tree.sample_cut_point();
    // time = branch.lower_node->time + 1e-6;
    return {pos, branch, time};
}
 */

/*
tuple<int, Branch, float> ARG::sample_internal_cut() {
    int pos = 0;
    Branch branch;
    float time;
    float q = random()*(recombination_info.size() - 2);
    for (auto x : recombination_info) {
        Recombination r = x.second;
        if (r.pos != 0) {
            q -= 1;
        }
        if (q < 0) {
            branch = r.recombined_branch;
            pos = r.pos;
            break;
        }
    }
    assert(pos < bin_num and pos > 0);
    Tree tree = get_tree_at(pos);
    time = branch.lower_node->time + 1e-6;
    return {pos, branch, time};
}
 */

tuple<int, Branch, float> ARG::sample_terminal_cut() {
    Branch branch;
    float time = 1e-10;
    vector<Node *> nodes = vector<Node *>(sample_nodes.begin(), sample_nodes.end());
    int index = rand() % nodes.size();
    // index = 0;
    Node *terminal_node = nodes[index];
    Tree start_tree = get_tree_at(0);
    for (Branch b : start_tree.branches) {
        if (b.lower_node == terminal_node) {
            branch = b;
            break;
        }
    }
    return {0, branch, time};
}

tuple<int, Branch, float> ARG::get_terminal_cut(int i) {
    Branch branch;
    float time = 1e-10;
    Tree start_tree = get_tree_at(0);
    for (Branch b : start_tree.branches) {
        if (b.lower_node->index == i) {
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
