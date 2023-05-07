//
//  Recombination.cpp
//  SINGER
//
//  Created by Yun Deng on 4/9/22.
//

#include "Recombination.hpp"

Recombination::Recombination() {
}

Recombination::Recombination(set<Branch> db, set<Branch> ib) {
    deleted_branches = db;
    inserted_branches = ib;
    simplify_branches();
    find_nodes();
}

void Recombination::set_pos(float x) {
    pos = x;
}

bool Recombination::affect(Branch b) {
    if (deleted_branches.count(b) > 0) {
        return true;
    }
    return false;
}

bool Recombination::create(Branch b) {
    if (inserted_branches.count(b) > 0) {
        return true;
    }
    return false;
}

Branch Recombination::trace_forward(float t, Branch curr_branch) {
    if (pos == 0 or pos == INT_MAX) {
        return Branch();
    }
    if (!affect(curr_branch)) {
        return curr_branch;
    }
    if (curr_branch == source_branch) {
        if (t >= start_time) {
            return Branch();
        } else {
            return recombined_branch;
        }
    } else if (curr_branch == target_branch) {
        if (t >= inserted_node->time) {
            return upper_transfer_branch;
        } else {
            return lower_transfer_branch;
        }
    } else {
        return merging_branch;
    }
}

Branch Recombination::trace_backward(float t, Branch curr_branch) {
    if (deleted_branches.size() == 0) {
        return Branch();
    }
    if (!create(curr_branch)) {
        return curr_branch;
    }
    if (curr_branch == recombined_branch) {
        if (t >= start_time) {
            return Branch();
        } else {
            return source_branch;
        }
    } else if (curr_branch != merging_branch) {
        return target_branch;
    } else {
        if (t >= deleted_node->time) {
            return search_lower_node(deleted_node);
        } else {
            return search_upper_node(deleted_node);
        }
    }
}

void Recombination::remove(Branch prev_removed_branch, Branch next_removed_branch, Branch prev_split_branch, Branch next_split_branch, Node *cut_node) {
    if (deleted_branches.size() == 0 and inserted_branches.size() == 0) {
        return;
    }
    if (prev_removed_branch == Branch()) {
        break_front(next_removed_branch, next_split_branch, cut_node);
        return;
    } else if (next_removed_branch == Branch()) {
        break_end(prev_removed_branch, prev_split_branch, cut_node);
        return;
    }
    add_deleted_branch(prev_split_branch);
    add_deleted_branch(next_removed_branch);
    add_deleted_branch(Branch(next_split_branch.lower_node, next_removed_branch.upper_node));
    add_deleted_branch(Branch(next_removed_branch.upper_node, next_split_branch.upper_node));
    add_inserted_branch(next_split_branch);
    add_inserted_branch(prev_removed_branch);
    add_inserted_branch(Branch(prev_split_branch.lower_node, prev_removed_branch.upper_node));
    add_inserted_branch(Branch(prev_removed_branch.upper_node, prev_split_branch.upper_node));
    add_deleted_branch(Branch(prev_removed_branch.lower_node, cut_node));
    add_inserted_branch(Branch(next_removed_branch.lower_node, cut_node));
    simplify_branches();
    if (source_branch == Branch(prev_split_branch.lower_node, prev_removed_branch.upper_node) or source_branch == Branch(prev_removed_branch.upper_node, prev_split_branch.upper_node)) { // when the previous source branch was destroyed
        source_branch = prev_split_branch;
    }
    find_nodes();
    find_target_branch();
    if (deleted_branches.size() > 0) {
        find_recomb_info();
    }
    assert(deleted_branches.size() == 0 or deleted_branches.count(source_branch) > 0);
}

void Recombination::remove(Branch prev_removed_branch, Branch next_removed_branch, Branch prev_split_branch, Branch next_split_branch) {
    add_deleted_branch(prev_split_branch);
    add_deleted_branch(next_removed_branch);
    add_deleted_branch(Branch(next_split_branch.lower_node, next_removed_branch.upper_node));
    add_deleted_branch(Branch(next_removed_branch.upper_node, next_split_branch.upper_node));
    add_inserted_branch(next_split_branch);
    add_inserted_branch(prev_removed_branch);
    add_inserted_branch(Branch(prev_split_branch.lower_node, prev_removed_branch.upper_node));
    add_inserted_branch(Branch(prev_removed_branch.upper_node, prev_split_branch.upper_node));
    simplify_branches();
    if (deleted_branches.size() == 0 and inserted_branches.size() == 0) {
        return;
    }
    if (source_branch == Branch(prev_split_branch.lower_node, prev_removed_branch.upper_node) or source_branch == Branch(prev_removed_branch.upper_node, prev_split_branch.upper_node)) { // when the previous source branch was destroyed
        source_branch = prev_split_branch;
    }
    find_nodes();
    find_target_branch();
    find_recomb_info();
}

void Recombination::add(Branch prev_added_branch, Branch next_added_branch, Branch prev_joining_branch, Branch next_joining_branch, Node *cut_node) {
    if (prev_added_branch == next_added_branch and prev_joining_branch == next_joining_branch) {
        return;
    }
    if (next_added_branch != Branch()) {
        add_inserted_branch(next_added_branch);
        add_inserted_branch(Branch(next_joining_branch.lower_node, next_added_branch.upper_node));
        add_inserted_branch(Branch(next_added_branch.upper_node, next_joining_branch.upper_node));
        add_deleted_branch(next_joining_branch);
        if (cut_node != nullptr) {
            add_deleted_branch(Branch(next_added_branch.lower_node, cut_node));
        }
    }
    if (prev_added_branch != Branch()) {
        add_deleted_branch(prev_added_branch);
        add_deleted_branch(Branch(prev_joining_branch.lower_node, prev_added_branch.upper_node));
        add_deleted_branch(Branch(prev_added_branch.upper_node, prev_joining_branch.upper_node));
        add_inserted_branch(prev_joining_branch);
        if (cut_node != nullptr) {
            add_inserted_branch(Branch(prev_added_branch.lower_node, cut_node));
        }
    }
    /*
    add_deleted_branch(prev_added_branch);
    add_deleted_branch(next_joining_branch);
    add_deleted_branch(Branch(prev_joining_branch.lower_node, prev_added_branch.upper_node));
    add_deleted_branch(Branch(prev_added_branch.upper_node, prev_joining_branch.upper_node));
    add_inserted_branch(next_added_branch);
    add_inserted_branch(prev_joining_branch);
    add_inserted_branch(Branch(next_joining_branch.lower_node, next_added_branch.upper_node));
    add_inserted_branch(Branch(next_added_branch.upper_node, next_joining_branch.upper_node));
    if (cut_node != nullptr) {
        if (next_added_branch != Branch()) {
            add_deleted_branch(Branch(next_added_branch.lower_node, cut_node));
        }
        if (prev_added_branch != Branch()) {
            add_inserted_branch(Branch(prev_added_branch.lower_node, cut_node));
        }
    }
     */
    simplify_branches();
    if (pos == 0) {
        return;
    }
    assert(deleted_branches.size() == inserted_branches.size());
    find_nodes();
    // when joining the source branch, depending on whether it joins above the start time or below, determine the new source branch
    if (prev_joining_branch == source_branch) {
        float t = prev_added_branch.upper_node->time;
        if (t > start_time) {
            source_branch = Branch(source_branch.lower_node, prev_added_branch.upper_node);
        } else {
            source_branch = Branch(prev_added_branch.upper_node, source_branch.upper_node);
        }
    } else {
        source_branch = search_lower_node(source_branch.lower_node);
    }
    find_target_branch();
    find_recomb_info();
    for (Branch b : inserted_branches) {
        assert(b.upper_node->time > b.lower_node->time);
    }
    assert(start_time < inserted_node->time);
    assert(merging_branch != Branch());
    assert(deleted_branches.count(source_branch) > 0);
    assert(deleted_branches.size() == 3 or deleted_branches.size() == 4);
}

Branch Recombination::next_added_branch(Branch prev_joining_branch, Branch prev_added_branch, Node *base_node) {
    Node *prev_node = prev_added_branch.upper_node;
    Node *next_node = nullptr;
    if (!affect(prev_joining_branch)) {
        next_node = prev_node;
    }
    if (prev_joining_branch == source_branch and prev_node->time > start_time) {
        next_node = deleted_node;
    } else {
        next_node = prev_node;
    }
    return Branch(base_node, next_node);
}

// private methods:

void Recombination::break_front(Branch next_removed_branch, Branch next_split_branch, Node *cut_node) {
    add_deleted_branch(next_removed_branch);
    add_deleted_branch(Branch(next_split_branch.lower_node, next_removed_branch.upper_node));
    add_deleted_branch(Branch(next_removed_branch.upper_node, next_split_branch.upper_node));
    add_inserted_branch(next_split_branch);
    add_inserted_branch(Branch(next_removed_branch.lower_node, cut_node));
    simplify_branches();
}

void Recombination::break_end(Branch prev_removed_branch, Branch prev_split_branch, Node *cut_node) {
    add_inserted_branch(prev_removed_branch);
    add_inserted_branch(Branch(prev_split_branch.lower_node, prev_removed_branch.upper_node));
    add_inserted_branch(Branch(prev_removed_branch.upper_node, prev_split_branch.upper_node));
    add_deleted_branch(prev_split_branch);
    add_deleted_branch(Branch(prev_removed_branch.lower_node, cut_node));
    simplify_branches();
}

void Recombination::break_front(Branch next_removed_branch, Branch next_split_branch) {
    add_deleted_branch(next_removed_branch);
    add_deleted_branch(Branch(next_split_branch.lower_node, next_removed_branch.upper_node));
    add_deleted_branch(Branch(next_removed_branch.upper_node, next_split_branch.upper_node));
    add_inserted_branch(next_split_branch);
    simplify_branches();
}

void Recombination::break_end(Branch prev_removed_branch, Branch prev_split_branch) {
    add_inserted_branch(prev_removed_branch);
    add_inserted_branch(Branch(prev_split_branch.lower_node, prev_removed_branch.upper_node));
    add_inserted_branch(Branch(prev_removed_branch.upper_node, prev_split_branch.upper_node));
    add_deleted_branch(prev_split_branch);
    simplify_branches();
}

void Recombination::fix_front(Branch next_added_branch, Branch next_joining_branch, Node *cut_node) {
    add_inserted_branch(next_added_branch);
    add_inserted_branch(Branch(next_joining_branch.lower_node, next_added_branch.upper_node));
    add_inserted_branch(Branch(next_added_branch.upper_node, next_joining_branch.upper_node));
    add_deleted_branch(next_joining_branch);
    add_deleted_branch(Branch(next_added_branch.lower_node, cut_node));
    simplify_branches();
    if (deleted_branches.size() == 0) {
        return;
    }
    source_branch = search_lower_node(source_branch.lower_node);
    find_nodes();
    find_target_branch();
    find_recomb_info();
}

void Recombination::fix_end(Branch prev_added_branch, Branch prev_joining_branch, Node *cut_node) {
    if (deleted_branches.size() == inserted_branches.size() and deleted_branches.size() == 0) {
        return;
    }
    add_deleted_branch(prev_added_branch);
    add_deleted_branch(Branch(prev_joining_branch.lower_node, prev_added_branch.upper_node));
    add_deleted_branch(Branch(prev_added_branch.upper_node, prev_joining_branch.upper_node));
    add_inserted_branch(prev_joining_branch);
    add_inserted_branch(Branch(prev_added_branch.lower_node, cut_node));
    simplify_branches();
    if (deleted_branches.size() == 0) {
        return;
    }
    source_branch = search_lower_node(source_branch.lower_node);
    find_nodes();
    find_target_branch();
    find_recomb_info();
}

void Recombination::simplify_branches() {
    /*
    set<Branch> simplified_deleted_branches;
    set<Branch> simplified_inserted_branches;
    for (Branch b : deleted_branches) {
        if (inserted_branches.count(b) == 0) {
            simplified_deleted_branches.insert(b);
        }
    }
    for (Branch b : inserted_branches) {
        if (deleted_branches.count(b) == 0) {
            simplified_inserted_branches.insert(b);
        }
    }
    deleted_branches = simplified_deleted_branches;
    inserted_branches = simplified_inserted_branches;
     */
    for (auto it = deleted_branches.begin(); it != deleted_branches.end();) {
        if (inserted_branches.count(*it) > 0) {
            inserted_branches.erase(*it);
            it = deleted_branches.erase(it);
        } else {
            ++it;
        }
    }
    // assert(deleted_branches.size() == inserted_branches.size() or deleted_branches.size() == 0);
    // assert(deleted_branches.size() == 3 or deleted_branches.size() == 4);
}

void Recombination::add_deleted_branch(Branch b) {
    if (b.upper_node != nullptr and b.lower_node != nullptr) {
        deleted_branches.insert(b);
    }
}

void Recombination::add_inserted_branch(Branch b) {
    if (b.upper_node != nullptr and b.lower_node != nullptr) {
        inserted_branches.insert(b);
    }
}

void Recombination::find_nodes() {
    // find deleted and inserted nodes by simply comparing nodes
    set<Node*> prev_nodes = {};
    set<Node*> next_nodes = {};
    for (Branch b : deleted_branches) {
        // prev_nodes.insert(b.lower_node);
        prev_nodes.insert(b.upper_node);
    }
    for (Branch b : inserted_branches) {
        // next_nodes.insert(b.lower_node);
        next_nodes.insert(b.upper_node);
    }
    for (Node* n : prev_nodes) {
        if (next_nodes.count(n) == 0) {
            deleted_node = n;
        }
    }
    for (Node* n : next_nodes) {
        if (prev_nodes.count(n) == 0) {
            inserted_node = n;
        }
    }
}

void Recombination::find_target_branch() {
    Branch lower_branch;
    Branch upper_branch;
    for (Branch b : deleted_branches) {
        if (b.lower_node->time >= inserted_node->time or b.upper_node->time <= inserted_node->time or b == source_branch) {
            continue;
        }
        lower_branch = Branch(b.lower_node, inserted_node);
        upper_branch = Branch(inserted_node, b.upper_node);
        if (inserted_branches.count(lower_branch) > 0 and inserted_branches.count(upper_branch) > 0) {
            target_branch = b;
            return;
        }
    }
    for (Branch b : deleted_branches) {
        if (b.lower_node->time >= inserted_node->time or b.upper_node->time <= inserted_node->time or b == source_branch) {
            continue;
        }
        lower_branch = Branch(b.lower_node, inserted_node);
        upper_branch = Branch(inserted_node, b.upper_node);
        if (inserted_branches.count(lower_branch) > 0 or inserted_branches.count(upper_branch) > 0) {
            target_branch = b;
            return;
        }
    }
    target_branch = Branch();
    return;
}

void Recombination::find_recomb_info() {
    if (pos == 0 or pos == INT_MAX) { // no need to process the pseudo terminal recombinations
        return;
    }
    Node *l = nullptr;
    Node *u = nullptr;
    // find merging branch by looking for deleted node in deleted branches
    for (Branch b : deleted_branches) {
        if (b == source_branch) {
            continue;
        }
        if (b.upper_node == deleted_node) {
            if (b == target_branch) {
                l = inserted_node;
            } else {
                l = b.lower_node;
            }
        } else if (b.lower_node == deleted_node) {
            if (b == target_branch) {
                u = inserted_node;
            } else {
                u = b.upper_node;
            }
        }
    }
    merging_branch = Branch(l, u);
    recombined_branch = Branch(source_branch.lower_node, inserted_node); // recombined branch is source lower node to inserted node;
    source_sister_branch = search_upper_node(deleted_node); // find sister branch of source branch in a naive way
    source_parent_branch = search_lower_node(deleted_node); // find parent branch of source branch in a naive way
    // find transfer branches
    Branch candidate_lower_transfer = Branch(target_branch.lower_node, inserted_node);
    if (create(candidate_lower_transfer)) {
        lower_transfer_branch = candidate_lower_transfer;
    } else {
        lower_transfer_branch = merging_branch;
    }
    Branch candidate_upper_transfer = Branch(inserted_node, target_branch.upper_node);
    if (create(candidate_upper_transfer)) {
        upper_transfer_branch = candidate_upper_transfer;
    } else {
        upper_transfer_branch = merging_branch;
    }
}

Branch Recombination::search_upper_node(Node *n) {
    Branch branch = Branch();
    for (Branch b : deleted_branches) {
        if (b != source_branch and b.upper_node == n) {
            branch = b;
        }
    }
    return branch;
}

Branch Recombination::search_lower_node(Node *n) {
    Branch branch = Branch();
    for (Branch b : deleted_branches) {
        if (b.lower_node == n) {
            branch = b;
        }
    }
    return branch;
}
