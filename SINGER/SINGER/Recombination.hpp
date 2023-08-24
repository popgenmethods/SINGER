//
//  Recombination.hpp
//  SINGER
//
//  Created by Yun Deng on 4/9/22.
//

#ifndef Recombination_hpp
#define Recombination_hpp

#include <stdio.h>
#include "Branch.hpp"

class Recombination {
    
public:
    
    float pos;
    Branch source_branch;
    Branch target_branch;
    Branch source_sister_branch;
    Branch source_parent_branch;
    Branch recombined_branch;
    Branch merging_branch;
    Branch lower_transfer_branch;
    Branch upper_transfer_branch;
    float start_time = -1;
    Node_ptr deleted_node;
    Node_ptr inserted_node;
    set<Branch> deleted_branches = {};
    set<Branch> inserted_branches = {};
    
    Recombination();
    
    Recombination(set<Branch> db, set<Branch> ib);
    
    void set_pos(float x);
    
    bool affect(const Branch &b);
    
    bool create(const Branch &b);
    
    void find_nodes(); // find deleted/inserted nodes;
    
    void find_target_branch(); // find target branch;
    
    void find_recomb_info(); // find other branches;
    
    Branch trace_forward(float t, Branch curr_branch);
    
    Branch trace_backward(float t, Branch curr_branch);
    
    Branch next_joining_branch(Branch removed_branch, Branch joining_branch);
    
    Branch prev_joining_branch(Branch removed_branch, Branch joining_branch);
    
    void remove(Branch prev_removed_branch, Branch next_removed_branch, Branch prev_split_branch, Branch next_split_branch, Node_ptr cut_node);
    
    void remove(Branch prev_removed_branch, Branch next_removed_branch, Branch prev_split_branch, Branch next_split_branch);
    
    void add(Branch prev_added_branch, Branch next_added_branch, Branch prev_joining_branch, Branch next_joining_branch, Node_ptr cut_node);
    
    void break_front(Branch next_removed_branch, Branch next_split_branch, Node_ptr cut_node);
    
    void break_end(Branch next_removed_branch, Branch next_split_branch, Node_ptr cut_node);
    
    void break_front(Branch next_removed_branch, Branch next_split_branch);
    
    void break_end(Branch next_removed_branch, Branch next_split_branch);
    
    void fix_front(Branch next_added_branch, Branch next_joining_branch, Node_ptr cut_node);
    
    void fix_end(Branch prev_added_branch, Branch prev_joining_branch, Node_ptr cut_node);
    
    Branch next_added_branch(Branch prev_joining_branch, Branch prev_added_branch, Node_ptr base_node);
    
// private:
    
    void simplify_branches();
    
    void add_deleted_branch(Branch b);
    
    void add_inserted_branch(Branch b);
    
    Branch search_upper_node(Node_ptr n); // search in the deleted branches a non-source branch with n as upper node
    
    Branch search_lower_node(Node_ptr n); // search in the deleted branches a non-source branch with n as lower node
    
};

#endif /* Recombination_hpp */
