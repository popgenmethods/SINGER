//
//  ARG.hpp
//  SINGER
//
//  Created by Yun Deng on 4/14/22.
//

#ifndef ARG_hpp
#define ARG_hpp

#include <stdio.h>
#include <map>
#include "Recombination.hpp"
#include "Tree.hpp"
#include "RSP.hpp"
#include "Fitch_reconstruction.hpp"
#include "Reconstruction.hpp"

class ARG {
    
public:
    
    Node *root = new Node(numeric_limits<float>::infinity());
    Node *cut_node = nullptr;
    float cut_time = -1;
    map<int, Node*> base_nodes = {};
    map<int, Node *> upper_nodes = {};
    map<int, set<float>> mutation_info;
    map<int, Recombination> recombination_info;
    int bin_num;
    double rho_unit;
    vector<float> bin_sizes;
    vector<float> rhos;
    vector<float> thetas;
    set<Node *, compare_node> sample_nodes;
    set<Node *, compare_node> node_set;
    
    float sequence_length = 0;
    float bin_size = 0;
    float recomb_rate;
    float mut_rate;
    
    ARG();
    
    ARG(Node *n, float s, float l);
    
    ARG(float s, float l);
    
    ~ARG();
    
    void set_rates(float r, float m);
    
    void add_sample(Node *n);
    
    void add_node(Node *n);
    
    Node *select_node(int index);
    
    Node *add_new_node(float t);
    
    Tree get_tree_at(int x);
    
    map<int, pair<Branch, Node *>> remove(tuple<int, Branch, float> cut_point);
    
    void add(map<int, pair<Branch, Node*>> joining_points);
    
    void sample_recombinations();
    
    void check_incompatibility();
    
    int count_incompatibility();
    
    void write_mutations(string mut_file);
    
    void write(float Ne, string node_file, string branch_file, string recomb_file);
    
    void read(string node_file, string branch_file);
    
    void read(string node_file, string branch_file, string recomb_file);
    
    float get_arg_length();
    
    float get_arg_length(int x, int y);
    
    tuple<int, Branch, float> sample_internal_cut();
    
    tuple<int, Branch, float> sample_terminal_cut();
    
    tuple<int, Branch, float> get_terminal_cut(int i);
    
    void impute_nodes(int x, int y);
    
    void clear_memory();
    
    void clear_memory(map<int, pair<Branch, Node *>> joining_points);
    
    void check_mapping();
    
    void clear_upper_nodes();
    
    void clear_base_nodes();
    
    float prior_likelihood();
    
    float data_likelihood();
    
    float likelihood();
    
    int get_update_length();
    
    bool node_check();
    
    set<int> get_check_points();
    
    bool check_disjoint_nodes(int x, int y);
    
    vector<float> get_observed_distances(map<int, pair<Branch, Node *>> joining_points);
    
    vector<float> get_expected_distances(map<int, pair<Branch, Node *>> joining_points);
    
    // private:
    
    float random();
    
    void new_recombination(int pos, Branch prev_added_branch, Branch prev_joining_branch, Branch next_added_branch, Branch next_joining_branch);
    
    void remove_empty_recombinations();
    
    int count_incompatibility(Tree tree, float x);
    
    void sort_nodes();
    
    void write_nodes(string filename);
    
    void write_branches(string filename);
    
    void write_recombs(string filename);
    
    void read_nodes(string filename);
    
    void read_branches(string filename);
    
    void read_recombs(string filename);
    
};

bool compare_edge(const tuple<int, int, float, float>& edge1, const tuple<int, int, float, float>& edge2);


#endif /* ARG_hpp */
