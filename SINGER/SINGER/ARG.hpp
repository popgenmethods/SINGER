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
#include "RSP_smc.hpp"
#include "Reconstruction.hpp"
#include "Fitch_reconstruction.hpp"

class ARG {
    
public:
    
    float Ne = 1;
    Node *root = new Node(numeric_limits<float>::infinity());
    Node *cut_node = nullptr;
    float cut_time = -1;
    set<float> mutation_sites;
    map<float, Recombination> recombinations;
    int bin_num = 0;
    float sequence_length = 0;
    vector<float> coordinates = {};
    set<Node *, compare_node> sample_nodes = {};
    set<Node *, compare_node> node_set = {};
    map<float, Node*> base_nodes = {};
    
    ARG(float N, float l);
    
    ~ARG();
    
    void discretize(float s);
    
    void init_arg(Node *n);
    
    void add_sample(Node *n);
    
    void add_node(Node *n);
    
    Node *add_new_node(float t);
    
    Tree get_tree_at(float x);
    
    map<float, pair<Branch, Node *>> remove(tuple<float, Branch, float> cut_point);
    
    void add(map<float, pair<Branch, Node*>> joining_points);
    
    void smc_sample_recombinations();
    
    // void smcprime_sample_recombinations();
    
    int count_incompatibility();
    
    void write(string node_file, string branch_file, string recomb_file);
    
    void read(string node_file, string branch_file);
    
    void read(string node_file, string branch_file, string recomb_file);
    
    float get_arg_length();
    
    float get_arg_length(float x, float y);
    
    tuple<float, Branch, float> sample_internal_cut();
    
    tuple<float, Branch, float> sample_terminal_cut();
    
    void impute_nodes(float x, float y);
    
    void clear_memory();
    
    void clear_memory(map<int, pair<Branch, Node *>> joining_points);
    
    void check_mapping();
    
    void clear_base_nodes();
    
    float smc_prior_likelihood(float r);
    
    // float smcprime_prior_likelihood();
    
    float data_likelihood(float m);
    
    float smc_likelihood(float r, float m);
    
    // float smcprime_likelihood();
    
    set<int> get_check_points();
    
    bool check_disjoint_nodes(float x, float y);
    
    // private:
    
    float random();
    
    void new_recombination(float pos, Branch prev_added_branch, Branch prev_joining_branch, Branch next_added_branch, Branch next_joining_branch);
    
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
