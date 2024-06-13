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
#include "Rate_map.hpp"

class ARG {
    
public:
    
    double Ne = 1;
    Node_ptr root = new_node(numeric_limits<double>::infinity());
    Node_ptr cut_node = nullptr;
    double cut_time = 0;
    set<double> mutation_sites = {};
    map<double, set<Branch>> mutation_branches = {};
    map<double, Recombination> recombinations = {};
    int bin_num = 0;
    double sequence_length = 0;
    double bin_size = 0;
    vector<double> coordinates = {};
    vector<double> rhos = {};
    vector<double> thetas = {};
    set<Node_ptr, compare_node> sample_nodes = {};
    set<Node_ptr, compare_node> node_set = {};
    map<double, Branch> joining_branches = {};
    map<double, Branch> removed_branches = {};
    map<double, Tree> tree_map = {};
    
    double start = 0;
    double end = 0;
    double cut_pos = 0;
    Tree cut_tree;
    Tree start_tree;
    Tree end_tree;
    
    ARG();
    
    ARG(double N, double l);
    
    ~ARG();
    
    void discretize(double s);
    
    int get_index(double x);
    
    void compute_rhos_thetas(double r, double m);
    
    void compute_rhos_thetas(Rate_map &rm, Rate_map &mm);
    
    void build_singleton_arg(Node_ptr n);
    
    void add_sample(Node_ptr n);
    
    void add_node(Node_ptr n);
    
    void add_new_node(double t);
    
    Tree get_tree_at(double x);
    
    Node_ptr get_query_node_at(double x);
    
    Tree modify_tree_to(double x, Tree &reference_tree, double x0);
    
    Tree internal_modify_tree_to(double x, Tree &reference_tree, double x0);

    void remove(tuple<double, Branch, double> cut_point);
    
    void remove(map<double, Branch> seed_branches);
    
    void remove_leaf(int index);
    
    double get_updated_length();
    
    void add(map<double, Branch> &new_joining_branches, map<double, Branch> &added_branches);
    
    void smc_sample_recombinations();
    
    void approx_sample_recombinations();
    
    void adjust_recombinations();
    
    int count_incompatibility();
    
    int count_flipping();
    
    void read_coordinates(string filename);
    
    void write_coordinates(string filename);
    
    void write(string node_file, string branch_file);
    
    void write(string node_file, string branch_file, string recomb_file);
    
    void write(string node_file, string branch_file, string recomb_file, string mutation_file);
    
    void read(string node_file, string branch_file);
    
    void read(string node_file, string branch_file, string recomb_file);
    
    void read(string node_file, string branch_file, string recomb_file, string mut_file);
    
    double get_arg_length();
    
    double get_arg_length(double x, double y);
    
    double get_arg_length(map<double, Branch> &new_joining_branches, map<double, Branch> &new_added_branches);
    
    tuple<double, Branch, double> sample_internal_cut();
    
    tuple<double, Branch, double> sample_terminal_cut();
    
    tuple<double, Branch, double> sample_recombination_cut();
    
    tuple<double, Branch, double> sample_mutation_cut();
    
    void impute_nodes(double x, double y);
    
    void impute(map<double, Branch> &new_joining_branches, map<double, Branch> &added_branches);
    
    void map_mutations(double x, double y);
    
    void remap_mutations();
    
    void map_mutation(double x, Branch joining_branch, Branch added_branch);
    
    void map_mutation(Tree tree, double x);

    void check_mapping();
    
    int num_unmapped();
    
    void check_incompatibility();
    
    void clear_remove_info();
    
    double smc_prior_likelihood(double r);
    
    double data_likelihood(double m);
    
    double smc_likelihood(double r, double m);
    
    set<double> get_check_points();
    
    bool check_disjoint_nodes(double x, double y);
    
    // private:
    
    double random();
    
    void new_recombination(double pos, Branch prev_added_branch, Branch prev_joining_branch, Branch next_added_branch, Branch next_joining_branch);
    
    void remove_empty_recombinations();
    
    int count_incompatibility(Tree tree, double x);
    
    void create_node_set();
    
    void write_nodes(string filename);
    
    void write_branches(string filename);
    
    void write_recombs(string filename);
    
    void write_mutations(string filename);
    
    void read_nodes(string filename);
    
    void read_branches(string filename);
    
    void read_recombs(string filename);
    
    void read_muts(string filename);
    
    // void normalize(double t, Distribution &d);
    
    // void normalize();
    
};

bool compare_edge(const tuple<int, int, double, double>& edge1, const tuple<int, int, double, double>& edge2);


#endif /* ARG_hpp */
