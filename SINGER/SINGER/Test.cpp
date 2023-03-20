//
//  Test.cpp
//  SINGER
//
//  Created by Yun Deng on 3/17/23.
//

#include "Test.hpp"

void test_read_arg() {
    ARG a = ARG(2e4, 1e7);
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_branches.txt");
    a.discretize(100);
    for (Node *n : a.sample_nodes) {
        int index = n->index;
        string mutation_file = "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_sample_" + to_string(index) + ".txt";
        n->read_mutation(mutation_file);
        a.add_sample(n);
    }
    a.impute_nodes(0, 1e7);
    cout << "Number of incompatibilities: " << a.count_incompatibility() <<  endl;
}

void test_pruner() {
    ARG a = ARG(2e4, 1e7);
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_branches.txt");
    a.discretize(100);
    a.smc_sample_recombinations();
    for (Node *n : a.sample_nodes) {
        int index = n->index;
        string mutation_file = "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_sample_" + to_string(index) + ".txt";
        n->read_mutation(mutation_file);
        a.add_sample(n);
    }
    a.impute_nodes(0, 1e7);
    Node *query_node = new Node(0.0);
    query_node->read_mutation("/Users/yun_deng/Desktop/SINGER/arg_files/continuous_sample_99.txt");
    query_node->set_index(99);
    map<float, Node *> base_nodes = {{0, query_node}, {INT_MAX, query_node}};
    Parsimony_pruner parsimony_pruner = Parsimony_pruner();
    parsimony_pruner.extend(a, query_node);
}
