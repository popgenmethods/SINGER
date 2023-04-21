//
//  Test.cpp
//  SINGER
//
//  Created by Yun Deng on 3/17/23.
//

#include "Test.hpp"

void test_read_arg() {
    ARG a = ARG(2e4, 1e6);
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_branches.txt");
    a.discretize(10);
    for (Node *n : a.sample_nodes) {
        int index = n->index;
        string mutation_file = "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_sample_" + to_string(index) + ".txt";
        n->read_mutation(mutation_file);
        a.add_sample(n);
    }
    a.impute_nodes(0, 1e6);
    a.map_mutations(0, 1e6);
    cout << "Number of incompatibilities: " << a.count_incompatibility() <<  endl;
}

void test_pruner() {
    Sampler sampler = Sampler(1, 1, 1);
    srand(23491256);
    ARG a = ARG(2e4, 1e6);
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_continuous_ts_10_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_continuous_ts_10_branches.txt");
    a.discretize(100);
    a.smc_sample_recombinations();
    for (Node *n : a.sample_nodes) {
        int index = n->index;
        string mutation_file = "/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_continuous_sample_10_" + to_string(index) + ".txt";
        n->read_mutation(mutation_file);
        a.add_sample(n);
    }
    a.impute_nodes(0, 1e6);
    a.map_mutations(0, 1e6);
    a.remove_leaf(9);
    Parsimony_pruner parsimony_pruner = Parsimony_pruner();
    cout << sampler.get_time() << endl;
    parsimony_pruner.prune_arg(a);
    cout << sampler.get_time() << endl;
    // parsimony_pruner.write_reduction_distance(a, "/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_reduction_distance_100.txt");
    // parsimony_pruner.write_reduction_size("/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_reduction_size_100.txt");
}

void test_pruner_efficiency() {
    srand(93723823);
    Sampler sampler = Sampler(2e4, 2e-9, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(1);
    sampler.set_sequence_length(1e7);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/low_rho_smc_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/low_rho_smc_sample0");
    sampler.set_log_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/check");
    sampler.iterative_start();
    Node *n = sampler.build_node(1, 0.0);
    sampler.arg.add_sample(n);
    Parsimony_pruner pp = Parsimony_pruner();
    cout << sampler.get_time() << endl;
    pp.prune_arg(sampler.arg);
    cout << sampler.get_time() << endl;
}

void test_iterative_start() {
    srand(93723823);
    Sampler sampler = Sampler(2e4, 2e-9, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(8);
    sampler.set_sequence_length(1e7);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/low_rho_smc_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/low_rho_smc_sample0");
    sampler.set_log_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/check");
    sampler.iterative_start();
}

void test_fast_iterative_start() {
    srand(93723823);
    Sampler sampler = Sampler(2e4, 2e-9, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(8);
    sampler.set_sequence_length(1e7);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/low_rho_smc_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/low_rho_smc_sample0");
    sampler.set_log_file_prefix("/Users/yun_deng/Desktop/conditional-coalescent/arg_files/check");
    sampler.fast_iterative_start();
}
