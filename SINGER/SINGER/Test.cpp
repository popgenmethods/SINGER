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
    for (Node_ptr n : a.sample_nodes) {
        int index = n->index;
        string mutation_file = "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_sample_" + to_string(index) + ".txt";
        n->read_mutation(mutation_file);
        a.add_sample(n);
    }
    a.impute_nodes(0, 1e6);
    a.map_mutations(0, 1e6);
    cout << "Number of incompatibilities: " << a.count_incompatibility() <<  endl;
}

void test_parsimony_pruner() {
    srand(23491256);
    ARG a = ARG(2e4, 1e6);
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_continuous_ts_10_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_continuous_ts_10_branches.txt");
    a.discretize(100);
    a.smc_sample_recombinations();
    for (Node_ptr n : a.sample_nodes) {
        int index = n->index;
        string mutation_file = "/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_continuous_sample_10_" + to_string(index) + ".txt";
        n->read_mutation(mutation_file);
        a.add_sample(n);
    }
    a.impute_nodes(0, 1e6);
    a.map_mutations(0, 1e6);
    a.remove_leaf(9);
    Parsimony_pruner parsimony_pruner = Parsimony_pruner();
    cout << get_time() << endl;
    parsimony_pruner.prune_arg(a);
    cout << get_time() << endl;
    parsimony_pruner.write_reductions(a);
    parsimony_pruner.write_reduction_distance(a, "/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_parsimony_reduction_distance_10.txt");
    parsimony_pruner.write_reduction_size("/Users/yun_deng/Desktop/SINGER/arg_files/high_mu_parsimony_reduction_size_10.txt");
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
    Node_ptr n = sampler.build_node(1, 0.0);
    sampler.arg.add_sample(n);
    Parsimony_pruner pp = Parsimony_pruner();
    cout << get_time() << endl;
    pp.prune_arg(sampler.arg);
    cout << get_time() << endl;
}

void test_trace_pruner() {
    Sampler sampler = Sampler(1, 1, 1);
    srand(23491256);
    set_seed(23491256);
    ARG a = ARG(2e4, 1e6);
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_100_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_100_branches.txt");
    a.discretize(100);
    a.smc_sample_recombinations();
    for (Node_ptr n : a.sample_nodes) {
        int index = n->index;
        string mutation_file = "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_sample_100_" + to_string(index) + ".txt";
        n->read_mutation(mutation_file);
        a.add_sample(n);
    }
    a.impute_nodes(0, 1e6);
    a.map_mutations(0, 1e6);
    a.remove_leaf(99);
    Trace_pruner trace_pruner = Trace_pruner();
    cout << get_time() << endl;
    trace_pruner.prune_arg(a);
    cout << get_time() << endl;
    trace_pruner.write_reductions(a);
    trace_pruner.write_reduction_distance(a, "/Users/yun_deng/Desktop/SINGER/arg_files/trace_reduction_distance_100.txt");
    trace_pruner.write_reduction_size("/Users/yun_deng/Desktop/SINGER/arg_files/trace_reduction_size_100.txt");
}

void test_iterative_start() {
    set_seed(93723823);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0.vcf");
    sampler.iterative_start();
}

void test_fast_iterative_start() {
    set_seed(93723823);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0.vcf");
    sampler.fast_iterative_start();
}

void test_fast_larger_iterative_start() {
    set_seed(842651);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_1");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_1.vcf");
    sampler.fast_iterative_start();
}

void test_terminal_sampling() {
    srand(93723823);
    set_seed(38);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(50);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0.vcf");
    sampler.iterative_start();
    sampler.terminal_sample(200);
}

void test_fast_terminal_sampling() {
    set_seed(6328);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_100_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_100_0.vcf");
    sampler.fast_iterative_start();
    sampler.fast_terminal_sample(1000);
}

void test_internal_sampling() {
    set_seed(93723823);
    Sampler sampler = Sampler(2e4, 1e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(100);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0.vcf");
    sampler.iterative_start();
    // sampler.recombination_climb(300, 1);
    sampler.internal_sample(500, 1);
}

void test_fast_internal_sampling() {
    set_seed(72354);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0.vcf");
    sampler.fast_iterative_start();
    sampler.fast_internal_sample(500, 1);
    // sampler.fast_recombination_climb(1000, 1);
}

void test_fast_larger_internal_sampling() {
    set_seed(842651);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.01);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_500_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_500_0.vcf");
    sampler.fast_iterative_start();
    // sampler.fast_recombination_climb(1000, 1);
    sampler.fast_internal_sample(1000, 1);
}

void test_optimal_ordering() {
    srand(93723823);
    set_seed(93723823);
    // set_seed(38);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(50);
    sampler.set_sequence_length(1e6);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/sample_ts");
    sampler.optimal_ordering();
}

void test_load_vcf() {
    set_seed(93723823);
    // set_seed(38);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0.vcf");
    sampler.optimal_ordering();
}

void test_no_recomb() {
    set_seed(8);
    Sampler sampler = Sampler(2e4, 0, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e3);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/no_recomb");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/no_recomb.vcf");
    sampler.iterative_start();
    sampler.internal_sample(1000, 1);
}

void test_no_mut() {
    set_seed(38);
    // set_seed(93723823);
    Sampler sampler = Sampler(2e4, 2e-8, 0);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/no_mut_8");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/no_mut_8.vcf");
    sampler.iterative_start();
    // sampler.internal_sample(500, 1);
}

void test_normalization() {
    set_seed(93723823);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0.vcf");
    sampler.iterative_start();
    Distribution d = Distribution();
    d.read("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_quantiles.txt");
    // sampler.arg.normalize(0.1, d);
    sampler.arg.normalize();
    sampler.arg.normalize();
    sampler.arg.normalize();
    sampler.arg.write("/Users/yun_deng/Desktop/SINGER/arg_files/normalized_ts_node.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/normalized_ts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/normalized_ts_recombs.txt");
    // sampler.internal_sample(500, 1);
}

void test_normalizer() {
    set_seed(93723823);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_0.vcf");
    sampler.fast_iterative_start();
    Normalizer nm = Normalizer();
    nm.normalize(sampler.arg);
    sampler.arg.write("/Users/yun_deng/Desktop/SINGER/arg_files/normalized_200_node.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/normalized_200_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/normalized_200_recombs.txt");
}

void test_tsp() {
    set_seed(2423);
    ARG a = ARG(2e4, 1e6);
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_10_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/continuous_ts_10_branches.txt");
    a.discretize(10);
    a.smc_sample_recombinations();
    a.remove_leaf(9);
    a.compute_rhos_thetas(4e-4, 0.0);
    shared_ptr<Binary_emission> e = make_shared<Binary_emission>();
    Threader_smc threader = Threader_smc(0.01, 0.05, e);
    Node_ptr n = a.removed_branches.begin()->second.lower_node;
    threader.end_index = (int) a.coordinates.size();
    threader.new_joining_branches = a.joining_branches;
    threader.bsp.simplify(threader.new_joining_branches);
    threader.run_TSP(a);
    threader.sample_joining_points(a);
    a.add(threader.new_joining_branches, threader.added_branches);
    a.smc_sample_recombinations();
    a.write("/Users/yun_deng/Desktop/SINGER/arg_files/new_continuous_ts_10_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/new_continuous_ts_10_branches.txt");
}
