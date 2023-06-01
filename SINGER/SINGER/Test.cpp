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
    a.remove_leaf(9);
    Trace_pruner trace_pruner = Trace_pruner();
    cout << get_time() << endl;
    trace_pruner.prune_arg(a);
    cout << get_time() << endl;
    trace_pruner.write_reductions(a);
    trace_pruner.write_reduction_distance(a, "/Users/yun_deng/Desktop/SINGER/arg_files/trace_reduction_distance_100.txt");
    trace_pruner.write_reduction_size("/Users/yun_deng/Desktop/SINGER/arg_files/trace_reduction_size_100.txt");
}

void test_iterative_start() {
    srand(93723823);
    set_seed(93723823);
    // set_seed(31);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(50);
    sampler.set_sequence_length(1e6);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.iterative_start();
}

void test_fast_iterative_start() {
    srand(38);
    // set_seed(38);
    set_seed(93723823);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(50);
    sampler.set_sequence_length(1e6);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.fast_iterative_start();
}

void test_terminal_sampling() {
    srand(93723823);
    set_seed(15);
    Sampler sampler = Sampler(2e4, 1.9e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(50);
    sampler.set_sequence_length(1e6);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/sample_ts");
    sampler.iterative_start();
    sampler.terminal_sample(50);
}

void test_internal_sampling() {
    srand(93723823);
    // set_seed(93723823);
    set_seed(15);
    Sampler sampler = Sampler(2e4, 1.9e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(50);
    sampler.set_sequence_length(1e6);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/sample_ts");
    sampler.iterative_start();
    // sampler.recombination_climb(300, 1);
    sampler.internal_sample(50, 1);
}

void test_fast_internal_sampling() {
    srand(93723823);
    set_seed(93723823);
    // set_seed(38);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(50);
    sampler.set_sequence_length(1e6);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/sample_ts");
    sampler.fast_iterative_start();
    sampler.fast_internal_sample(50, 1);
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

void test_succint_bsp() {
    srand(93723823);
    // set_seed(93723823);
    set_seed(15);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_num_samples(3);
    sampler.set_sequence_length(1e6);
    sampler.set_input_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_hap0");
    sampler.iterative_start();
    Node_ptr n = sampler.build_node(3, 0.0);
    sampler.arg.add_sample(n);
    BSP_smc bsp = BSP_smc();
    succint_BSP sbsp = succint_BSP();
    bsp.reserve_memory(100000);
    sbsp.reserve_memory(100000);
    bsp.set_cutoff(0.01);
    sbsp.set_cutoff(0.01);
    bsp.set_emission(sampler.eh);
    sbsp.set_emission(sampler.eh);
    bsp.start(sampler.arg.start_tree.branches, 0.005);
    sbsp.start(sampler.arg.start_tree.branches, 0.005);
    auto recomb_it = sampler.arg.recombinations.upper_bound(0);
    auto mut_it = sampler.arg.mutation_sites.lower_bound(0);
    auto query_it = sampler.arg.removed_branches.begin();
    vector<float> mutations;
    set<float> mut_set = {};
    set<Branch> deletions = {};
    set<Branch> insertions = {};
    Node_ptr query_node = nullptr;
    for (int i = 0; i < 100000; i++) {
        if (sampler.arg.coordinates[i] == query_it->first) {
            query_node = query_it->second.lower_node;
            query_it++;
        }
        if (sampler.arg.coordinates[i] == recomb_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            bsp.transfer(r);
            sbsp.transfer(r);
        } else if (sampler.arg.coordinates[i] != 0) {
            bsp.forward(sampler.arg.rhos[i - 1]);
            sbsp.forward(sampler.arg.rhos[i - 1]);
        }
        vector<float> curr_bsp_probs = bsp.forward_probs.back();
        vector<float> curr_sbsp_probs = sbsp.forward_probs.back();
        sort(curr_bsp_probs.begin(), curr_bsp_probs.end());
        sort(curr_sbsp_probs.begin(), curr_sbsp_probs.end());
        assert(curr_bsp_probs.size() == curr_sbsp_probs.size());
        mut_set = {};
        while (*mut_it < sampler.arg.coordinates[i + 1]) {
            mut_set.insert(*mut_it);
            mut_it++;
        }
        if (mut_set.size() > 0) {
            bsp.mut_emit(sampler.arg.thetas[i], sampler.arg.coordinates[i + 1] - sampler.arg.coordinates[i], mut_set, query_node);
            sbsp.mut_emit(sampler.arg.thetas[i], sampler.arg.coordinates[i + 1] - sampler.arg.coordinates[i], mut_set, query_node);
        } else {
            bsp.null_emit(sampler.arg.thetas[i], query_node);
            sbsp.null_emit(sampler.arg.thetas[i], query_node);
        }
        vector<float> bsp_null_probs = bsp.null_emit_probs;
        vector<float> sbsp_null_probs = sbsp.null_emit_probs;
        sort(bsp_null_probs.begin(), bsp_null_probs.end());
        sort(sbsp_null_probs.begin(), sbsp_null_probs.end());
        assert(bsp_null_probs == sbsp_null_probs);
    }
}
