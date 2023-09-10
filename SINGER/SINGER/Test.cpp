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
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0.vcf", 0, 1e6);
    sampler.iterative_start();
}

void test_fast_iterative_start() {
    set_seed(93723823);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0.vcf", 0, 1e6);
    sampler.fast_iterative_start();
}

void test_fast_larger_iterative_start() {
    set_seed(842651);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_1");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_200_1.vcf", 0, 1e6);
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
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0.vcf", 0, 1e6);
    sampler.iterative_start();
    sampler.terminal_sample(200);
}

void test_fast_terminal_sampling() {
    set_seed(6328);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_100_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_100_0.vcf", 0, 1e6);
    sampler.fast_iterative_start();
    sampler.fast_terminal_sample(1000);
}

void test_internal_sampling() {
    set_seed(15);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/approx_smc_50_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0", 0, 1e6);
    sampler.iterative_start();
    sampler.internal_sample(2000, 1);
}

void test_fast_internal_sampling() {
    // srand(73);
    set_seed(7254);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/approx_smc_50_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0.vcf", 0, 1e6);
    sampler.fast_iterative_start();
    sampler.fast_internal_sample(50, 1);
}

void test_fast_larger_internal_sampling() {
    set_seed(842651);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.01);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_500_0");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_500_0.vcf", 0, 1e6);
    sampler.fast_iterative_start();
    sampler.fast_internal_sample(1000, 1);
}

void test_optimal_ordering() {
    srand(93723823);
    set_seed(93723823);
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
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0.vcf", 0, 1e6);
    sampler.optimal_ordering();
}

void test_no_recomb() {
    set_seed(8);
    Sampler sampler = Sampler(2e4, 0, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e3);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/no_recomb");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/no_recomb.vcf", 0, 1e6);
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
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/no_mut_8.vcf", 0, 1e6);
    sampler.iterative_start();
    // sampler.internal_sample(500, 1);
}

/*
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
 */

void test_normalizer() {
    ARG a = ARG(8e4, 1e6);
    // a.read("/Users/yun_deng/Desktop/SINGER/arg_files/african_16_start_nodes_0.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/african_16_start_branches_0.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/african_16_start_recombs_0.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/african_16_start_muts_0.txt");
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/debug_fts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_fts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_fts_recombs.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_fts_muts.txt");
    for (int i = 0; i < 10; i++) {
        Normalizer nm = Normalizer();
        nm.normalize(a, 8e4*1.25e-8);
        a.write("/Users/yun_deng/Desktop/SINGER/arg_files/normalized_african_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/normalized_african_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/normalized_african_recombs.txt");
    }
}

void test_read_muts() {
    ARG a = ARG(2e4, 1e6);
    a.read("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0_nodes_internal_1999.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0_branches_internal_1999.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0_recombs_internal_1999.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0_muts_internal_1999.txt");
    a.read_coordinates("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0_coordinates.txt");
    a.check_incompatibility();
}

void test_resume() {
    set_seed(15);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0");
    // sampler.resume_internal_sample(1000, 1, 1999);
}

void test_resume_fast_larger_internal_sampling() {
    set_seed(15);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_500_0");
    // sampler.resume_fast_internal_sample(100, 1, 3914);
}

void test_fast_resume() {
    set_seed(15);
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/approx_smc_50_0");
    // sampler.resume_fast_internal_sample(1000, 1, 1032);
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
    Threader_smc threader = Threader_smc(0.01, 0.05);
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

void test_debug_resume() {
    Sampler sampler = Sampler(2e4, 2e-8, 2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/smc_50_0");
    // sampler.resume_fast_internal_sample(500, 1, 1640, 1066149215);
}

void test_african_dataset() {
    set_seed(72);
    // set_seed(15);
    Sampler sampler = Sampler(7e4, 1.2e-8, 1.2e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/african_subsample_16.0_1");
    sampler.load_vcf("/Users/yun_deng/Desktop/SINGER/arg_files/16.subsampled.0_999999.vcf", 0, 1e6);
    sampler.fast_iterative_start();
    sampler.fast_internal_sample(500, 1);
    // sampler.mutation_climb(500, 1);
    // sampler.resume_internal_sample(500, 1, 0);
}

void test_start_african_sampling() {
    set_seed(72);
    Sampler sampler = Sampler(8e4, 1e-8, 1.25e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/african_16_full");
    sampler.start_fast_internal_sample(500, 1);
}

void test_resume_african_dataset() {
    Sampler sampler = Sampler(6e4, 1e-8, 1.25e-8);
    sampler.set_precision(0.01, 0.05);
    sampler.set_sequence_length(1e6);
    sampler.set_output_file_prefix("/Users/yun_deng/Desktop/SINGER/arg_files/african_16");
    // sampler.resume_fast_internal_sample(500, 1, 1079, 913090935);
}
