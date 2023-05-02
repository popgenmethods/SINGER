//
//  Sampler.cpp
//  SINGER
//
//  Created by Yun Deng on 3/31/23.
//

#include "Sampler.hpp"

Sampler::Sampler(float pop_size, float r, float m) {
    Ne = pop_size;
    mut_rate = m*pop_size;
    recomb_rate = r*pop_size;
}

string Sampler::get_time() {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
    auto timer = system_clock::to_time_t(now);

    std::tm bt = *std::localtime(&timer);
    std::ostringstream oss;
    oss << "[" << std::put_time(&bt, "%H:%M:%S"); // HH:MM:SS
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count() << "]";
    return oss.str();
}

void Sampler::set_precision(float c, float q) {
    bsp_c = c;
    tsp_q = q;
}

void Sampler::set_pop_size(float n) {
    Ne = n;
}

void Sampler::set_input_file_prefix(string f) {
    input_prefix = f;
}

void Sampler::set_output_file_prefix(string f) {
    output_prefix = f;
}

void Sampler::set_log_file_prefix(string f) {
    log_prefix = f;
}

void Sampler::set_sequence_length(float x) {
    sequence_length = x;
}

void Sampler::set_num_samples(int n) {
    num_samples = n;
}

Node *Sampler::build_node(int index, float time) {
    Node *n = new Node(time);
    n->index = index;
    string mutation_file = input_prefix + "_" + to_string(index) + ".txt";
    n->read_mutation(mutation_file);
    return n;
}

void Sampler::build_singleton_arg() {
    float bin_size = rho_unit/recomb_rate;
    Node *n = build_node(0, 0.0);
    arg = ARG(Ne, sequence_length);
    arg.discretize(bin_size);
    arg.build_singleton_arg(n);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
}

void Sampler::iterative_start() {
    build_singleton_arg();
    for (int i = 1; i < num_samples; i++) {
        random_seed = rand();
        srand(random_seed);
        Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
        Node *n = build_node(i, 0.0);
        threader.thread(arg, n);
        // arg.check_incompatibility();
    }
    arg.write("/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/arg_files/debug_ts_recombs.txt");
    // arg.check_incompatibility();
}

void Sampler::fast_iterative_start() {
    build_singleton_arg();
    for (int i = 1; i < num_samples; i++) {
        random_seed = rand();
        srand(random_seed);
        Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
        Node *n = build_node(i, 0.0);
        if (arg.sample_nodes.size() > 1) {
            threader.fast_thread(arg, n);
        } else {
            threader.thread(arg, n);
        }
        arg.write("/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/arg_files/debug_ts_recombs.txt");
    }
}
