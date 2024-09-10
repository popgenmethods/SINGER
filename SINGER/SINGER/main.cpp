//
//  main.cpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#include <iostream>
#include "Test.hpp"

int main(int argc, const char * argv[]) {
    bool fast = false;
    bool resume = false;
    bool debug = false;
    double r = -1, m = -1, Ne = -1;
    int num_iters = 0;
    int spacing = 1;
    double start_pos = -1, end_pos = -1;
    string input_filename = "", output_prefix = "";
    string recomb_map_filename = "", mut_map_filename = "";
    double penalty = 0.01;
    double polar = 0.5;
    double epsilon_hmm = 0.1;
    double epsilon_psmc = 0.05;
    int seed = 42;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-fast") {
            if (i + 1 < argc && argv[i+1][0] != '-') {
                cerr << "Error: -fast flag doesn't take any value. " << endl;
                exit(1);
            }
            fast = true;
        }
        else if (arg == "-resume") {
            if (i + 1 < argc && argv[i+1][0] != '-') {
                cerr << "Error: -resume flag doesn't take any value. " << endl;
                exit(1);
            }
            resume = true;
        }
        else if (arg == "-debug") {
            if (i + 1 < argc && argv[i+1][0] != '-') {
                cerr << "Error: -debug flag doesn't take any value. " << endl;
                exit(1);
            }
            debug = true;
        }
        else if (arg == "-Ne") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -Ne flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                Ne = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -Ne flag expects a number. " << endl;
                exit(1);
            }
            Ne = 2*Ne;
        }
        else if (arg == "-r") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -r flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                r = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -r flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-m") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -m flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                m = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -r flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-penalty") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -penalty flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                penalty = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -penalty flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-polar") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -polar flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                polar = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -polar flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-hmm_epsilon") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -hmm_epsilon flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                epsilon_hmm = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -hmm_epsilon flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-psmc_bins") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -psmc_epsilon flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                epsilon_psmc = 1.0/stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -p flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-start") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -start flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                start_pos = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -start flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-end") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -end flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                end_pos = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -end flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-input") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -input flag cannot be empty. " << endl;
                exit(1);
            }
            input_filename = argv[++i];
        }
        else if (arg == "-output") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -output flag cannot be empty. " << endl;
                exit(1);
            }
            output_prefix = argv[++i];
        }
        else if (arg == "-recomb_map") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -recomb_map flag cannot be empty. " << endl;
                exit(1);
            }
            recomb_map_filename = argv[++i];
        }
        else if (arg == "-mut_map") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -mut_map flag cannot be empty. " << endl;
                exit(1);
            }
            mut_map_filename = argv[++i];
        }
        else if (arg == "-n") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -n flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                num_iters = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -n flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-thin") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -thin flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                spacing = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -thin flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-seed") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -seed flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                seed = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -seed flag expects a number. " << endl;
                exit(1);
            }
        }
        else {
            cerr << "Error: Unknown flag. " << arg << endl;
            exit(1);
        }
    }
    if (r < 0) {
        cerr << "-r flag missing or invalid value. " << endl;
        exit(1);
    }
    if (m < 0) {
        cerr << "-m flag missing or invalid value. " << endl;
        exit(1);
    }
    if (Ne < 0) {
        cerr << "-Ne flag missing or invalid value. " << endl;
        exit(1);
    }
    if (input_filename.size() == 0) {
        cerr << "-input flag missing or invalid value. " << endl;
        exit(1);
    }
    if (output_prefix.size() == 0) {
        cerr << "-output flag missing or invalid value. " << endl;
        exit(1);
    }
    if (num_iters < 0) {
        cerr << "-num_iters flag is invalid. " << endl;
        exit(1);
    }
    if (spacing < 1) {
        cerr << "-thin flag is invalid. " << endl;
        exit(1);
    }
    Sampler sampler;
    if (r > 0 and m > 0) {
        Sampler sampler = Sampler(Ne, r, m);
    } else {
        Rate_map recomb_map = Rate_map();
        recomb_map.load_map(recomb_map_filename);
        Rate_map mut_map = Rate_map();
        mut_map.load_map(mut_map_filename);
        Sampler sampler = Sampler(Ne, recomb_map, mut_map);
    }
    sampler.penalty = penalty;
    sampler.polar = polar;
    sampler.set_precision(epsilon_hmm, epsilon_psmc);
    sampler.set_input_file_prefix(input_filename);
    sampler.set_output_file_prefix(output_prefix);
    sampler.fast_mode = fast;
    sampler.random_seed = seed;
    sampler.start = start_pos;
    sampler.end = end_pos;
    if (resume) {
        sampler.sequence_length = end_pos - start_pos;
        sampler.resume_internal_sample(num_iters, spacing);
        return 0;
    } else if (debug) {
        sampler.sequence_length = end_pos - start_pos;
        sampler.debug_resume_internal_sample(num_iters, spacing);
        return 0;
    }
    sampler.load_vcf(input_filename, start_pos, end_pos);
    sampler.iterative_start();
    sampler.internal_sample(num_iters, spacing);
    return 0;
}
