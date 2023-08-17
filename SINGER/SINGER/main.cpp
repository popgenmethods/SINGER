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
    float r = -1, m = -1, Ne = -1;
    int num_iters = 0;
    int spacing = 1;
    float start_pos = -1, end_pos = -1;
    string input_filename = "", output_prefix = "";
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-fast") {
            if (i + 1 < argc && argv[i+1][0] != '-') {
                cerr << "Error: -fast flag doesn't take any value. " << endl;
                return 1;
            }
            fast = true;
        }
        else if (arg == "-Ne") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -Ne flag cannot be empty. " << endl;
                return 1;
            }
            try {
                Ne = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -Ne flag expects a number. " << endl;
                return 1;
            }
            Ne = 2*Ne;
        }
        else if (arg == "-r") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -r flag cannot be empty. " << endl;
                return 1;
            }
            try {
                r = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -r flag expects a number. " << endl;
                return 1;
            }
        }
        else if (arg == "-m") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -m flag cannot be empty. " << endl;
                return 1;
            }
            try {
                m = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -r flag expects a number. " << endl;
                return 1;
            }
        }
        else if (arg == "-start") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -start flag cannot be empty. " << endl;
                return 1;
            }
            try {
                start_pos = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -start flag expects a number. " << endl;
                return 1;
            }
        }
        else if (arg == "-end") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -end flag cannot be empty. " << endl;
                return 1;
            }
            try {
                end_pos = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -end flag expects a number. " << endl;
                return 1;
            }
        }
        else if (arg == "-input") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -input flag cannot be empty. " << endl;
            }
            input_filename = argv[++i];
        }
        else if (arg == "-output") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -output flag cannot be empty. " << endl;
            }
            output_prefix = argv[++i];
        }
        else if (arg == "-n") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -n flag cannot be empty. " << endl;
            }
            try {
                num_iters = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -n flag expects a number. " << endl;
                return 1;
            }
        }
        else if (arg == "-thinning") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -thinning flag cannot be empty. " << endl;
            }
            try {
                spacing = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -thinning flag expects a number. " << endl;
                return 1;
            }
        }
        else {
            cerr << "Error: Unknown flag. " << arg << endl;
            return 1;
        }
    }
    if (r < 0) {
        cerr << "-r flag missing or invalid value. " << endl;
    }
    if (m < 0) {
        cerr << "-m flag missing or invalid value. " << endl;
    }
    if (Ne < 0) {
        cerr << "-Ne flag missing or invalid value. " << endl;
    }
    if (input_filename.size() == 0) {
        cerr << "-input flag missing or invalid value. " << endl;
    }
    if (output_prefix.size() == 0) {
        cerr << "-output flag missing or invalid value. " << endl;
    }
    if (num_iters < 0) {
        cerr << "-num_iters flag is invalid. " << endl;
    }
    if (spacing < 1) {
        cerr << "-thinning flag is invalid. " << endl;
    }
    Sampler sampler = Sampler(Ne, r, m);
    sampler.set_precision(0.01, 0.05);
    sampler.set_output_file_prefix(output_prefix);
    sampler.load_vcf(input_filename, start_pos, end_pos);
    if (fast) {
        sampler.fast_iterative_start();
    } else {
        sampler.iterative_start();
    }
    if (fast) {
        sampler.internal_sample(num_iters, spacing);
    } else {
        sampler.internal_sample(num_iters, spacing);
    }
    return 0;
}
