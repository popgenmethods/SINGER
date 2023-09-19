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

/*
void Sampler::load_vcf(string vcf_file, float start_pos, float end_pos) {
    ifstream file(vcf_file);
    string line;
    int num_individuals = 0;
    int prev_pos = -1;
    vector<Node_ptr> nodes = {};
    int valid_mutation = 0;
    int removed_mutation = 0;
    vector<float> genotypes = {};
    while (getline(file, line)) {
        if (line.substr(0, 6) == "#CHROM") {
            istringstream iss(line);
            vector<string> fields;
            string field;
            while (iss >> field) {
                fields.push_back(field);
            }
            num_individuals = (int) fields.size() - 9;
            nodes.resize(2*num_individuals);
            for (int i = 0; i < 2*num_individuals; i++) {
                nodes[i] = new_node(0.0);
                nodes[i]->set_index(i);
                sample_nodes.insert(nodes[i]);
            }
            genotypes.resize(2*num_individuals);
        } else if (line[0] == '#') {
            continue; // skip these header lines
        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        int pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        
        if (pos == prev_pos) {continue;} // skip multi-allelic sites
        if (ref.size() > 1 or alt.size() > 1) {
            removed_mutation += 1;
            continue;
        } // skip multi-allelic sites or structural variant
        
        streampos old_pos = file.tellg();
        string next_line;
        if (getline(file, next_line)) {
            istringstream next_iss(next_line);
            string next_chrom;
            int next_pos;
            next_iss >> next_chrom >> next_pos;
            if (next_pos == pos) {
                removed_mutation += 1;
                prev_pos = pos;
                continue;
            }
            file.seekg(old_pos);
        }
        int individual_index = 0;
        while (iss >> genotype) {
            if (genotype[0] == '1') {
                genotypes[2*individual_index] = 1;
            } else {
                genotypes[2*individual_index] = 0;
            }
            if (genotype[2] == '1') {
                genotypes[2*individual_index + 1] = 1;
            } else {
                genotypes[2*individual_index + 1] = 0;
            }
            individual_index += 1;
        }
        int genotype_sum = accumulate(genotypes.begin(), genotypes.end(), 0.0);
        if (genotype_sum >= 1 and genotype_sum < genotypes.size()) {
            valid_mutation += 1;
            for (int i = 0; i < genotypes.size(); i++) {
                if (genotypes[i] == 1) {
                    nodes[i]->add_mutation(pos - start_pos);
                }
            }
        }
    }
    num_samples = (int) sample_nodes.size();
    ordered_sample_nodes = vector<Node_ptr>(sample_nodes.begin(), sample_nodes.end());
    shuffle(ordered_sample_nodes.begin(), ordered_sample_nodes.end(), random_engine);
    sequence_length = end_pos - start_pos;
    cout << "valid mutations: " << valid_mutation << endl;
    cout << "removed mutations: " << removed_mutation << endl;
}
 */

void Sampler::load_vcf(string prefix, float start, float end) {
    random_engine.seed(random_seed);
    string vcf_file = prefix + ".vcf";
    string index_file = prefix + ".index";
    ifstream idx_stream(index_file);
    if (!idx_stream.is_open()) {
        cerr << "Index file not found: " + index_file << endl;
        exit(1);
    }
    string line;
    long byte_offset = -1;
    while (getline(idx_stream, line)) {
        istringstream iss(line);
        float segment_start;
        long offset;
        iss >> segment_start >> offset;
        if (segment_start == start) {
            byte_offset = offset;
            break;
        }
    }
    if (byte_offset == -1) {
        cerr << "Start position not found in index file: " + index_file << endl;
        exit(1);
    }
    ifstream vcf_stream(vcf_file, ios::binary);
    if (!vcf_stream.is_open()) {
        cerr << "VCF file not found: " + vcf_file << endl;
    }
    vcf_stream.seekg(byte_offset, ios::beg);
    int prev_pos = -1;
    vector<Node_ptr> nodes = {};
    int valid_mutation = 0;
    int removed_mutation = 0;
    vector<float> genotypes = {};
    while (getline(vcf_stream, line)) {
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        int pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        if (pos == prev_pos) {continue;} // skip multi-allelic sites
        if (pos >= end) {break;} // variant out of scope
        if (ref.size() > 1 or alt.size() > 1) {
            removed_mutation += 1;
            continue;
        } // skip multi-allelic sites or structural variant
        streampos old_pos = vcf_stream.tellg();
        string next_line;
        if (getline(vcf_stream, next_line)) {
            istringstream next_iss(next_line);
            string next_chrom;
            int next_pos;
            next_iss >> next_chrom >> next_pos;
            if (next_pos == pos) {
                removed_mutation += 1;
                prev_pos = pos;
                continue;
            }
            vcf_stream.seekg(old_pos);
        }
        int individual_index = 0;
        while (iss >> genotype) {
            if (genotypes.size() < 2*individual_index + 2) {
                genotypes.resize(2*individual_index + 2);
            }
            if (genotype[0] == '1') {
                genotypes[2*individual_index] = 1;
            } else {
                genotypes[2*individual_index] = 0;
            }
            if (genotype[2] == '1') {
                genotypes[2*individual_index + 1] = 1;
            } else {
                genotypes[2*individual_index + 1] = 0;
            }
            individual_index += 1;
        }
        if (nodes.size() == 0) {
            nodes.resize(genotypes.size());
            for (int i = 0; i < nodes.size(); i++) {
                nodes[i] = new_node(0.0);
                nodes[i]->set_index(i);
                sample_nodes.insert(nodes[i]);
            }
        } else {
            assert(nodes.size() == genotypes.size());
        }
        int genotype_sum = accumulate(genotypes.begin(), genotypes.end(), 0.0);
        if (genotype_sum >= 1 and genotype_sum < genotypes.size()) {
            valid_mutation += 1;
            for (int i = 0; i < genotypes.size(); i++) {
                if (genotypes[i] == 1) {
                    nodes[i]->add_mutation(pos - start);
                }
            }
        }
    }
    if (valid_mutation < 3) {
        cerr << "there are too few variants in this region, algorithm not run" << endl;
    }
    num_samples = (int) sample_nodes.size();
    ordered_sample_nodes = vector<Node_ptr>(sample_nodes.begin(), sample_nodes.end());
    shuffle(ordered_sample_nodes.begin(), ordered_sample_nodes.end(), random_engine);
    sequence_length = end - start;
    cout << "valid mutations: " << valid_mutation << endl;
    cout << "removed mutations: " << removed_mutation << endl;
}

void Sampler::optimal_ordering() {
    ordered_sample_nodes.clear();
    set<Node_ptr, compare_node> covered_nodes = {};
    set<float> covered_mutations = {};
    while (mutation_sets.size() > 0) {
        cout << "Curr number of nodes: " << ordered_sample_nodes.size() << endl;
        cout << "Number of covered mutations: " << covered_mutations.size() << endl;
        auto it = min_element(mutation_sets.begin(), mutation_sets.end(), [](const auto& l, const auto& r) {return l.second.size() < r.second.size();});
        Node_ptr n = it->first;
        set<float> curr_mutations = it->second;
        mutation_sets.erase(n);
        for (auto &x : mutation_sets) {
            for (float m : curr_mutations) {
                x.second.erase(m);
            }
        }
        for (float m : curr_mutations) {
            covered_mutations.insert(m);
        }
        ordered_sample_nodes.push_back(n);
        covered_nodes.insert(n);
    }
    cout << "Finished ordering" << endl;
}

Node_ptr Sampler::build_node(int index, float time) {
    Node_ptr n = new_node(time);
    n->index = index;
    string mutation_file = input_prefix + "_" + to_string(index) + ".txt";
    n->read_mutation(mutation_file);
    return n;
}

void Sampler::build_all_nodes() {
    for (int i = 0; i < num_samples; i++) {
        Node_ptr n = build_node(i, 0.0);
        sample_nodes.insert(n);
    }
}

void Sampler::build_singleton_arg() {
    // float bin_size = max(rho_unit/recomb_rate, 10.0f);
    float bin_size = rho_unit/recomb_rate;
    Node_ptr n = *ordered_sample_nodes.begin();
    arg = ARG(Ne, sequence_length);
    arg.discretize(bin_size);
    arg.build_singleton_arg(n);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
}

void Sampler::build_void_arg() {
    float bin_size = rho_unit/recomb_rate;
    arg = ARG(Ne, sequence_length);
    arg.discretize(bin_size);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
}

void Sampler::iterative_start() {
    start_log();
    build_singleton_arg();
    auto it = ordered_sample_nodes.begin();
    it++;
    while (it != ordered_sample_nodes.end()) {
        random_engine.seed(random_seed);
        Threader_smc threader = Threader_smc(bsp_c, tsp_q);
        threader.pe->penalty = penalty;
        threader.pe->ancestral_prob = polar;
        Node_ptr n = *it;
        threader.thread(arg, n);
        arg.check_incompatibility();
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        it++;
        random_seed = random_engine();
        write_iterative_start();
    }
    cout << "orignal ARG length: " << arg.get_arg_length() << endl;
    normalize();
    cout << "rescaled ARG length: " << arg.get_arg_length() << endl;
    string node_file = output_prefix + "_start_nodes_" + to_string(sample_index) + ".txt";
    string branch_file= output_prefix + "_start_branches_" + to_string(sample_index) + ".txt";
    string recomb_file = output_prefix + "_start_recombs_" + to_string(sample_index) + ".txt";
    string mut_file = output_prefix + "_start_muts_" + to_string(sample_index) + ".txt";
    arg.write(node_file, branch_file, recomb_file, mut_file);
    string coord_file = output_prefix + "_coordinates.txt";
    arg.write_coordinates(coord_file);
}

void Sampler::fast_iterative_start() {
    build_singleton_arg();
    auto it = ordered_sample_nodes.begin();
    it++;
    while (it != ordered_sample_nodes.end()) {
        Threader_smc threader = Threader_smc(bsp_c, tsp_q);
        threader.pe->penalty = penalty;
        threader.pe->ancestral_prob = polar;
        Node_ptr n = *it;
        if (arg.sample_nodes.size() > 1) {
            threader.fast_thread(arg, n);
        } else {
            threader.thread(arg, n);
        }
        if (arg.sample_nodes.size() % 5 == 4) {
            normalize();
        }
        arg.check_incompatibility();
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        it++;
    }
    normalize();
    string node_file = output_prefix + "_fast_start_nodes_" + to_string(sample_index) + ".txt";
    string branch_file= output_prefix + "_fast_start_branches_" + to_string(sample_index) + ".txt";
    string recomb_file = output_prefix + "_fast_start_recombs_" + to_string(sample_index) + ".txt";
    string mut_file = output_prefix + "_fast_start_muts_" + to_string(sample_index) + ".txt";
    arg.write(node_file, branch_file, recomb_file, mut_file);
    string coord_file = output_prefix + "_fast_coordinates.txt";
    arg.write_coordinates(coord_file);
}

void Sampler::recombination_climb(int num_iters, int spacing) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        random_seed = rand();
        srand(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_recombination_cut();
            threader.internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        arg.check_incompatibility();
        string node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
    }
}

void Sampler::mutation_climb(int num_iters, int spacing) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        random_seed = rand();
        srand(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_mutation_cut();
            threader.internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        arg.check_incompatibility();
        string node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
    }
}

void Sampler::fast_recombination_climb(int num_iters, int spacing) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(sample_index) << endl;
        float updated_length = 0;
        random_seed = rand();
        srand(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_recombination_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        arg.check_incompatibility();
        string node_file = output_prefix + "_fast_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
    }
}

void Sampler::fast_mutation_climb(int num_iters, int spacing) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(sample_index) << endl;
        float updated_length = 0;
        random_seed = rand();
        srand(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_mutation_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        arg.check_incompatibility();
        string node_file = output_prefix + "_fast_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
    }
}

void Sampler::terminal_sample(int num_iters) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        random_seed = rand();
        srand(random_seed);
        Threader_smc threader = Threader_smc(bsp_c, tsp_q);
        threader.pe->penalty = penalty;
        threader.pe->ancestral_prob = polar;
        tuple<int, Branch, float> cut_point = arg.sample_terminal_cut();
        threader.terminal_rethread(arg, cut_point);
        arg.clear_remove_info();
        string node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
        arg.check_incompatibility();
        arg.write(node_file, branch_file, recomb_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
    }
}

void Sampler::fast_terminal_sample(int num_iters) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        random_seed = rand();
        srand(random_seed);
        Threader_smc threader = Threader_smc(bsp_c, tsp_q);
        threader.pe->penalty = penalty;
        threader.pe->ancestral_prob = polar;
        tuple<float, Branch, float> cut_point = arg.sample_terminal_cut();
        threader.fast_terminal_rethread(arg, cut_point);
        updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
        arg.clear_remove_info();
        arg.check_incompatibility();
        string node_file = output_prefix + "_fast_nodes_" + to_string(i) + ".txt";
        string branch_file= output_prefix + "_fast_branches_" + to_string(i) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_" + to_string(i) + ".txt";
        arg.write(node_file, branch_file, recomb_file);
        cout << "Number of trees: " << arg.recombinations.size() << endl;
    }
}

void Sampler::internal_sample(int num_iters, int spacing) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        cout << "Random seed: " << random_seed << endl;
        random_engine.seed(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
            // write_cut(cut_point);
        }
        normalize();
        random_seed = random_engine();
        write_sample();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_" + to_string(sample_index) + ".txt";
        sample_index += 1;
        arg.write(node_file, branch_file, recomb_file, mut_file);
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
    }
}

void Sampler::fast_internal_sample(int num_iters, int spacing) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(sample_index) << endl;
        float updated_length = 0;
        random_seed = rand();
        set_seed(random_seed);
        cout << "Random seed: " << random_seed << endl;
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_fast_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_fast_muts_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        // cout << "Data likelihood: " << arg.data_likelihood(2e-8) << endl;
    }
}

void Sampler::start_fast_internal_sample(int num_iters, int spacing) {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_fast_start_nodes_0.txt";
    string branch_file= output_prefix + "_fast_start_branches_0.txt";
    string recomb_file = output_prefix + "_fast_start_recombs_0.txt";
    string mut_file = output_prefix + "_fast_start_muts_0.txt";
    string coord_file = output_prefix + "_fast_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
    arg.start_tree = arg.get_tree_at(arg.start);
    sample_index = 0;
    bsp_c = 0.1;
    for (int i = sample_index; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_fast_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_fast_muts_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
    }
}

void Sampler::resume_internal_sample(int num_iters, int spacing) {
    string log_file = output_prefix + ".log";
    read_resume_point(log_file);
    sample_index += 1;
    arg.check_incompatibility();
    cout << "Number of trees: " << arg.recombinations.size() << endl;
    cout << "Number of flippings: " << arg.count_flipping() << endl;
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(sample_index) << endl;
        float updated_length = 0;
        cout << "Random seed: " << random_seed << endl;
        random_engine.seed(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
            // write_cut(cut_point);
        }
        normalize();
        random_seed = random_engine();
        write_sample();
        arg.check_incompatibility();
        string node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
    }
}

void Sampler::resume_fast_internal_sample(int num_iters, int spacing) {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_fast_nodes_" + to_string(sample_index) + ".txt";
    string branch_file= output_prefix + "_fast_branches_" + to_string(sample_index) + ".txt";
    string recomb_file = output_prefix + "_fast_recombs_" + to_string(sample_index) + ".txt";
    string mut_file = output_prefix + "_fast_muts_" + to_string(sample_index) + ".txt";
    string coord_file = output_prefix + "_fast_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
    arg.start_tree = arg.get_tree_at(arg.start);
    sample_index += 1;
    bsp_c = 0.1;
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(sample_index) << endl;
        float updated_length = 0;
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_fast_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_fast_muts_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
    }
}

/*
void Sampler::resume_internal_sample(int num_iters, int spacing, int resume_point, int seed, float cut_pos) {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_nodes_" + to_string(resume_point) + ".txt";
    string branch_file= output_prefix + "_branches_" + to_string(resume_point) + ".txt";
    string recomb_file = output_prefix + "_recombs_" + to_string(resume_point) + ".txt";
    string mut_file = output_prefix + "_muts_" + to_string(resume_point) + ".txt";
    string coord_file = output_prefix + "_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
    arg.end = cut_pos;
    arg.end_tree = arg.get_tree_at(arg.end);
    sample_index = resume_point + 1;
    for (int i = 0; i < resume_point + num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        cout << "Random seed: " << random_seed << endl;
        random_engine.seed(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        random_seed = random_engine();
        arg.check_incompatibility();
        string node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
    }
}

void Sampler::resume_fast_internal_sample(int num_iters, int spacing, int resume_point, int seed, float cut_pos) {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_fast_nodes_" + to_string(resume_point) + ".txt";
    string branch_file= output_prefix + "_fast_branches_" + to_string(resume_point) + ".txt";
    string recomb_file = output_prefix + "_fast_recombs_" + to_string(resume_point) + ".txt";
    string mut_file = output_prefix + "_fast_muts_" + to_string(resume_point) + ".txt";
    string coord_file = output_prefix + "_fast_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
    arg.start = cut_pos;
    arg.start_tree = arg.get_tree_at(arg.start);
    sample_index = resume_point + 1;
    bsp_c = 0.1;
    for (int i = sample_index; i < resume_point + num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        if (i == resume_point + 1) {
            random_seed = seed;
            set_seed(random_seed);
        } else {
            random_seed = rand();
            set_seed(random_seed);
        }
        cout << "Random seed: " << random_seed << endl;
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q);
            threader.pe->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            tuple<float, Branch, float> cut_point = arg.sample_mutation_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_fast_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_fast_muts_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
    }
}
*/
 
void Sampler::normalize() {
    Normalizer nm = Normalizer();
    nm.normalize(arg, mut_rate);
}

void Sampler::start_log() {
    string filename = output_prefix + ".log";
    ofstream file(filename, ios::out|ios::trunc);
    if (!file) {
        cerr << "Error opening the file: " << filename << endl;
        return;
    }
    file << "Time" << "\t"
    << "Iteration:" << "\t"
    << "Threading_type" << "\t"
    << "#Recombinations" << "\t"
    << "#Mutations_not_uniquely_mapped" << "\t"
    << "Last_updated_pos" << "\t"
    << "Random_seed" << "\t"
    << "Counter" << endl;
    file.close();
}

void Sampler::write_iterative_start() {
    string filename = output_prefix + ".log";
    ofstream file(filename, ios::out|ios::app);
    if (!file) {
        cerr << "Error opening the file: " << filename << endl;
        return;
    }
    file << get_time() << "\t"
    << arg.sample_nodes.size() << "\t"
    << "initial_thread" << "\t"
    << arg.recombinations.size() - 2 << "\t"
    << arg.num_unmapped() << "\t"
    << setprecision(numeric_limits<float>::max_digits10)
    << arg.end << "\t"
    << random_seed << "\t"
    << TSP_smc::counter << endl;
}

void Sampler::write_sample() {
    string filename = output_prefix + ".log";
    ofstream file(filename, ios::out|ios::app);
    if (!file) {
        cerr << "Error opening the file: " << filename << endl;
        return;
    }
    file << get_time() << "\t"
    << sample_index << "\t"
    << "rethread" << "\t"
    << arg.recombinations.size() - 2 << "\t"
    << arg.num_unmapped() << "\t"
    << setprecision(numeric_limits<float>::max_digits10)
    << arg.end << "\t"
    << random_seed << "\t"
    << TSP_smc::counter << endl;
}

void Sampler::write_cut(tuple<float, Branch, float> cut_point) {
    string filename = output_prefix + "_cut.log";
    ofstream file(filename, ios::out|ios::app);
    if (!file) {
        cerr << "Error opening the file: " << filename << endl;
        return;
    }
    file << get<0>(cut_point) << "\t"
    << get<1>(cut_point).lower_node->index << "\t"
    << get<1>(cut_point).upper_node->index << "\t"
    << get<2>(cut_point) << "\t"
    << arg.count_incompatibility() << "\t"
    << arg.recombinations.size() << "\t"
    << arg.count_flipping() << "\t"
    << endl;
}

void Sampler::load_resume_arg() {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
    string branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
    string recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
    string mut_file = output_prefix + "_muts_" + to_string(sample_index) + ".txt";
    string coord_file = output_prefix + "_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
}

void Sampler::read_resume_point(string filename) {
    ifstream file(filename, ios::in);
    vector<string> words;
    if (!file.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        exit(1);
    }
    file.seekg(-1, ios_base::end);
    bool found = false;
    int newlines_encountered = 0;
    char ch;
    while (file.get(ch)) {
        if (ch == '\n') {
            newlines_encountered++;
            if (newlines_encountered == 2) {
                found = true;
                break;
            }
        }

        // Move one character backward
        file.seekg(-2, ios_base::cur);
    }

    if (!found) {
        cerr << "Didn't find the desired line in the file!" << endl;
        exit(1);
    }

    string last_line_with_content;
    getline(file, last_line_with_content);

    stringstream ss(last_line_with_content);
    string word;
    while (ss >> word) {
        words.push_back(word);
    }
    int log_length = (int) words.size();
    // Fill the resume point information
    TSP_smc::counter = stoi(words[log_length - 1]);
    random_seed = stoi(words[log_length - 2]);
    file.seekg(0, ios::beg);
    sample_index = stoi(words[1]);
    // sample_index = 403;
    load_resume_arg();
    arg.end = stof(words[log_length - 3]);
    arg.end_tree = arg.get_tree_at(arg.end);
}
