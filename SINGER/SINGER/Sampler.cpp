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

void Sampler::load_vcf(string vcf_file) {
    ifstream file(vcf_file);
    string line;
    int num_individuals = 0;
    int prev_pos = -1;
    vector<Node_ptr> nodes = {};
    int valid_mutation = 0;
    int removed_mutation = 0;
    while (getline(file, line)) {
        if (line.substr(0, 6) == "#CHROM") {
            std::istringstream iss(line);
            std::vector<std::string> fields;
            std::string field;
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
        valid_mutation += 1;
        int freq = 0;
        while (iss >> genotype) {
            if (nodes.size() < 2*individual_index + 2) {
                Node_ptr left_node = new_node(0);
                Node_ptr right_node = new_node(0);
                left_node->set_index(2*individual_index);
                right_node->set_index(2*individual_index + 1);
                nodes.push_back(left_node);
                nodes.push_back(right_node);
                sample_nodes.insert(left_node);
                sample_nodes.insert(right_node);
            }
            if (genotype[0] == '1') {
                freq += 1;
                Node_ptr node = nodes[2*individual_index];
                node->add_mutation(pos);
                carriers[pos].insert(node);
                mutation_sets[node].insert(pos);
            }
            if (genotype[2] == '1') {
                freq += 1;
                Node_ptr node = nodes[2*individual_index + 1];
                node->add_mutation(pos);
                carriers[pos].insert(node);
                mutation_sets[node].insert(pos);
            }
            individual_index += 1;
        }
    }
    num_samples = (int) sample_nodes.size();
    ordered_sample_nodes = vector<Node_ptr>(sample_nodes.begin(), sample_nodes.end());
    cout << "valid mutations: " << valid_mutation << endl;
    cout << "removed mutations: " << removed_mutation << endl;
}

/*
void Sampler::load_vcf(string vcf_file) {
    ifstream file(vcf_file);
    string line;
    int num_individuals = 0;
    int prev_pos = -1;
    vector<Node_ptr> nodes = {};
    bool first_line = true;
    while (getline(file, line)) {
        if (line[0] == '#') {
            continue; // skip these header lines
        }
        if (first_line) {
            std::istringstream iss(line);
            std::vector<std::string> fields;
            std::string field;
            for (int i = 0; i < 9; i++) {
                iss >> field; // Read the first 9 fixed fields
            }
            while (iss >> field) {
                fields.push_back(field);
            }
            num_individuals = (int)fields.size();
            nodes.resize(2 * num_individuals);
            for (int i = 0; i < 2 * num_individuals; i++) {
                nodes[i] = new_node(0.0);
                nodes[i]->set_index(i);
                sample_nodes.insert(nodes[i]);
            }
            first_line = false;
        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        int pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        
        if (pos == prev_pos) {continue;} // skip multi-allelic sites
        if (ref.size() > 1 or alt.size() > 1) {continue;} // skip multi-allelic sites or structural variant
        
        streampos old_pos = file.tellg();
        string next_line;
        if (getline(file, next_line)) {
            istringstream next_iss(next_line);
            string next_chrom;
            int next_pos;
            next_iss >> next_chrom >> next_pos;
            if (next_pos == pos) {
                prev_pos = pos;
                continue;
            }
            file.seekg(old_pos);
        }
        int individual_index = 0;
        while (iss >> genotype) {
            if (nodes.size() < 2*individual_index + 2) {
                Node_ptr left_node = new_node(0);
                Node_ptr right_node = new_node(0);
                left_node->set_index(2*individual_index);
                right_node->set_index(2*individual_index + 1);
                nodes.push_back(left_node);
                nodes.push_back(right_node);
                sample_nodes.insert(left_node);
                sample_nodes.insert(right_node);
            }
            if (genotype[0] == '1') {
                nodes[2*individual_index]->add_mutation(pos);
                carriers[pos].insert(nodes[2*individual_index]);
            }
            if (genotype[2] == '1') {
                nodes[2*individual_index + 1]->add_mutation(pos);
                carriers[pos].insert(nodes[2*individual_index + 1]);
            }
            individual_index += 1;
        }
    }
    num_samples = (int) sample_nodes.size();
}
*/
 
/*
void Sampler::optimal_ordering() {
    set<Node_ptr, compare_node> covered_nodes = {};
    while (carriers.size() > 0) {
        cout << "Curr number of nodes: " << ordered_sample_nodes.size() << endl;
        cout << "Number of remaining mutations: " << carriers.size() << endl;
        auto it = max_element(carriers.begin(), carriers.end(), [](const auto& l, const auto& r) {return l.second.size() < r.second.size();});
        Node_ptr n = *(it->second.begin());
        for (float m : n->mutation_sites) {
            carriers.erase(m);
        }
        ordered_sample_nodes.push_back(n);
        covered_nodes.insert(n);
    }
    set_difference(sample_nodes.begin(), sample_nodes.end(), covered_nodes.begin(), covered_nodes.end(), back_inserter(ordered_sample_nodes));
    cout << "Finished ordering" << endl;
}
 */

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
    float bin_size = max(rho_unit/recomb_rate, 10.0f);
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
    build_singleton_arg();
    auto it = sample_nodes.begin();
    it++;
    while (it != sample_nodes.end()) {
        bsp_c = min(0.01, 0.05/arg.sample_nodes.size());
        Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
        Node_ptr n = *it;
        threader.thread(arg, n);
        arg.check_incompatibility();
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        it++;
    }
    arg.write("/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_recombs.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_muts.txt");
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

/*
void Sampler::fast_iterative_start() {
    build_singleton_arg();
    auto it = sample_nodes.begin();
    it++;
    while (it != sample_nodes.end()) {
        bsp_c = min(0.01, 0.05/arg.sample_nodes.size());
        Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
        Node_ptr n = *it;
        if (arg.sample_nodes.size() > 1) {
            threader.fast_thread(arg, n);
        } else {
            threader.thread(arg, n);
        }
        arg.check_incompatibility();
        arg.write("/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_recombs.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_ts_muts.txt");
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
*/

void Sampler::fast_iterative_start() {
    build_singleton_arg();
    auto it = ordered_sample_nodes.begin();
    it++;
    while (it != ordered_sample_nodes.end()) {
        bsp_c = min(0.01, 0.05/arg.sample_nodes.size());
        Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
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
    arg.write("/Users/yun_deng/Desktop/SINGER/arg_files/debug_fts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_fts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_fts_recombs.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/debug_fts_muts.txt");
    normalize();
    string node_file = output_prefix + "_fast_start_nodes_" + to_string(sample_index) + ".txt";
    string branch_file= output_prefix + "_fast_start_branches_" + to_string(sample_index) + ".txt";
    string recomb_file = output_prefix + "_fast_start_recombs_" + to_string(sample_index) + ".txt";
    string mut_file = output_prefix + "_fast_start_muts_" + to_string(sample_index) + ".txt";
    arg.write(node_file, branch_file, recomb_file, mut_file);
    string coord_file = output_prefix + "_fast_coordinates.txt";
    arg.write_coordinates(coord_file);
}

void Sampler::multiple_iterative_start(int num_iters) {
    
}

void Sampler::multiple_fast_iterative_start(int num_iters) {
    
}

void Sampler::recombination_climb(int num_iters, int spacing) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        random_seed = rand();
        srand(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
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
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
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
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
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
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
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
        Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
        tuple<int, Branch, float> cut_point = arg.sample_terminal_cut();
        threader.terminal_rethread(arg, cut_point);
        arg.clear_remove_info();
        string node_file = output_prefix + "_nodes_terminal_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_terminal_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_terminal_" + to_string(sample_index) + ".txt";
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
        Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
        tuple<float, Branch, float> cut_point = arg.sample_terminal_cut();
        threader.fast_terminal_rethread(arg, cut_point);
        updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
        arg.clear_remove_info();
        arg.check_incompatibility();
        string node_file = output_prefix + "_fast_nodes_terminal_" + to_string(i) + ".txt";
        string branch_file= output_prefix + "_fast_branches_terminal_" + to_string(i) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_terminal_" + to_string(i) + ".txt";
        arg.write(node_file, branch_file, recomb_file);
        cout << "Number of trees: " << arg.recombinations.size() << endl;
    }
}

void Sampler::internal_sample(int num_iters, int spacing) {
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        random_seed = rand();
        set_seed(random_seed);
        cout << "Random seed: " << random_seed << endl;
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_nodes_internal_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_internal_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_internal_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_internal_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        // cout << "Data likelihood: " << arg.data_likelihood(2e-8) << endl;
    }
}

void Sampler::fast_internal_sample(int num_iters, int spacing) {
    bsp_c = min(0.01, 0.05/arg.sample_nodes.size());
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(sample_index) << endl;
        float updated_length = 0;
        random_seed = rand();
        set_seed(random_seed);
        cout << "Random seed: " << random_seed << endl;
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_fast_nodes_internal_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_internal_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_internal_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_fast_muts_internal_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        // cout << "Data likelihood: " << arg.data_likelihood(2e-8) << endl;
    }
}

void Sampler::resume_internal_sample(int num_iters, int spacing, int resume_point) {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_nodes_internal_" + to_string(resume_point) + ".txt";
    string branch_file= output_prefix + "_branches_internal_" + to_string(resume_point) + ".txt";
    string recomb_file = output_prefix + "_recombs_internal_" + to_string(resume_point) + ".txt";
    string mut_file = output_prefix + "_muts_internal_" + to_string(resume_point) + ".txt";
    string coord_file = output_prefix + "_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
    sample_index = resume_point + 1;
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        random_seed = rand();
        set_seed(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        string node_file = output_prefix + "_nodes_internal_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_internal_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_internal_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_internal_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        // cout << "Data likelihood: " << arg.data_likelihood(2e-8) << endl;
    }
}

void Sampler::resume_fast_internal_sample(int num_iters, int spacing, int resume_point) {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_fast_nodes_internal_" + to_string(resume_point) + ".txt";
    string branch_file= output_prefix + "_fast_branches_internal_" + to_string(resume_point) + ".txt";
    string recomb_file = output_prefix + "_fast_recombs_internal_" + to_string(resume_point) + ".txt";
    string mut_file = output_prefix + "_fast_muts_internal_" + to_string(resume_point) + ".txt";
    string coord_file = output_prefix + "_fast_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
    arg.start_tree = arg.get_tree_at(arg.start);
    sample_index = resume_point + 1;
    bsp_c = min(0.01, 0.05/arg.sample_nodes.size());
    for (int i = sample_index; i < resume_point + num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_fast_nodes_internal_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_internal_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_internal_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_fast_muts_internal_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        // cout << "Data likelihood: " << arg.data_likelihood(2e-8) << endl;
    }
}

void Sampler::resume_internal_sample(int num_iters, int spacing, int resume_point, int seed) {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_nodes_internal_" + to_string(resume_point) + ".txt";
    string branch_file= output_prefix + "_branches_internal_" + to_string(resume_point) + ".txt";
    string recomb_file = output_prefix + "_recombs_internal_" + to_string(resume_point) + ".txt";
    string mut_file = output_prefix + "_muts_internal_" + to_string(resume_point) + ".txt";
    string coord_file = output_prefix + "_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
    arg.start = 0;
    arg.end = 0;
    arg.start_tree = arg.get_tree_at(arg.start);
    sample_index = resume_point + 1;
    for (int i = 0; i < num_iters; i++) {
        cout << get_time() << " Iteration: " << to_string(i) << endl;
        float updated_length = 0;
        random_seed = rand();
        set_seed(random_seed);
        while (updated_length < spacing*arg.sequence_length) {
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
            tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            threader.internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        string node_file = output_prefix + "_nodes_internal_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_internal_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_internal_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_internal_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        // cout << "Data likelihood: " << arg.data_likelihood(2e-8) << endl;
    }
}

void Sampler::resume_fast_internal_sample(int num_iters, int spacing, int resume_point, int seed) {
    arg = ARG(Ne, sequence_length);
    string node_file = output_prefix + "_fast_nodes_internal_" + to_string(resume_point) + ".txt";
    string branch_file= output_prefix + "_fast_branches_internal_" + to_string(resume_point) + ".txt";
    string recomb_file = output_prefix + "_fast_recombs_internal_" + to_string(resume_point) + ".txt";
    string mut_file = output_prefix + "_fast_muts_internal_" + to_string(resume_point) + ".txt";
    string coord_file = output_prefix + "_fast_coordinates.txt";
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    arg.compute_rhos_thetas(recomb_rate, mut_rate);
    arg.start = 140720;
    arg.end = 718690.5;
    arg.start_tree = arg.get_tree_at(arg.start);
    sample_index = resume_point + 1;
    bsp_c = min(0.01, 0.05/arg.sample_nodes.size());
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
            Threader_smc threader = Threader_smc(bsp_c, tsp_q, eh);
            // tuple<float, Branch, float> cut_point = arg.sample_internal_cut();
            tuple<float, Branch, float> cut_point = arg.sample_mutation_cut();
            threader.fast_internal_rethread(arg, cut_point);
            updated_length += arg.coordinates[threader.end_index] - arg.coordinates[threader.start_index];
            arg.clear_remove_info();
        }
        normalize();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_fast_nodes_internal_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_fast_branches_internal_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_fast_recombs_internal_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_fast_muts_internal_" + to_string(sample_index) + ".txt";
        arg.write(node_file, branch_file, recomb_file, mut_file);
        sample_index += 1;
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        // cout << "Data likelihood: " << arg.data_likelihood(2e-8) << endl;
    }
}

void Sampler::normalize() {
    Normalizer nm = Normalizer();
    nm.normalize(arg, mut_rate);
    // nm.normalize(arg);
    // nm.randomized_normalize(arg);
}
