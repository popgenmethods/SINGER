//
//  Threader_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 4/4/23.
//

#include "Threader_smc.hpp"

Threader_smc::Threader_smc(float c, float q) {
    cutoff = c;
    gap = q;
}

Threader_smc::~Threader_smc() {
}

void Threader_smc::thread(ARG &a, Node_ptr n) {
    cout << "Iteration: " << a.sample_nodes.size() << endl;
    cut_time = 0;
    a.cut_time = cut_time;
    a.add_sample(n);
    get_boundary(a);
    cout << get_time() << " : begin BSP" << endl;
    run_BSP(a);
    cout << "BSP avg num of states: " << bsp.avg_num_states() << endl;
    cout << get_time() << " : begin sampling branches" << endl;
    sample_joining_branches(a);
    cout << get_time() << " : begin TSP" << endl;
    run_TSP(a);
    cout << get_time() << " : begin sampling points" << endl;
    sample_joining_points(a);
    cout << get_time() << " : begin adding" << endl;
    a.add(new_joining_branches, added_branches);
    cout << get_time() << " : begin sampling recombination" << endl;
    a.approx_sample_recombinations();
    a.clear_remove_info();
    cout << get_time() << " : finish" << endl;
    cout << a.recombinations.size() << endl;
}

void Threader_smc::fast_thread(ARG &a, Node_ptr n) {
    cout << "Iteration: " << a.sample_nodes.size() << endl;
    cut_time = 0;
    a.cut_time = cut_time;
    a.add_sample(n);
    get_boundary(a);
    cout << get_time() << " : begin pruner" << endl;
    run_pruner(a);
    cout << get_time() << " : begin BSP" << endl;
    run_fast_BSP(a);
    cout << "BSP avg num of states: " << fbsp.avg_num_states() << endl;
    cout << get_time() << " : begin sampling branches" << endl;
    sample_fast_joining_branches(a);
    cout << get_time() << " : begin TSP" << endl;
    run_TSP(a);
    cout << get_time() << " : begin sampling points" << endl;
    sample_joining_points(a);
    cout << get_time() << " : begin adding" << endl;
    a.add(new_joining_branches, added_branches);
    cout << get_time() << " : begin sampling recombination" << endl;
    a.approx_sample_recombinations();
    a.clear_remove_info();
    cout << get_time() << " : finish" << endl;
    cout << a.recombinations.size() << endl;
}

void Threader_smc::internal_rethread(ARG &a, tuple<float, Branch, float> cut_point) {
    cut_time = get<2>(cut_point);
    // a.write("/Users/yun_deng/Desktop/SINGER/arg_files/full_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/full_ts_branches.txt");
    a.remove(cut_point);
    // a.write("/Users/yun_deng/Desktop/SINGER/arg_files/partial_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/partial_ts_branches.txt");
    get_boundary(a);
    set_check_points(a);
    run_BSP(a);
    // boundary_check(a);
    sample_joining_branches(a);
    run_TSP(a);
    sample_joining_points(a);
    float ar = acceptance_ratio(a);
    float q = random();
    if (q < ar) {
        a.add(new_joining_branches, added_branches);
    } else {
        a.add(a.joining_branches, a.removed_branches);
    }
    a.approx_sample_recombinations();
    a.clear_remove_info();
}


void Threader_smc::terminal_rethread(ARG &a, tuple<float, Branch, float> cut_point) {
    cut_time = get<2>(cut_point);
    a.remove(cut_point);
    get_boundary(a);
    set_check_points(a);
    run_BSP(a);
    // boundary_check(a);
    cout << "BSP avg num states: " << bsp.avg_num_states() << endl;
    sample_joining_branches(a);
    run_TSP(a);
    sample_joining_points(a);
    a.add(new_joining_branches, added_branches);
    a.smc_sample_recombinations();
    a.clear_remove_info();
}

void Threader_smc::fast_internal_rethread(ARG &a, tuple<float, Branch, float> cut_point) {
    cut_time = get<2>(cut_point);
    // a.write("/Users/yun_deng/Desktop/SINGER/arg_files/full_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/full_ts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/full_ts_recombs.txt");
    a.remove(cut_point);
    // a.write("/Users/yun_deng/Desktop/SINGER/arg_files/partial_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/partial_ts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/partial_ts_recombs.txt");
    get_boundary(a);
    set_check_points(a);
    run_pruner(a);
    run_fast_BSP(a);
    // boundary_check(a);
    sample_fast_joining_branches(a);
    run_TSP(a);
    sample_joining_points(a);
    float ar = acceptance_ratio(a);
    float q = random();
    if (q < ar) {
        a.add(new_joining_branches, added_branches);
    } else {
        a.add(a.joining_branches, a.removed_branches);
    }
    a.approx_sample_recombinations();
    a.clear_remove_info();
    // a.write("/Users/yun_deng/Desktop/SINGER/arg_files/full_ts_nodes.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/full_ts_branches.txt", "/Users/yun_deng/Desktop/SINGER/arg_files/full_ts_recombs.txt");
}

void Threader_smc::fast_terminal_rethread(ARG &a, tuple<float, Branch, float> cut_point) {
    cut_time = get<2>(cut_point);
    a.remove(cut_point);
    get_boundary(a);
    set_check_points(a);
    run_pruner(a);
    run_fast_BSP(a);
    // boundary_check(a);
    sample_fast_joining_branches(a);
    run_TSP(a);
    sample_joining_points(a);
    a.add(new_joining_branches, added_branches);
    a.smc_sample_recombinations();
    a.clear_remove_info();
}

void Threader_smc::get_boundary(ARG &a) {
    start = a.start;
    end = a.end;
    start_index = a.get_index(start);
    end_index = a.get_index(end);
}

void Threader_smc::set_check_points(ARG &a) {
    set<float> check_points = a.get_check_points();
    pruner.set_check_points(check_points);
    bsp.set_check_points(check_points);
    fbsp.set_check_points(check_points);
    tsp.set_check_points(check_points);
}

/*
void Threader_smc::boundary_check(ARG &a) {
    float end_pos = a.end;
    if (bsp.check_points.count(end_pos) > 0) {
        Recombination &r = a.recombinations[end_pos];
        bsp.sanity_check(r);
        return;
    }
    if (fbsp.check_points.count(end_pos) > 0) {
        Recombination &r = a.recombinations[end_pos];
        fbsp.sanity_check(r);
        return;
    }
    if (tsp.check_points.count(end_pos) > 0) {
        Recombination &r = a.recombinations[end_pos];
        tsp.sanity_check(r);
    }
}
 */

void Threader_smc::run_pruner(ARG &a) {
    pruner.prune_arg(a);
}

void Threader_smc::run_BSP(ARG &a) {
    bsp.reserve_memory(end_index - start_index);
    bsp.set_cutoff(cutoff);
    bsp.set_emission(pe);
    bsp.start(a.start_tree, cut_time);
    auto recomb_it = a.recombinations.upper_bound(start);
    auto mut_it = a.mutation_sites.lower_bound(start);
    auto query_it = a.removed_branches.begin();
    vector<float> mutations;
    set<float> mut_set = {};
    set<Branch> deletions = {};
    set<Branch> insertions = {};
    Node_ptr query_node = nullptr;
    for (int i = start_index; i < end_index; i++) {
        if (a.coordinates[i] == query_it->first) {
            query_node = query_it->second.lower_node;
            query_it++;
        }
        if (a.coordinates[i] == recomb_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            bsp.transfer(r);
        } else if (a.coordinates[i] != start) {
            bsp.forward(a.rhos[i - 1]);
        }
        mut_set = {};
        while (*mut_it < a.coordinates[i + 1]) {
            mut_set.insert(*mut_it);
            mut_it++;
        }
        if (mut_set.size() > 0) {
            bsp.mut_emit(a.thetas[i], a.coordinates[i + 1] - a.coordinates[i], mut_set, query_node);
        } else {
            bsp.null_emit(a.thetas[i], query_node);
        }
    }
    if (bsp.check_points.count(end) > 0) {
        Recombination &r = a.recombinations[end];
        bsp.sanity_check(r);
    }
}


void Threader_smc::run_fast_BSP(ARG &a) {
    fbsp.reserve_memory(end_index - start_index);
    fbsp.set_cutoff(cutoff);
    fbsp.set_emission(pe);
    set<Interval_info> start_intervals = pruner.insertions.begin()->second;
    fbsp.start(a.start_tree, start_intervals, cut_time);
    auto recomb_it = a.recombinations.upper_bound(start);
    auto mut_it = a.mutation_sites.lower_bound(start);
    auto query_it = a.removed_branches.begin();
    auto delete_it = pruner.deletions.upper_bound(start);
    auto insert_it = pruner.insertions.upper_bound(start);
    vector<float> mutations;
    set<float> mut_set = {};
    Node_ptr query_node = nullptr;
    for (int i = start_index; i < end_index; i++) {
        if (a.coordinates[i] == query_it->first) {
            query_node = query_it->second.lower_node;
            query_it++;
        }
        while (delete_it->first <= a.coordinates[i]) {
            fbsp.update_states(delete_it->second, insert_it->second);
            delete_it++;
            insert_it++;
        }
        if (a.coordinates[i] == recomb_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            fbsp.transfer(r);
        } else if (a.coordinates[i] != start) {
            fbsp.forward(a.rhos[i - 1]);
        }
        mut_set = {};
        while (*mut_it < a.coordinates[i+1]) {
            mut_set.insert(*mut_it);
            mut_it++;
        }
        if (mut_set.size() > 0) {
            fbsp.mut_emit(a.thetas[i], a.coordinates[i+1] - a.coordinates[i], mut_set, query_node);
        } else {
            fbsp.null_emit(a.thetas[i], query_node);
        }
    }
    if (fbsp.check_points.count(end) > 0) {
        Recombination &r = a.recombinations[end];
        fbsp.sanity_check(r);
    }
}

void Threader_smc::run_TSP(ARG &a) {
    tsp.reserve_memory(end_index - start_index);
    tsp.set_gap(gap);
    tsp.set_emission(be);
    Branch start_branch = new_joining_branches.begin()->second;
    tsp.start(start_branch, cut_time);
    auto recomb_it = a.recombinations.upper_bound(start);
    auto join_it = new_joining_branches.upper_bound(start);
    auto mut_it = a.mutation_sites.lower_bound(start);
    auto query_it = a.removed_branches.lower_bound(start);
    Branch prev_branch = start_branch;
    Branch next_branch = start_branch;
    Node_ptr query_node = nullptr;
    set<float> mut_set = {};
    for (int i = start_index; i < end_index; i++) {
        if (a.coordinates[i] == query_it->first) {
            query_node = query_it->second.lower_node;
            query_it++;
        }
        if (a.coordinates[i] == join_it->first) {
            next_branch = join_it->second;
            join_it++;
        }
        if (a.coordinates[i] == recomb_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            tsp.transfer(r, prev_branch, next_branch);
            prev_branch = next_branch;
        } else if (prev_branch != next_branch) {
            tsp.recombine(prev_branch, next_branch);
            prev_branch = next_branch;
        } else if (a.coordinates[i] != start) {
            float rho = a.rhos[i];
            tsp.forward(rho);
        }
        mut_set.clear();
        while (*mut_it < a.coordinates[i+1]) {
            mut_set.insert(*mut_it);
            mut_it++;
        }
        if (mut_set.size() > 0) {
            tsp.mut_emit(a.thetas[i], a.coordinates[i+1] - a.coordinates[i], mut_set, query_node);
        } else {
            tsp.null_emit(a.thetas[i], query_node);
        }
    }
    if (tsp.check_points.count(end) > 0) {
        Recombination &r = a.recombinations[end];
        tsp.sanity_check(r);
    }
}

void Threader_smc::sample_joining_branches(ARG &a) {
    new_joining_branches = bsp.sample_joining_branches(start_index, a.coordinates);
}

void Threader_smc::sample_fast_joining_branches(ARG &a) {
    new_joining_branches = fbsp.sample_joining_branches(start_index, a.coordinates);
}

void Threader_smc::sample_joining_points(ARG &a) {
    map<float, Node_ptr> added_nodes = tsp.sample_joining_nodes(start_index, a.coordinates);
    auto add_it = added_nodes.begin();
    auto end_it = added_nodes.end();
    Node_ptr query_node = nullptr;
    Node_ptr added_node = nullptr;
    float x;
    while (add_it != end_it) {
        x = add_it->first;
        added_node = add_it->second;
        query_node = a.get_query_node_at(x);
        added_branches[x] = Branch(query_node, added_node);
        add_it++;
    }
}

float Threader_smc::acceptance_ratio(ARG &a) {
    float cut_height = a.cut_tree.parents.rbegin()->first->time;
    float old_height = cut_height;
    float new_height = cut_height;
    auto old_join_it = a.joining_branches.upper_bound(a.cut_pos);
    old_join_it--;
    auto new_join_it = new_joining_branches.upper_bound(a.cut_pos);
    new_join_it--;
    auto old_add_it = a.removed_branches.upper_bound(a.cut_pos);
    old_add_it--;
    auto new_add_it = added_branches.upper_bound(a.cut_pos);
    new_add_it--;
    if (old_join_it->second.upper_node == a.root) {
        old_height = old_add_it->second.upper_node->time;
    }
    if (new_join_it->second.upper_node == a.root) {
        new_height = new_add_it->second.upper_node->time;
    }
    return old_height/new_height;
}

float Threader_smc::random() {
    return uniform_random();
}

vector<float> Threader_smc::expected_diff(float m) {
    vector<float> diff(3);
    auto join_it = new_joining_branches.begin();
    auto add_it = added_branches.begin();
    float span;
    while (add_it != prev(added_branches.end())) {
        span = next(add_it)->first - add_it->first;
        diff[0] += m*span*(add_it->second.upper_node->time - join_it->second.lower_node->time);
        if (join_it->second.upper_node->index != -1) {
            diff[1] += m*span*(join_it->second.upper_node->time - add_it->second.upper_node->time);
        }
        diff[2] += m*span*(add_it->second.upper_node->time - add_it->second.lower_node->time);
        add_it++;
        if (add_it->first == next(join_it)->first) {
            join_it++;
        }
    }
    return diff;
}

vector<float> Threader_smc::observed_diff(ARG &a) {
    vector<float> diff(3);
    auto join_it = new_joining_branches.begin();
    auto add_it = added_branches.begin();
    auto mut_it = a.mutation_sites.begin();
    float sl, su, sm, s0;
    while (add_it != prev(added_branches.end())) {
        Branch &joining_branch = join_it->second;
        Branch &added_branch = add_it->second;
        while (*mut_it < next(add_it)->first) {
            float m = *mut_it;
            sl = joining_branch.lower_node->get_state(m);
            su = joining_branch.upper_node->get_state(m);
            sm = added_branch.upper_node->get_state(m);
            s0 = added_branch.lower_node->get_state(m);
            diff[0] += abs(sm - sl);
            if (join_it->second.upper_node->index != -1) {
                diff[1] += abs(su - sm);
            }
            diff[2] += abs(sm - s0);
            mut_it++;
        }
        add_it++;
        if (next(join_it)->first == add_it->first) {
            join_it++;
        }
    }
    return diff;
}
