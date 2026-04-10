import argparse
import numpy as np
import tskit
import os
import tszip

def read_long_ARG(node_files, branch_files, mutation_files, block_coordinates):
    if len(node_files) != len(branch_files):
        raise ValueError("Lengths of node_files and branch_files must be the same.")
    
    if len(node_files) != len(block_coordinates):
        raise ValueError("Lengths of node_files and coordinates must be the same.")    
    
    tables = tskit.TableCollection(sequence_length=0)
    node_table = tables.nodes
    branch_table = tables.edges
    
    time_zero_nodes_added = False
    node_num = 0
    sample_num = 0
    
    for node_file_index, (node_file, branch_file, mutation_file) in enumerate(zip(node_files, branch_files, mutation_files)):
        print(f"Processing segment {node_file_index}")
        node_time = np.loadtxt(node_file)
        node_num = node_table.num_rows - sample_num
        min_time = 0
        
        for t in node_time:
            if t == 0:
                if node_file_index == 0:  # Only add time 0 nodes from the first file
                    node_table.add_row(flags=tskit.NODE_IS_SAMPLE)
                    sample_num += 1
            else:
                assert t >= min_time 
                t = max(min_time + 1e-7, t)
                node_table.add_row(time=t)
                min_time = t

        if node_file_index == 0:
            time_zero_nodes_added = True
        
        edge_span = np.loadtxt(branch_file)
        edge_span = edge_span[edge_span[:, 2] >= 0, :]
        
        length = max(edge_span[:, 1])
        tables.sequence_length = length + block_coordinates[node_file_index]

        parent_indices = np.array(edge_span[:, 2], dtype=np.int32)
        child_indices = np.array(edge_span[:, 3], dtype=np.int32)
        
        parent_indices[parent_indices >= sample_num] += node_num
        child_indices[child_indices >= sample_num] += node_num
        
        branch_table.append_columns(
            left=edge_span[:, 0] + block_coordinates[node_file_index],
            right=edge_span[:, 1] + block_coordinates[node_file_index],
            parent=parent_indices,
            child=child_indices
        )
        mutations = np.loadtxt(mutation_file)
        mut_num = mutations.shape[0]
        mut_pos = 0
        for i in range(mut_num):
            if mutations[i, 0] != mut_pos and mutations[i, 0] < length:
                tables.sites.add_row(position=mutations[i, 0] + block_coordinates[node_file_index], ancestral_state='0')
                mut_pos = mutations[i, 0]
            site_id = tables.sites.num_rows - 1
            mut_node = int(mutations[i, 1])
            if (mut_node < sample_num):
                tables.mutations.add_row(site=site_id, node=int(mutations[i, 1]), derived_state=str(int(mutations[i, 3]))) 
            else:
                tables.mutations.add_row(site=site_id, node=int(mutations[i, 1]) + node_num, derived_state=str(int(mutations[i, 3])))    
    
    tables.sort()
    ts = tables.tree_sequence()
    
    return ts

def load_file_lists(file_list_path):
    node_files = []
    branch_files = []
    mutation_files = []
    block_coordinates = []

    with open(file_list_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                raise ValueError(f"Invalid line: {line}")
            node_files.append(parts[0])
            branch_files.append(parts[1])
            mutation_files.append(parts[2])
            block_coordinates.append(float(parts[3]))  # or int if desired

    return node_files, branch_files, mutation_files, block_coordinates

def sort_nodes_by_time(ts):
    tables = ts.dump_tables()
    times = tables.nodes.time
    sort_order = np.argsort(times, kind='stable')  # from most recent to most ancient (smaller to larger time)
    
    # Remap all references
    node_map = np.full(ts.num_nodes, tskit.NULL, dtype=int)
    for new_id, old_id in enumerate(sort_order):
        node_map[old_id] = new_id
    
    # Reorder nodes
    tables.nodes.set_columns(
        flags=tables.nodes.flags[sort_order],
        time=tables.nodes.time[sort_order],
        population=tables.nodes.population[sort_order],
        individual=tables.nodes.individual[sort_order],
    )
    
    # Remap edges
    edges = tables.edges
    edges.set_columns(
        left=edges.left,
        right=edges.right,
        parent=node_map[edges.parent].astype(np.int32),
        child=node_map[edges.child].astype(np.int32),
    )

    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    ts = tables.tree_sequence()     
    return tables.tree_sequence()

def write_output_ts(ts, output):
    print(f"Save to {output}")
    ts.dump(output)

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Generate tskit format for a long ARG.")
    
    # Add arguments with prefixes
    parser = argparse.ArgumentParser(description="Generate tskit format for a long ARG using file list.")
    parser.add_argument("--file_table", required=True, help="Sub file table")
    parser.add_argument("--output", required=True, help="Output file name")       
 
    args = parser.parse_args()

    # Generate file lists
    node_files, branch_files, mutation_files, block_coordinates = load_file_lists(args.file_table)
    # Apply the function
    output_ts_filename = args.output
    ts = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates)
    sorted_ts = sort_nodes_by_time(ts)
    write_output_ts(sorted_ts, args.output)    

if __name__  == "__main__":
    main()
