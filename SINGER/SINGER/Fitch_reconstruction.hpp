//
//  Fitch_reconstruction.hpp
//  SINGER
//
//  Created by Yun Deng on 8/8/22.
//

#ifndef Fitch_reconstruction_hpp
#define Fitch_reconstruction_hpp

#include <stdio.h>
#include <map>
#include "Reconstruction.hpp"

class Fitch_reconstruction : public Reconstruction {
    
public:
    
    Fitch_reconstruction(Tree tree);
    
    void reconstruct(float pos);
    
    void update(Recombination& r);
    
private:
    
    Tree base_tree = Tree();
    set<Node_ptr> node_set = {};
    map<Node_ptr, pair<Node_ptr, Node_ptr>> children_nodes = {};
    map<Node_ptr, Node_ptr> parent_node = {};
    map<Node_ptr, float> pruning_node_states = {};
    map<Node_ptr, float> peeling_node_states = {};
    float recon_pos = 0.0f;
    
    void fill_tree_info(Tree tree);
    
    void update_tree_info(Recombination &r);
    
    void fill_node_states();
    
    void Fitch_up(Node_ptr c1, Node_ptr c2, Node_ptr p);
    
    void Fitch_down(Node_ptr u, Node_ptr l);
    
    void pruning_pass(Node_ptr n);
    
    void peeling_pass(Node_ptr n);
    
};

#endif /* Fitch_reconstruction_hpp */
