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
#include "Tree.hpp"
#include "Reconstruction.hpp"

class Fitch_reconstruction : public Reconstruction {
    
public:
    
    Fitch_reconstruction(Tree tree);
    
    void reconstruct(float pos);
    
    void update(Recombination& r);
    
private:
    
    Tree base_tree = Tree();
    set<Node*> node_set = {};
    map<Node*, pair<Node*, Node*>> children_nodes = {};
    map<Node*, Node*> parent_node = {};
    map<Node*, float> pruning_node_states = {};
    map<Node*, float> peeling_node_states = {};
    float recon_pos = 0.0f;
    
    void fill_tree_info(Tree tree);
    
    void update_tree_info(Recombination &r);
    
    void fill_node_states();
    
    void Fitch_up(Node *c1, Node *c2, Node *p);
    
    void Fitch_down(Node *u, Node *l);
    
    void pruning_pass(Node *n);
    
    void peeling_pass(Node *n);
    
};

#endif /* Fitch_reconstruction_hpp */
