//
//  Fitch_reconstruction.cpp
//  SINGER
//
//  Created by Yun Deng on 8/8/22.
//

#include "Fitch_reconstruction.hpp"

Fitch_reconstruction::Fitch_reconstruction(Tree tree) {
    base_tree = tree;
    fill_tree_info(base_tree);
}

void Fitch_reconstruction::reconstruct(float pos) {
    recon_pos = pos;
    pruning_node_states.clear();
    peeling_node_states.clear();
    for (Node_ptr n : node_set) {
        pruning_pass(n);
    }
    for (Node_ptr n : node_set) {
        peeling_pass(n);
    }
    for (auto x : peeling_node_states) {
        x.first->write_state(pos, x.second);
    }
}

void Fitch_reconstruction::update(Recombination &r) {
    // assert(r.deleted_branches.size() == r.inserted_branches.size());
    update_tree_info(r);
}


// private methods:

void Fitch_reconstruction::fill_tree_info(Tree tree) {
    node_set.clear();
    children_nodes.clear();
    parent_node.clear();
    for (auto &x : tree.parents) {
        Node_ptr u = x.second;
        Node_ptr l = x.first;
        node_set.insert(l);
        parent_node.insert({l, u});
        if (children_nodes.count(u) > 0) {
            children_nodes.at(u).second = l;
        } else {
            children_nodes.insert({u, {l, nullptr}});
        }
    }
    for (auto x : children_nodes) {
        assert(x.first->index < 0 or x.second.second != nullptr);
    }
}

void Fitch_reconstruction::update_tree_info(Recombination &r) {
    base_tree.forward_update(r);
    fill_tree_info(base_tree);
}

void Fitch_reconstruction::Fitch_up(Node_ptr c1, Node_ptr c2, Node_ptr p) {
    if (c1 == nullptr or c2 == nullptr) {
        return;
    }
    float s;
    float s1 = pruning_node_states.at(c1);
    float s2 = pruning_node_states.at(c2);
    if (s1 == 0.5) {
        s = s2;
    } else if (s2 == 0.5) {
        s = s1;
    } else if (s1 != s2) {
        s = 0.5;
    } else {
        s = s1;
    }
    pruning_node_states.insert({p, s});
    return;
}

void Fitch_reconstruction::Fitch_down(Node_ptr u, Node_ptr l) {
    if (u->index == -1) {
        float top_state = pruning_node_states.at(l);
        if (top_state == 0.5) {
            peeling_node_states.insert({l, 0});
        } else {
            peeling_node_states.insert({l, top_state});
        }
        return;
    }
    float sp = peeling_node_states.at(u);
    float sc = pruning_node_states.at(l);
    if (sp == 0 or sp == 1) {
        if (sc == 0.5) {
            peeling_node_states.insert({l, sp});
        } else {
            peeling_node_states.insert({l, sc});
        }
    } else {
        peeling_node_states.insert({l, sc});
    }
}

// going up

void Fitch_reconstruction::pruning_pass(Node_ptr n) {
    if (pruning_node_states.count(n) > 0) {
        return;
    }
    float s;
    if (children_nodes.count(n) == 0) {
        s = n->get_state(recon_pos);
        pruning_node_states.insert({n, s});
        return;
    }
    Node_ptr c1;
    Node_ptr c2;
    tie(c1, c2) = children_nodes.at(n);
    pruning_pass(c1);
    pruning_pass(c2);
    Fitch_up(c1, c2, n);
}

// going down
void Fitch_reconstruction::peeling_pass(Node_ptr n) {
    if (peeling_node_states.count(n) > 0) {
        return;
    }
    float s;
    if (parent_node.count(n) == 0) {
        s = pruning_node_states.at(n);
        peeling_node_states.insert({n, s});
        return;
    }
    Node_ptr p = parent_node.at(n);
    if (p->index == -1) {
        float top_state = pruning_node_states.at(n);
        if (top_state == 0.5) {
            s = 0;
        } else {
            s = top_state;
        }
        peeling_node_states.insert({n, s});
    } else if (p->index == -2) {
        s = pruning_node_states.at(n);
        peeling_node_states.insert({n, s});
    } else {
        peeling_pass(p);
        Fitch_down(p, n);
    }
    return;
}
