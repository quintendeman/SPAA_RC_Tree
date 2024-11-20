//implementation of LCA that doesn't use RC trees, as a baseline
#ifndef VANILLA_LCA_H
#define VANILLA_LCA_H


#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include <time.h>
#include <cmath>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../include/parlay/internal/get_time.h"

//test
//tree is in format of (i,tree[i]) gives (child, parent) edge
//return true if u is a descendant of v (v is ancestor), false otherwise
template<typename T>
bool is_descendant(parlay::sequence<T>& tree, T u, T v) {
    T uold = u;
    T unew = tree[u];
    if (u==v) return true; //needed for edge case, u=v
    while (uold != unew) {
        if (unew==v) return true;
        uold=unew;
        unew = tree[unew];
        
    }
    return false;
}

//return true if u is an ancestor of v, false otherwise
template<typename T>
bool is_ancestor(parlay::sequence<T>& tree, T u, T v) {
    return is_descendant(tree,v,u);
}

//returns true if c is a common ancestor of u and v, false otherwise
template <typename T>
bool common_ancestor(parlay::sequence<T>& tree, T u, T v, T c) {
    return is_ancestor(tree,c,u) && is_ancestor(tree,c,v);
}

//find the root of a (connected tree)
//root = where the parent is itself (by our construction)
template<typename T>
T get_root(parlay::sequence<T>& tree) {
    T u = tree[0]; //we pick an arbitrary element of the tree to start climbing up from. 
    T next = tree[u];
    while (u != next) {
        u=next;
        next=tree[next];
    }
    return u;
}

//given nodes u and v in a tree rooted at r, find the LCA. Note that this does not use RC trees or anything fancy, this is using standard techniques. 
template<typename T>
T vanilla_lca(parlay::sequence<T>& tree, T u, T v, T r) {

    //the LCA of the root and anything is the root itself
    if (u == r || v == r) return r;

    //keep track of the paths u and v take up the tree
    parlay::sequence<T> u_path_up;
    parlay::sequence<T> v_path_up;

    T uprime = u;
    while (uprime != r) {
        u_path_up.push_back(uprime);
        uprime=tree[uprime]; //gets out the parent
    }
    T vprime = v;
    while (vprime != r) {
        v_path_up.push_back(vprime);
        vprime = tree[vprime];
    }
    T lca = r;
    
    //while the ancestor is the same, go down the tree, stop when the paths diverge
    T u_counter = u_path_up.size() - 1;
    T v_counter = v_path_up.size() - 1;
    while (u_counter >= 0 && v_counter >= 0 && u_path_up[u_counter]==v_path_up[v_counter]) {
        lca = u_path_up[u_counter];
        u_counter--;
        v_counter--;

    }
    return lca;

}

#endif