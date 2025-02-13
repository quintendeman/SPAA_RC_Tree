//implementation of LCA that doesn't use RC trees, as a baseline
#ifndef VANILLA_LCA_H
#define VANILLA_LCA_H

#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "random_trees.h"

//give tree featured in Vishkin & Schieber
parlay::sequence<int> give_example_tree() {
    parlay::sequence<int> parent_tree(22);

    parent_tree[0]=0;
    parent_tree[1]=0;
    parent_tree[2]=0;
    parent_tree[3]=0;
    parent_tree[4]=1;
    parent_tree[5]=2;
    parent_tree[6]=3;
    parent_tree[7]=4;
    parent_tree[8]=4;
    parent_tree[9]=4;
    parent_tree[10]=6;
    parent_tree[11]=7;
    parent_tree[12]=7;
    parent_tree[13]=8;
    parent_tree[14]=10;
    parent_tree[15]=10;
    parent_tree[16]=10;
    parent_tree[17]=10;
    parent_tree[18]=13;
    parent_tree[19]=16;
    parent_tree[20]=19;
    parent_tree[21]=19;
    return parent_tree;

}

//a tree with bugs (at one point)
parlay::sequence<int> give_example_tree2() {
    parlay::sequence<int> parent_tree(7);
    parent_tree[1]=1;
    parent_tree[3]=1;
    parent_tree[6]=3;
    parent_tree[2]=3;
    parent_tree[4]=6;
    parent_tree[0]=2;
    parent_tree[5]=0;
    return parent_tree;

}


//a tree with bugs
//chain of length 9, permuted
parlay::sequence<int> give_example_tree3() {
    parlay::sequence<int> parent_tree(9);
    parent_tree[0]=3;
    parent_tree[1]=5;
    parent_tree[2]=1;
    parent_tree[3]=7;
    parent_tree[4]=4; //root
    parent_tree[5]=8;
    parent_tree[6]=2;
    parent_tree[7]=6;
    parent_tree[8]=4;
    return parent_tree;

}


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


//given nodes u and v in a tree rooted at r, find the LCA. Note that this does not use RC trees or anything fancy, this is using standard techniques. 
//if u and v not in same tree, return -1
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
    u_path_up.push_back(r);
    T vprime = v;
    while (vprime != r) {
        v_path_up.push_back(vprime);
        vprime = tree[vprime];
    }
    v_path_up.push_back(r);
    T lca = -1;
    
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


//given 3 LCAs, LCA(u,v) LCA(u,r) LCA(v,r), find LCA(u,v,r)
template<typename T>
T logic_lca(T lca_uv, T lca_urp, T lca_vrp) {
    if (lca_uv == lca_urp) {
        return lca_vrp;
    }
    else if (lca_urp == lca_vrp) {
        return lca_uv;
    }
    else if (lca_uv==lca_vrp) {
        return lca_urp;
    }
    else {
        std::cout << "error abort 3 lca vals are different" << std::endl;
        exit(7003);
    }

}

//get arbitrary root from base root ok
//r is the true root of the tree
//rprime is the root around which we wish to orient our answer
//this version only for a single connected tree*
template<typename T>
T unrooted_lca(parlay::sequence<T>& parent_tree, T u, T v, T r, T rprime) {
    auto lca_uv = vanilla_lca(parent_tree,u,v,r);
    auto lca_urp = vanilla_lca(parent_tree,u,rprime,r);
    auto lca_vrp = vanilla_lca(parent_tree,v,rprime,r);

    // std::cout << "uv lca " << lca_uv << std::endl;
    // std::cout << "ur'  lca " << lca_urp << std::endl;
    // std::cout << "vr' lca " << lca_vrp << std::endl;

    //return the lca that is different (which also happens to be deeper in the original tree)
    return logic_lca(lca_uv,lca_urp,lca_vrp);
}

template<typename T>
T unrooted_forest_lca(parlay::sequence<T>& parent_forest, T u, T v, T rprime) {
    T uroot = get_root(parent_forest,u);
    T vroot = get_root(parent_forest,v);
    T rprime_root = get_root(parent_forest,rprime);
    //std::cout << "uroot " << uroot << std::endl;
    
    if (uroot != vroot || uroot != rprime_root) {
        return -1; //disconnected
    }
    //if connected, can run unrooted lca as normal
    return unrooted_lca(parent_forest,u,v,uroot,rprime);

}


//given a truly unrooted tree (undirected graph no cycles), find the LCA directly, not with the above fixed root LCA trick
//do so by directing all of the edges, then running the standard LCA (somewhat time expensive)
template<typename T>
T true_unrooted_lca(parlay::sequence<parlay::sequence<T>>& unrooted_tree,T u, T v, T rprime) {
    parlay::sequence<T> parent_tree(unrooted_tree.size(),-1);
    make_directed_tree(unrooted_tree,parent_tree,rprime);
    // std::cout << "finished direct tree make" << std::endl;
    // std::cout << "r': " << rprime << std::endl;
    // print_parent_tree(parent_tree,"directed tree");
    return vanilla_lca(parent_tree,u,v,rprime);
  
}


#endif