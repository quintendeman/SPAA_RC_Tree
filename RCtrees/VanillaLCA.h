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


void parse_input(int argc, char* argv[], int& n, int& NUM_TRIALS,int& seed, int& NUM_TREES) {
    for (int i = 1; i < argc; i++) {
        std::string arg=argv[i];
        if (arg=="-n" && i+1 < argc) {
            n=std::stoi(argv[i+1]);

        }
        if (arg=="-trials" && i+1 < argc) {
            NUM_TRIALS=std::stoi(argv[i+1]);
        }
        if (arg=="-seed" && i+1 < argc) {
            seed=std::stoi(argv[i+1]);
        }
        if (arg=="-trees" && i+1 < argc) {
            NUM_TREES=std::stoi(argv[i+1]);
        }

    }
}

template<typename T>
void print_parent_tree(parlay::sequence<T>& tree, std::string message) {
    std::cout << message << std::endl;
    for (int j = 0 ; j < tree.size(); j++) {
         printf("(%d,%d) ", j , tree[j]);

    }
    std::cout << std::endl;

}

template<typename T>
void print_child_tree(parlay::sequence<parlay::sequence<T>>& child_tree, std::string message) {
    std::cout << message << std::endl;
    for (int i = 0 ; i < child_tree.size(); i++) {
        std::cout << i << ":";
        for (int j = 0; j < child_tree[i].size(); j++) {
            std::cout << child_tree[i][j] << " ";
        }
        std::cout << std::endl;
    }

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

//get arbitrary root from base root
//r is the true root of the tree
//rprime is the root around which we wish to orient our answer
template<typename T>
T unrooted_lca(parlay::sequence<T>& parent_tree, T u, T v, T r, T rprime) {
    auto lca_uv = vanilla_lca(parent_tree,u,v,r);
    auto lca_urp = vanilla_lca(parent_tree,u,rprime,r);
    auto lca_vrp = vanilla_lca(parent_tree,v,rprime,r);
    //return the lca that is different (which also happens to be deeper in the original tree)
    if (lca_uv == lca_urp) {
        return lca_vrp;
    }
    else if (lca_urp == lca_vrp) {
        return lca_uv;
    }
    else {
        return lca_urp;
    }
}

//note that we are NOT passing child_tree by reference, to get a new copy as the undirected graph
template<typename T>
parlay::sequence<parlay::sequence<T>> make_undirected_tree(parlay::sequence<parlay::sequence<T>> child_tree, parlay::sequence<T>& parent_tree) {
    for (int i = 0; i < child_tree.size(); i++) {
        child_tree[i].push_back(parent_tree[i]); //just push back the parent
    }
    return child_tree;

}

//given an unrooted tree, make it into a parent tree rooted at rprime
//pass by reference breaking for some reason?! //continue here TODO*
template<typename T>
void make_directed_tree(parlay::sequence<parlay::sequence<T>>& unrooted_tree, parlay::sequence<T>& parent_tree, T rprime) {

    std::deque<T> stack;
    parent_tree[rprime]=rprime;
    stack.push_back(rprime); 

    //int count = 0;
    while (stack.size() > 0) {
        // count += 1;
        // if (count > 50) exit(700); //for printing
        
        auto s = stack.back();
        stack.pop_back();

        // print_parent_tree(parent_tree,"intermediate directed tree");
        // std::cout << "s: " << s << std::endl;
        
        //the parent of each other edge in s is s
        for (int i = 0; i < unrooted_tree[s].size(); i++) {
            T child = unrooted_tree[s][i];
            //issue: The old root has itself as an edge
            if (child != parent_tree[s] &&  s != child) {//don't count the back edge //&&  TODO* print out child
                parent_tree[child]=s;
                stack.push_back(child);
            }
        
        }
    }

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

//level ancestors structure

//O(d) depth, where d is height of tree
//O(nd) work
//la = level ancestors
parlay::sequence<parlay::sequence<T>> preprocess_la(parlay::sequence<T>& parent_tree,parlay::sequence<LCAnode<T>>& av, T root) {
    parlay::sequence<parlay::sequence<T>> table(parent_tree.size(),parlay::sequence<T>(av.level+1,-1));
    
    parlay::parallel_for(0,parent_tree.size(),[&] (size_t i) {
        int count = 0;
        T val = i;
        while (val != root) {
            table[i][count]=val;
            count += 1;
            val=parent_tree[val];
        }

    });
    return table;
}

//use the lookup table to get the i^{th} ancestor of node
T query_la(parlay::sequence<parlay::sequence<T>>& table, T node, int ith) {
    return table[node][ith];
}


#endif