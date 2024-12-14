//file to test static sequential LCA

//file to test LCA with

//header copied from RC.cpp
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
#include "../examples/samplesort.h"
#include<stdio.h>

//the functions here
#include "static_seqn_LCA.h"
#include "RC.h" //for generate random tree

#include "VanillaLCA.h"

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

void test_lca(int n, int NUM_TRIALS, int seed) {
    //generate random parent tree
    parlay::sequence<int> parent_tree = give_example_tree();
    n = 22;

    //from the parent tree, get a child_tree
    parlay::sequence<parlay::sequence<int>> child_tree(n,parlay::sequence<int>());
    partree_to_childtree(parent_tree,child_tree);

    print_parent_tree(parent_tree,"parent tree");
    print_child_tree(child_tree,"child tree");

    //head needs to store all inlabels, must be slightly bigger? come back to TOD2
    parlay::sequence<int> head(2*n+1);
    parlay::sequence<LCAnode<int>> augmented_vertices(n);
    for (int i = 0; i < n; i++) {
        augmented_vertices[i].id=i; //id is index
    }

    preprocess(parent_tree,child_tree,0,augmented_vertices,head);

    std::cout << "printing augment" << std::endl;
    std::cout << "id,\t inl,\t ino,\t preo,\t#ri0,\tlvl,\tascd,\t size" << std::endl;
    for (int i =0 ; i < augmented_vertices.size(); i++) {
        augmented_vertices[i].print(); 
    }

    std::mt19937 rand_gen2 = get_rand_gen(seed+1); //to not give same randomness that generated the tree (lol should pass around same rand gen better) TOD2*

    std::uniform_int_distribution<int> dis(0,n-1);
    for (int  i  = 0; i < 20; i++) {
        int u = dis(rand_gen2);
        int v = dis(rand_gen2);
        auto ans = query(head,parent_tree,augmented_vertices,u,v);
        auto real_ans = vanilla_lca(parent_tree,u,v,0);
       
        std::cout << u << " " << v << " " << "LCA: " << ans << std::endl;

         if (ans != real_ans) {
            std::cout << "break abort" << std::endl;
            std::cout << "Said: " << ans << "real: " << real_ans << std::endl;
            std::exit(5);
        }
    }


}



//issue: two main functions -- move some things to headers? think about later. 
int main() {
    std::cout << "hello world!2" << std::endl;

    int n=10;
    int NUM_TRIALS=1;
    int seed = 42;
    //fixed seed for easier testing

    test_lca(n,NUM_TRIALS,seed);




}