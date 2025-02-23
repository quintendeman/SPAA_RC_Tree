//file to test static parallel LCA
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
#include "../examples/samplesort.h"
#include<stdio.h>

//the functions here
#include "RC.h" //for generate random tree
#include "static_seqn_LCA.h"

#include "static_par_LCA.h"

#include "VanillaLCA.h"
#include "random_trees.h"

void test_lca(int n, int NUM_TRIALS, int NUM_TREES, std::mt19937& gen) {


    for (int iter = 0; iter < NUM_TREES; iter++) {
    //generate random parent tree
    parlay::sequence<int> parent_tree = give_example_tree2();
    int root = 1; //set back to 0*

    //from the parent tree, get a child_tree
    parlay::sequence<parlay::sequence<int>> child_tree(n,parlay::sequence<int>());
    partree_to_childtree(parent_tree,child_tree);

    print_parent_tree(parent_tree,"parent tree");
    print_child_tree(child_tree,"child tree");

    //head needs to store all inlabels, must be slightly bigger than n
    parlay::sequence<int> head(n+1);
    parlay::sequence<LCAnode<int>> augmented_vertices(n);
    for (int i = 0; i < n; i++) {
        augmented_vertices[i].id=i; //id is index
    }

    preprocess_par(parent_tree,child_tree,root,augmented_vertices,head);

    std::cout << "printing augment" << std::endl;
    std::cout << "id,\t inl,\t preo,\t#ri0,\tlvl,\tascd,\t size" << std::endl;
    for (int i =0 ; i < augmented_vertices.size(); i++) {
        augmented_vertices[i].print(); 
    }
    std::cout << "starting queries" << std::endl;


    std::uniform_int_distribution<int> dis(0,n-1);
    for (int  i  = 0; i < NUM_TRIALS; i++) {
        int u = dis(gen);
        int v = dis(gen);
        auto real_ans = vanilla_lca(parent_tree,u,v,root);
        std::cout << "calculated real ans" << std::endl;

        auto ans = query(head,parent_tree,augmented_vertices,u,v);

        std::cout << "calculated query " << std::endl;
       
        std::cout << u << " " << v << " " << "LCA: " << ans << std::endl;

         if (ans != real_ans) {
            std::cout << "break abort" << std::endl;
            std::cout << "Said: " << ans << "real: " << real_ans << std::endl;
            std::exit(5);
        }
    }

    }


}


//issue: two main functions -- move some things to headers? think about later. 
int main(int argc, char* argv[]) {
    std::cout << "hello world!2" << std::endl;

    int n=10; //defaults
    int NUM_TRIALS=100;
    int NUM_TREES = 1;
    int seed = 42; //fixed seed for easier testing

    parse_input(argc,argv,n,NUM_TRIALS,seed,NUM_TREES);

    std::mt19937 gen(seed);

    test_lca(n,NUM_TRIALS,NUM_TREES,gen);




}