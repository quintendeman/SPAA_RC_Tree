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
#include "../examples/samplesort.h"
#include<stdio.h>

//the functions here
#include "static_seqn_LCA.h"
#include "RC.h" //for generate random tree
#include "random_trees.h"

#include "VanillaLCA.h"



void test_lca(int n, int NUM_TREES, int NUM_TRIALS, std::mt19937& gen) {

    for (int iter = 0; iter < NUM_TREES; iter++) {

        //generate random parent tree
        //parlay::sequence<int> parent_tree = // generate_random_tree(n,gen);
        parlay::sequence<int> parent_tree = give_example_tree2();
        int root = 1; //set back to 0*
        n = 7; //for this example

        //from the parent tree, get a child_tree
        parlay::sequence<parlay::sequence<int>> child_tree(n,parlay::sequence<int>());
        partree_to_childtree(parent_tree,child_tree);

        print_parent_tree(parent_tree,"parent tree");
        print_child_tree(child_tree,"child tree");

        //head needs to store all inlabels, must be n+1 not n
        parlay::sequence<int> head(n+1);
        parlay::sequence<LCAnode<int>> augmented_vertices(n);
        for (int i = 0; i < n; i++) {
            augmented_vertices[i].id=i; //id is index
        }

        preprocess(parent_tree,child_tree,root,augmented_vertices,head);//

        std::cout << "printing augment" << std::endl;
        std::cout << "id,\t inl,\t ino,\t preo,\t#ri0,\tlvl,\tascd,\t size" << std::endl;
        for (int i =0 ; i < augmented_vertices.size(); i++) {
            augmented_vertices[i].print(); 
        }


        std::uniform_int_distribution<int> dis(0,n-1);
        for (int  i  = 0; i < NUM_TRIALS; i++) {
            int u = dis(gen);
            int v = dis(gen);
            auto ans = query(head,parent_tree,augmented_vertices,u,v);
            auto real_ans = vanilla_lca(parent_tree,u,v,root);
        
            std::cout << u << " " << v << " " << "LCA: " << ans << std::endl;

            if (ans != real_ans) {
                std::cout << "break abort" << std::endl;
                std::cout << "Said: " << ans << "real: " << real_ans << std::endl;
                std::exit(5);
            }
        }

    }


}

//test the folklore 'trick' that goes from fixed root LCA to arbitrary LCA
void test_unrooted_lca(int n, int NUM_TREES, int NUM_TRIALS, std::mt19937& gen) {

    for (int iter = 0; iter < NUM_TREES; iter++) {
    //generate random parent tree
    parlay::sequence<int> parent_tree = generate_random_tree(n,gen,.3);

    //from the parent tree, get a child_tree
    parlay::sequence<parlay::sequence<int>> child_tree(n,parlay::sequence<int>());
    partree_to_childtree(parent_tree,child_tree);
    auto unrooted_tree = make_undirected_tree(child_tree,parent_tree);

    print_parent_tree(parent_tree,"parent tree");
    print_child_tree(child_tree,"child tree");
    print_child_tree(unrooted_tree,"unrooted tree");


    std::uniform_int_distribution<int> dis(0,n-1);
    for (int  i  = 0; i < NUM_TRIALS; i++) {
        int u = dis(gen);
        int v = dis(gen);
        int rprime = dis(gen);
        auto ans = unrooted_lca(parent_tree,u,v,0,rprime);
        //std::cout << "finished fixed root" << std::endl;
        auto real_ans = true_unrooted_lca(unrooted_tree,u,v,rprime);
       
        std::cout << u << " " << v << " " << rprime << " " << "LCA: " << ans << std::endl;

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
    int NUM_TRIALS=1;
    int NUM_TREES=1;
    int seed = 42; //fixed seed for easier testing

    parse_input(argc,argv,n,NUM_TRIALS,seed,NUM_TREES);

    std::mt19937 gen(seed);

    test_lca(n,NUM_TREES,NUM_TRIALS,gen);




}