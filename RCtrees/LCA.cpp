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

#include "RCdynamic.h"
#include "RC_test.h"

//need our LCA funs
#include "LCA.h"
#include "VanillaLCA.h"
#include "random_trees.h"


//TOD2* is there a better way to do this?
//Find the root of the RC tree, store in ans.
//U -- cluster in this RC tree (needed so we know what tree we are looking at).
//Takes O(log n) time (time in height of RC tree)
template<typename T, typename D>
void get_root_RC(parlay::sequence<cluster<T,D>>& clusters, cluster<T,D>*& ans) {
    
    ans=&clusters[0]; //just get a random cluster to start with //TODO change from 0? 

    while (ans-> parent != nullptr)  {
        ans = ans->parent;
      
    }
 }
//TOD2* adjust LCA code to work with forest (multiple roots!) (all tests so far on 1 tree only)
void test_lca(int n, int NUM_TRIALS, int NUM_TREES, int BATCH_SIZE, std::mt19937& gen) {

    std::uniform_int_distribution<int> dis(0,n-1);


    for (int iter = 0; iter < NUM_TREES; iter++) {
    //generate random parent tree (permuted)
    parlay::sequence<int> parent_tree = generate_random_tree_perm(n,gen); //TOD2* switch to random perm tree generate_random_tree(n,gen); //
    int init_root = get_root(parent_tree); //default for generate_random_tree is 0; since permuting, must find explicitly

    //from the parent tree, get a child_tree
    parlay::sequence<parlay::sequence<int>> child_tree(n,parlay::sequence<int>());
    partree_to_childtree(parent_tree,child_tree);

    print_parent_tree(parent_tree,"parent tree");
    print_child_tree(child_tree,"child tree");

    //create the RC tree
    parlay::sequence<parlay::sequence<int>> G;
    G=convert_parents_to_graph(G,parent_tree);
    parlay::sequence<cluster<int,int>> clusters;
    create_base_clusters(G, clusters, 3); //max deg=3; each vertex in the tree has maximum degree 3
    create_RC_tree(clusters, n, false,0,false); //randomized=false -- use the deterministic contraction method //defretval=0 -- default cluster val (ex -inf, inf, 0) //print=false -- don't print contraction sizes

    //sanity check that the RC tree is valid before main test
    test_rc_valid(parent_tree, clusters,false);

    //see the clusters
    printTree(clusters);

    cluster<int,int>* root;
    get_root_RC(clusters,root); //TOD2* I'm not counting this get root time (log n) against the alg time because I assume the RC tree knows its root (somewhere) for O(1) charge, I just don't know where to access it 


    std::cout << "root is " << root->index << std::endl;

    for (int iter2 = 0; iter2 < NUM_TRIALS; iter2++) {


        parlay::sequence<std::tuple<int,int,int>> queries;

        //TOD2* test using variable batch sizes on same tree?
        for (int i = 0; i < BATCH_SIZE; i++) {
            int u = dis(gen);
            int v = dis(gen);
            int r = dis(gen);
            queries.push_back(std::make_tuple(u,v,r));
        }

        parlay::sequence<cluster<int,int>*> answers;

        //std::cout << "about to batch" << std::endl;

        batchLCA(clusters,root,queries,answers);

        //std::cout << "done batching" << std::endl;


        for (int i  = 0; i < BATCH_SIZE; i++) {
            int u = std::get<0>(queries[i]);
            int v = std::get<1>(queries[i]);
            int rprime = std::get<2>(queries[i]);
            //this version uses the fixed to arbitrary LCA trick for the real ans (but this was tested so okay)
            auto real_ans = unrooted_lca(parent_tree,u,v,init_root,rprime);
        
            std::cout << u << " " << v << " " << "r'" << rprime << " LCA: " << real_ans << std::endl;

            if (answers[i]->index != real_ans) {
                std::cout << "break abort" << std::endl;
                std::cout << "Said: " << answers[i]->index << "real: " << real_ans << std::endl;
                std::cout << "Tree# " << iter << ", " << "Trial: " << iter2 << std::endl;

                //see the clusters
                printTree(clusters);

                print_parent_tree(parent_tree,"parent tree");
                print_child_tree(child_tree,"child tree");

                std::exit(5);
            }
        }
        //reset on each print debug
        std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;

        }

        deleteRCtree(clusters);


    }


}

void extensive_test_perm(int n, int NUM_TREES, std::mt19937& gen) {
    for (int i = 0; i < NUM_TREES; i++) {
        auto parents = generate_random_tree(n,gen);
        //auto perm_tree = generate_random_tree(n,gen); //fail example (example where trees not equal)
        auto perm_tree = permute_tree(n,parents,gen);
        auto root1 = get_root(parents);
        auto root2 = get_root(perm_tree);
        parlay::sequence<parlay::sequence<int>> child_tree1(n,parlay::sequence<int>());
        partree_to_childtree(parents,child_tree1);
        parlay::sequence<parlay::sequence<int>> child_tree2(n,parlay::sequence<int>());
        partree_to_childtree(perm_tree,child_tree2);
        if (!possibly_equal(child_tree1,child_tree2,root1,root2)) {
            std::cout << "Error abort, tree and permuted tree NOT equal" << std::endl;
            exit(3002);
        } 
        std::cout << "." << std::endl;


    }
    std::cout << "passed extensive perm test" << std::endl;

}

void test_perm() {
    int n = 10;
    int seed = 40;
    std::mt19937 gen(seed);
    auto parents = generate_random_tree_perm(n,gen);
    print_parent_tree(parents,"par tree");
    parlay::sequence<parlay::sequence<int>> child_tree(n,parlay::sequence<int>());
    partree_to_childtree(parents,child_tree);
    print_child_tree(child_tree,"child tree");
}

int main(int argc, char* argv[]) {
    std::cout << "hello world!3" << std::endl;

    int n = 10; //defaults
    int NUM_TRIALS = 1;
    int NUM_TREES = 1;
    int BATCH_SIZE = 1;
    int seed = 42; //fixed seed for easier testing

    parse_input(argc,argv,n,NUM_TRIALS,seed,NUM_TREES,BATCH_SIZE);

    std::mt19937 gen(seed);

    test_lca(n,NUM_TRIALS,NUM_TREES,BATCH_SIZE,gen);
    //extensive_test_perm(n,NUM_TREES,gen);

}