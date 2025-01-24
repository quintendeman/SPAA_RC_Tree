//file to test LCA with
//TOD2* what happens when running on -03 optimization flag? Would code still pass tests? 

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

const bool PRINT_B = false;

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

void get_RC_tree(parlay::sequence<cluster<int,int>>& clusters,  parlay::sequence<int>& parent_tree) {
    parlay::sequence<parlay::sequence<int>> G;
    G=convert_parents_to_graph(G,parent_tree);
    create_base_clusters(G, clusters, 3); //max deg=3; each vertex in the tree has maximum degree 3
    create_RC_tree(clusters, static_cast<int>(parent_tree.size()), 0, [] (int A, int B) {return A+B;}, false,false); //randomized=false -- use the deterministic contraction method //defretval=0 -- default cluster val (ex -inf, inf, 0) //print=false -- don't print contraction sizes
}

//TOD2* test on T types other than ints? 
void handle_answers(parlay::sequence<std::tuple<int,int,int>>& queries, parlay::sequence<cluster<int,int>*>& answers, int k, parlay::sequence<int>& parent_tree, parlay::sequence<cluster<int,int>>& clusters, int iter, int iter2, double forest_ratio, double chain_ratio) {

    std::cout << "y" << std::endl;
    parlay::sequence<int> real_answers = parlay::tabulate(k,[&] (size_t i) {
        int u = std::get<0>(queries[i]);
        int v = std::get<1>(queries[i]);
        int rprime = std::get<2>(queries[i]);
        return unrooted_forest_lca(parent_tree,u,v,rprime);
    });
    std::cout << "check" << std::endl;
    //std::cout << "<";
    for (int i  = 0; i < k; i++) {
        int u = std::get<0>(queries[i]);
        int v = std::get<1>(queries[i]);
        int rprime = std::get<2>(queries[i]);
        //this version uses the fixed to arbitrary LCA trick for the real ans (but this was tested so okay)
    
        // std::cout << u << " " << v << " " << rprime << " LCA: " << real_answers[i];
        // if (i == k-1) std::cout << ">";    
        // std::cout << std::endl;

        if (answers[i] == nullptr && real_answers[i]==-1 ) {
            if (PRINT_B) std::cout << "lca queried nodes not connected, hence -1" << std::endl;
        }

        else if (answers[i]==nullptr || answers[i]->index != real_answers[i]) {
            std::cout << "break abort" << std::endl;
            if (answers[i]==nullptr) {
                std::cout << "answers[i] is nullptr, real: " << real_answers[i] << std::endl;
            }
            else {
                std::cout << "Said: " << answers[i]->index << "real: " << real_answers[i] << std::endl;

            }
            
            std::cout << "Tree# " << iter << ", " << "Trial: " << iter2 << std::endl;
            std::cout << "Stats: n: " << parent_tree.size() << " k: " << k << " forest-ratio: " << forest_ratio << ", chain-ratio: " << chain_ratio << std::endl;

            //see the clusters
            //printTree(clusters);

            //print_parent_tree(parent_tree,"parent tree");

            std::cout << "bad exit" << std::endl;

            std::exit(5);
        }
    }
}

void run_tree(int NUM_TRIALS, parlay::sequence<cluster<int,int>>& clusters, int k,  parlay::sequence<int>& parent_tree, std::mt19937& gen,  std::uniform_int_distribution<int>& dis, int iter, double forest_ratio, double chain_ratio) {

    for (int iter2 = 0; iter2 < NUM_TRIALS; iter2++) {
        parlay::sequence<std::tuple<int,int,int>> queries;

        //TOD2* test using variable batch sizes on same tree?
        for (int i = 0; i < k; i++) {
            int u = dis(gen);
            int v = dis(gen);
            int r = dis(gen);
            queries.push_back(std::make_tuple(u,v,r));
        }

        //NOTE! answers size initialized here*
        parlay::sequence<cluster<int,int>*> answers(k);

        if (PRINT_B) std::cout << "about to batch" << std::endl;

        batchLCA(clusters,queries,answers);

        if (PRINT_B) std::cout << "done batching" << std::endl;

        handle_answers(queries,answers,k,parent_tree,clusters,iter,iter2,forest_ratio,chain_ratio);


        
        //reset on each print debug
        if (PRINT_B) std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;

        }


}


//TOD2* adjust LCA code to work with forest (multiple roots!) (all tests so far on 1 tree only)
void test_lca(int n, int NUM_TRIALS, int NUM_TREES, int k, std::mt19937& gen,parlay::random_generator& pgen, double forest_ratio, double chain_ratio) {

    std::uniform_int_distribution<int> dis(0,n-1);

    for (int iter = 0; iter < NUM_TREES; iter++) {
        //generate random parent tree (permuted)
        parlay::sequence<int> parent_tree = generate_random_tree_perm_par(n,pgen,chain_ratio); //TOD2* switch to random perm tree generate_random_tree(n,gen); // 
        ////generate_random_tree_perm(n,gen,chain_ratio);//

        divide_into_forests_par(n,parent_tree,forest_ratio,pgen);
        //see the clusters
        //print_parent_tree(parent_tree,"parent tree:");
        //for debugging
        // parlay::sequence<parlay::sequence<int>> child_tree(parent_tree.size(),parlay::sequence<int>());
        // partree_to_childtree(parent_tree,child_tree);
        // print_child_tree(child_tree,"child tree");

        //create the RC tree
        parlay::sequence<cluster<int,int>> clusters;
        get_RC_tree(clusters,parent_tree);
      
        //sanity check that the RC tree is valid before main test
        //TOD2* reinstall this check (subtree sums check currently failing)
        // test_rc_valid(parent_tree, clusters,true); //true for print

   
        //printTree(clusters);
       
        std::cout << "num roots: " << num_roots(parent_tree) << std::endl;


        run_tree(NUM_TRIALS,clusters,k,parent_tree,gen,dis,iter,forest_ratio,chain_ratio);

        deleteRCtree(clusters);


    }
}

    
// // //more extensive LCA test (lots of randomness)
// void grand_test_lca(std::mt19937& gen,parlay::random_generator& pgen, int iters) {
//     std::uniform_int_distribution<int> disn(1,100000000);
//     std::exponential_distribution<int> dis_tri(1000);
//     std::exponential_distribution<int> dis_tree(100);
//     std::normal_distribution<int> dis_k(10000,50000); //make sure to take abs value then add 1
//     std::uniform_real_distribution<double> dis_forest(0,1);
//     std::uniform_real_distribution<double> dis_chain(0,1);

//     for (int i = 0; i < iters; i++) {
//         int n = dis(gen);
//         int trials = dis_tri(gen);
//         int num_trees = dis_tree(gen);
//         int k = abs(dis_k(gen))+1;
//         double forest_ratio = dis_forest(gen);
//         double chain_ratio = dis_chain(gen);


//         test_lca(n,trials,num_trees,k,gen,pgen,forest_ratio,chain_ratio);
//     }

// }

// // //choose a VERY large forest to test on
// void large_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
//     test_lca(500'000'000, 100, 2, 10000, gen, pgen, 0, 0.3);

// }


// //very small trees
// void small_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
//     //test_lca(n, # trials, # trees, k, rand gen, rand ||gen, forest ratio, chain ratio);
//     test_lca(1,2,100,1,gen,pgen,0,.3);
//     test_lca(2,10,20,2,gen,pgen,0.3,.3); 
//     test_lca(3,10,20,3,gen,pgen,0.15,.3); 
//     test_lca(4,10,100,3,gen,pgen,0.5,.3); 
//     test_lca(5,10,100,5,gen,pgen,0.25,.3); 
//     test_lca(6,10,10,1000,gen,pgen,0.2,.3); //note that k >> n is okay
//     test_lca(7,10,100,7,gen,pgen,0.1,.3); 
//     test_lca(8,10,100,100,gen,pgen,0.6,.3);  
//     test_lca(9,10,100,20,gen,pgen,0.05,.3); 
//     test_lca(10,10,100,17,gen,pgen,0.01,.3); 

// }

// //a chain (large and mid sized)
// void chain_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
//     test_lca(5,100,3,10,gen,pgen,0,1);
//     test_lca(50,100,3,10,gen,pgen,0,1);
//     test_lca(500,100,3,50,gen,pgen,0,1);
//     test_lca(5000,100,3,50,gen,pgen,0.001,1);
//     test_lca(5000,100,3,50,gen,pgen,0,1);
//     test_lca(100000,100,3,1000,gen,pgen,0,1);
//     test_lca(1000000,10,3,1000,gen,pgen,0,1);
//     test_lca(100000000,10,3,1000,gen,pgen,0,1);

// }

// // //no chains, if a node x has a child y that is not a leaf, then x has 2 children
// //large and small
// void balanced_test_lca(std::mt19937& gen) {
//     test_lca(6,100,3,20,gen,pgen,0,0);
//     test_lca(51,100,3,15,gen,pgen,0,0);
//     test_lca(400,100,3,40,gen,pgen,0,0);
//     test_lca(3000,50,10,9000000,gen,pgen,0.01,0); //high k balanced few trees
//     test_lca(5000,100,3,1000000,gen,pgen,0,0); //high k balanced single tree
//     test_lca(5300,100,3,70,gen,pgen,0,0);
//     test_lca(5500,100,3,100,gen,pgen,0.001,0);
//     test_lca(90000,100,3,900,gen,pgen,0,0);
//     test_lca(1200000,10,3,3000,gen,pgen,0,0);
//     test_lca(1200000,10,3,3000,gen,pgen,0.5,0); //many forests
//     test_lca(90000000,11,3,3000,gen,pgen,0,0);

// }

// //large k testing
// void large_k_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
//     test_lca(10,10,10,2000,gen,pgen,.1,.3);
//     test_lca(11,5,5,20000,gen,pgen,.1,.3);
//     test_lca(1000,5,5,10001,gen,pgen,.01,.2);
//     test_lca(500,15,2,10001,gen,pgen,0,1); //all chains

// }

// //large k testing
// void small_k_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
//     test_lca(10,10,10,1,gen,pgen,.1,.3);
//     test_lca(200,5,5,2,gen,pgen,.1,.3);
//     test_lca(10000,5,5,1,gen,pgen,.01,.2);
//     test_lca(100000000,1000,10,2,gen,pgen,.001,.3);
//     test_lca(70000000,1000,10,1,gen,pgen,.05,.3);

// }

// // //test where all vertices isolated
// // //large and small
// void isolated_test_lca(std::mt19937& gen, parlay::random_generator& pgen)) {
//     test_lca(6,100,3,20,gen,pgen,1,0.3);
//     test_lca(51,100,3,15,gen,pgen,1,0.3);
//     test_lca(400,100,3,40,gen,pgen,1,0.3);
//     test_lca(5500,100,3,100,gen,pgen,1,0.3);
// }
// // //many forests, not much connection
// void mostly_forests_lca(std::mt19937& gen, parlay::random_generator& pgen) {
//     test_lca(12,50,3,3,gen,pgen,.8,0.3);
//     test_lca(101,120,10,1000,gen,pgen,.8,0.2);
//     test_lca(306,80,20,50,gen,pgen,.8,0.1);
//     test_lca(8777,110,50,80,gen,pgen,.8,0.4);
//     test_lca(150000,90,100,9000,gen,pgen,.9,0.5);
//     test_lca(3200000,1000,1000,100,gen,pgen,.95,0.25);
//     test_lca(10000000,11,10,30000,gen,pgen,0.999,0.32);
// }

// void few_forests_lca(std::mt19937& gen, parlay::random_generator& pgen)  {

//     test_lca(6,100,3,20,gen,pgen,.8,0.3);
//     test_lca(51,100,3,15,gen,pgen,.8,0.2);
//     test_lca(400,100,3,40,gen,pgen,.8,0.1);
//     test_lca(5500,100,3,100,gen,pgen,.8,0.4);
//     test_lca(90000,100,3,900,gen,pgen,.9,0.5);
//     test_lca(1200000,10,3,3000,gen,pgen,.95,0.25);
//     test_lca(1200000,10,3,3000,gen,pgen,.95,1); //all chains
//     test_lca(90000000,11,3,3000,gen,pgen,0.999,0.32);


// }

// //one huge tree
// void single_tree_lca(std::mt19937& gen, parlay::random_generator& pgen) {
//     test_lca(150000000,100,10,50000,gen,pgen,1,.2);
//     test_lca(120020409,90,20,10001,gen,pgen,1,.4);

// }


// void extensive_lca(std::mt19937& gen, parlay::random_generator& pgen) {
//     small_test_lca(gen, pgen); //test small #s (edge case)

//     chain_test_lca(gen,pgen); //vary balanced vs chain formation
//     balanced_test_lca(gen,pgen);

//     large_k_test_lca(gen,pgen);
//     small_k_test_lca(gen,pgen);

//     isolated_test_lca(gen,pgen); //isolated vs connected
//     mostly_forests_lca(gen,pgen);
//     few_forests_lca(gen,pgen);
//     single_tree_lca(gen,pgen);

//     grand_test_lca(gen,pgen,100); //just lots of random values

//     large_test_lca(gen,pgen);

//     //TOD2* test on T types other than int, Ex. custom T types (user defined structs?)

// }

void extensive_test_perm(int n, int NUM_TREES, std::mt19937& gen) {
    for (int i = 0; i < NUM_TREES; i++) {
        auto parents = generate_random_tree(n,gen,.3); //chain ratio is .3
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
    auto parents = generate_random_tree_perm(n,gen,.3);
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
    int k = 1;
    int seed = 42; //fixed seed for easier testing
    int pseed = 43;
    double forest_ratio = 5.0/n;//1.41 / n; //ratio scales with n to not disconnect too much

    parse_input(argc,argv,n,NUM_TRIALS,seed,NUM_TREES,k);

    std::mt19937 gen(seed);
    parlay::random_generator pgen(pseed);

    std::cout << "Batch size: " << k << std::endl;

    test_lca(n,NUM_TRIALS,NUM_TREES,k,gen,pgen,0,.3); //0 is forest ratio
    //extensive_test_perm(n,NUM_TREES,gen);

}