#ifndef LCA_TEST_H
#define LCA_TEST_H

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/internal/get_time.h"


//file to hold different lca tests
template<typename T, typename D>
void test_lca(int n, int NUM_TRIALS, int NUM_TREES, int k, std::mt19937& gen,parlay::random_generator& pgen, double forest_ratio, double chain_ratio, bool extra_print=false);

// //more extensive LCA test (lots of randomness)
void grand_test_lca(std::mt19937& gen,parlay::random_generator& pgen, int iters) {
    std::uniform_int_distribution<int> disn(1,100000000);
    std::uniform_int_distribution<int> dis_tri(100);
    std::uniform_int_distribution<int> dis_tree(10);
    std::uniform_int_distribution<int> dis_k(10000,50000); 
    std::uniform_real_distribution<double> dis_forest(0,1);
    std::uniform_real_distribution<double> dis_chain(0,1);

    for (int i = 0; i < iters; i++) {
        int n = disn(gen);
        int trials = 10; 
        int num_trees = 3; 
        int k = dis_k(gen);
        double forest_ratio = dis_forest(gen);
        double chain_ratio = dis_chain(gen);


        test_lca<int,int>(n,trials,num_trees,k,gen,pgen,forest_ratio,chain_ratio);
    }

}

// //choose a VERY large tree to test on
void large_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    test_lca<int,int>(300'000'000, 100, 2, 10000, gen, pgen, 0, 0.3);

}


//very small trees
//~6s of testing time
template<typename T, typename D>
void small_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    //test_lca(n, # trials, # trees, k, rand gen, rand ||gen, forest ratio, chain ratio);
   //std::cout << "1" << std::endl;
    test_lca<T,D>(1,2,100,1,gen,pgen,0,.3);
    //    std::cout << "1" << std::endl;
    test_lca<T,D>(2,10,20,2,gen,pgen,0.3,.3); 
      //  std::cout << "1" << std::endl;
        test_lca<T,D>(3,10,20,3,gen,pgen,0.15,.3); 
      //  std::cout << "1" << std::endl;

    test_lca<T,D>(4,10,100,3,gen,pgen,0.5,.3); 
     //  std::cout << "1" << std::endl;

    test_lca<T,D>(5,10,100,5,gen,pgen,0.25,.3); 
      //  std::cout << "1" << std::endl;

    test_lca<T,D>(6,10,10,10000,gen,pgen,0.2,.3); //note that k >> n is okay
      //  std::cout << "1" << std::endl;

    test_lca<T,D>(7,10,100,7,gen,pgen,0.1,.3); 
       // std::cout << "1" << std::endl;

    test_lca<T,D>(8,10,100,100,gen,pgen,0.6,.3);  
     //   std::cout << "1" << std::endl;

    test_lca<T,D>(9,10,200,20,gen,pgen,0.05,.3); 
    //   std::cout << "1" << std::endl;

   test_lca<T,D>(10,10,100,17,gen,pgen,0.01,.3); 
    //   std::cout << "1" << std::endl;
}

// //a chain (large and mid sized)
// ~ 300 mil tree budget
void chain_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    parlay::internal::timer t;
    t.start();
    std::cout << "about to start chain 5 " << t.next_time() << std::endl;
    test_lca<int,int>(5,100,3,10,gen,pgen,0,1,false);
        std::cout << "about to start chain 50 " << t.next_time() << std::endl;

    test_lca<int,int>(50,100,3,10,gen,pgen,0,1);
        std::cout << "about to start chain 500: " << t.next_time() << std::endl;

    test_lca<int,int>(500,100,3,50,gen,pgen,0,1);
        std::cout << "about to start chain 5000: " << t.next_time() << std::endl;

    test_lca<int,int>(5000,100,3,50,gen,pgen,0.001,1);

        std::cout << "about to start chain 5000'2: "  << t.next_time()  << std::endl;

    test_lca<int,int>(5000,100,3,50,gen,pgen,0,1);
        std::cout << "about to start chain 100K: "  << t.next_time() << std::endl;

    test_lca<int,int>(100000,100,3,1000,gen,pgen,0,1);
        std::cout << "about to start chain 1mil: "  << t.next_time()  << std::endl;

    test_lca<int,int>(1000000,10,3,1000,gen,pgen,0,1);
        std::cout << "about to start chain 10mil: "  << t.next_time() << std::endl;

    test_lca<int,int>(100000000,5,3,100,gen,pgen,0,1);
        std::cout << "finished chain 10mil: " << t.next_time()  << std::endl;


}

// // //no chains, if a node x has a child y that is not a leaf, then x has 2 children
// //large and small
// < 500 mil tree budget
void balanced_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    test_lca<int,int>(6,100,3,20,gen,pgen,0,0);
    test_lca<int,int>(51,100,3,15,gen,pgen,0,0);
    test_lca<int,int>(400,100,3,40,gen,pgen,0,0);
    test_lca<int,int>(3000,50,10,9000000,gen,pgen,0.01,0); //high k balanced few trees
    test_lca<int,int>(5000,100,3,1000000,gen,pgen,0,0); //high k balanced single tree
    test_lca<int,int>(5300,100,3,70,gen,pgen,0,0);
    test_lca<int,int>(5500,100,3,100,gen,pgen,0.001,0);
    test_lca<int,int>(90000,100,3,900,gen,pgen,0,0);
    test_lca<int,int>(1200000,10,3,3000,gen,pgen,0,0);
    test_lca<int,int>(1200000,10,3,3000,gen,pgen,0.5,0); //many forests
    test_lca<int,int>(90000000,11,3,3000,gen,pgen,0,0);

}

//k ~ n test
void n_k_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    test_lca<int,int>(10,10,10,10,gen,pgen,.1,.3);
    test_lca<int,int>(200,5,5,210,gen,pgen,.1,.3);
    test_lca<int,int>(10000,5,5,7000,gen,pgen,.01,.2);
    test_lca<int,int>(100000000,10,3,100005000,gen,pgen,.001,.3);
    test_lca<int,int>(70000000,10,3,50000000,gen,pgen,.05,.3);

}

// //large k testing
void large_k_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    test_lca<int,int>(10,10,10,2000,gen,pgen,.1,.3);
    test_lca<int,int>(11,5,5,20000,gen,pgen,.1,.3);
    test_lca<int,int>(1000,5,5,10001,gen,pgen,.01,.2);
    test_lca<int,int>(500,15,2,10001,gen,pgen,0,1); //all chains

}

// //large k testing
// ~ 1 bil tree budget
void small_k_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    test_lca<int,int>(10,10,10,1,gen,pgen,.1,.3);
    test_lca<int,int>(200,5,5,2,gen,pgen,.1,.3);
    test_lca<int,int>(10000,5,5,1,gen,pgen,.01,.2);
    test_lca<int,int>(100000000,1000,10,2,gen,pgen,.001,.3);
    test_lca<int,int>(70000000,1000,10,1,gen,pgen,.05,.3);

}

// //test where all vertices isolated
// //large and small
void isolated_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    test_lca<int,int>(6,100,3,20,gen,pgen,1,0.3);
    test_lca<int,int>(51,100,3,15,gen,pgen,1,0.3);
    test_lca<int,int>(400,100,3,40,gen,pgen,1,0.3);
    test_lca<int,int>(5500,100,3,100,gen,pgen,1,0.3);
}
// //many forests, not much connection
void mostly_forests_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    test_lca<int,int>(12,50,3,3,gen,pgen,.8,0.3);
    test_lca<int,int>(101,120,10,1000,gen,pgen,.8,0.2);
    test_lca<int,int>(306,80,20,50,gen,pgen,.8,0.1);
    test_lca<int,int>(8777,110,50,80,gen,pgen,.8,0.4);
    test_lca<int,int>(150000,90,100,9000,gen,pgen,.9,0.5);
    test_lca<int,int>(3200000,1000,1000,100,gen,pgen,.95,0.25);
    test_lca<int,int>(10000000,11,10,30000,gen,pgen,0.999,0.32);
}

void few_forests_lca(std::mt19937& gen, parlay::random_generator& pgen)  {

    test_lca<int,int>(6,100,3,20,gen,pgen,.0001,0.3);
    test_lca<int,int>(51,100,3,15,gen,pgen,.0001,0.2);
    test_lca<int,int>(400,100,3,40,gen,pgen,.0001,0.1);
    test_lca<int,int>(5500,100,3,100,gen,pgen,.0001,0.4);
    test_lca<int,int>(90000,100,3,900,gen,pgen,.0001,0.5);
    test_lca<int,int>(1200000,10,3,3000,gen,pgen,.000001,0.25);
    test_lca<int,int>(1200000,10,3,3000,gen,pgen,.00001,1); //all chains
    test_lca<int,int>(90000000,11,3,3000,gen,pgen,0.00000001,0.32);


}

// //one huge tree
//1.5 bil + 3 bil tree budget (fairly large)
void single_tree_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    test_lca<int,int>(150000000,100,10,50000,gen,pgen,0,.2);
    test_lca<int,int>(120020409,90,20,10001,gen,pgen,0,.4);

}
//benchmarking times TOD2*

//highlights of testing, for quick on laptop
template<typename T, typename D>
void mid_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    small_test_lca<T,D>(gen,pgen);
    test_lca<T,D>(5500,100,3,100,gen,pgen,1,0.3);
    test_lca<T,D>(10000001,100,3,50000,gen,pgen,.11,.2);
    test_lca<T,D>(200000,79,100,7999,gen,pgen,.82,0.6);
    test_lca<T,D>(10000,5,5,7000,gen,pgen,.01,.22);
    test_lca<T,D>(100000,100,3,1000,gen,pgen,0,1);

}

void short_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    small_test_lca<short,short>(gen,pgen);
    test_lca<short,short>(7000,100,3,555,gen,pgen,.2,.25);
    test_lca<short,short>(16000,100,10,203,gen,pgen,0,.2);
    test_lca<short,short>(200,10,10,40'000,gen,pgen,.01,.4);

}

void extensive_lca(std::mt19937& gen, parlay::random_generator& pgen,parlay::internal::timer& tim) {
    
    
    few_forests_lca(gen,pgen);
    std::cout << "Few forests done in " << tim.next_time() << std::endl;

    single_tree_lca(gen,pgen);
    std::cout << "Single tree done in " << tim.next_time() << std::endl;

    std::cout << "finished tests of interest, exit" << std::endl;
    exit(7);


    large_test_lca(gen,pgen);
    std::cout << "Large test done in " << tim.next_time() << std::endl;

    small_test_lca<int,int>(gen, pgen); //test small #s (edge case)
    std::cout << "Small test done in " << tim.next_time() << std::endl;

    mid_test_lca<long,long>(gen,pgen);
    std::cout << "mid test long done in " << tim.next_time() << std::endl;

    short_test_lca(gen,pgen);
    std::cout << "small test short type done in " << tim.next_time() << std::endl;

    grand_test_lca(gen,pgen,10); //just lots of random values
    std::cout << "Grand test done in " << tim.next_time() << std::endl;

    isolated_test_lca(gen,pgen); //isolated vs connected
        std::cout << "Isolated test done in " << tim.next_time() << std::endl;

    mostly_forests_lca(gen,pgen);
        std::cout << "Mostly forests test done in " << tim.next_time() << std::endl;

   

    chain_test_lca(gen,pgen); //vary balanced vs chain formation
        std::cout << "Chain test done in " << tim.next_time() << std::endl;

    balanced_test_lca(gen,pgen);
        std::cout << "Balanced test done in  " << tim.next_time() << std::endl;


    large_k_test_lca(gen,pgen);
        std::cout << "Large k test done in " << tim.next_time() << std::endl;

    small_k_test_lca(gen,pgen);

    std::cout << "Small k test done in " << tim.next_time() << std::endl;

    n_k_test_lca(gen,pgen);
    std::cout << "n~k test done in " << tim.next_time() << std::endl;


    //TOD2* test on T types other than int, Ex. custom T types (user defined structs?)

}

//designed to fail?
void fail_test(std::mt19937& gen, parlay::random_generator& pgen) {

    test_lca<int,int>(8,10,100,100,gen,pgen,0.6,.3,true);  

    //how would this line affect the next?? (and each parameter here matters*)
    test_lca<int,int>(9,10,3,20,gen,pgen,0,.3,true); 

}

#endif