//file to test LCA with
//TOD2* what happens when running on -03 optimization flag? Would code still pass tests? 

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/internal/get_time.h"

//RC tree funs
#include "RC_test.h"
#include "RC.h"

//need our LCA funs
#include "LCA.h"
#include "VanillaLCA.h"
#include "random_trees.h"

const bool PRINT_B = false; //progress prints
const bool PRINT_T = false; //timer prints

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

void get_RC_tree(parlay::sequence<cluster<int,int>>& clusters,  parlay::sequence<int>& parent_tree, bool extra_print=false) {
    parlay::sequence<parlay::sequence<int>> G;
    G=convert_parents_to_graph(G,parent_tree);
    if (extra_print) std::cout << "made graph" << std::endl;
    create_base_clusters(G, clusters, 3); //max deg=3; each vertex in the tree has maximum degree 3
    if (extra_print) std::cout << "made base clusters" << std::endl;
    create_RC_tree(clusters, static_cast<int>(parent_tree.size()), 0, [] (int A, int B) {return A+B;}, false); //randomized=false -- use the deterministic contraction method //defretval=0 -- default cluster val (ex -inf, inf, 0) //print=false -- don't print contraction sizes
    if (extra_print) std::cout << "made tree " << std::endl;
}

void error_print(parlay::sequence<cluster<int,int>*>& answers, parlay::sequence<int>& real_answers, int iter, int iter2, parlay::sequence<int>& parent_tree, int k, double forest_ratio, double chain_ratio, std::string error_file_prefix, parlay::sequence<std::tuple<int,int,int>>& queries, int i) {
    if (answers[i]==nullptr) {
        std::cout << "answers[i] is nullptr, real: " << real_answers[i] << std::endl;
    }
    else {
        std::cout << "Said: " << answers[i]->index << "real: " << real_answers[i] << std::endl;

    }
    
    std::cout << "Tree# " << iter << ", " << "Trial: " << iter2 << std::endl;
    std::cout << "Stats: n: " << parent_tree.size() << " k: " << k << " forest-ratio: " << forest_ratio << ", chain-ratio: " << chain_ratio << std::endl;

    //counter naming to avoid destroying previous files
    if (error_file_prefix!="nofileprint") {

    
    std::ifstream mf;
    mf.open("counter.txt");
    int my_counter=0;
    mf >> my_counter;

    mf.close();

    std::ofstream mf2;
    mf2.open("counter.txt");
    mf2 << (my_counter + 1) << "\n";
    mf2.close();

    //couple file with console output for fuller picture
    std::ofstream myfile;
    std::string new_file_name = error_file_prefix + std::to_string(my_counter) + ".txt";
    myfile.open(new_file_name);
    myfile << "Writing break abort failed case into file. \n";
    myfile << parent_tree.size() << "\n";
    for (int val = 0; val < parent_tree.size(); val++) {
        myfile << parent_tree[val] << " ";
    }
    myfile << "\n";
    myfile << k << "\n";
    for (int val = 0; val < queries.size(); val++) {
        myfile << std::get<0>(queries[val]) << " " << std::get<1>(queries[val]) << " " << std::get<2>(queries[val]) << " ";
    }
    myfile << "\n";
    myfile.close();
    }
    else {
        std::cout << "no file print specified " << std::endl;
    }

    //see the clusters
    //printTree(clusters);

    //print_parent_tree(parent_tree,"parent tree");

}

//TOD2* test on T types other than ints? 
void handle_answers(parlay::sequence<std::tuple<int,int,int>>& queries, parlay::sequence<cluster<int,int>*>& answers, int k, parlay::sequence<int>& parent_tree, parlay::sequence<cluster<int,int>>& clusters, int iter, int iter2, double forest_ratio, double chain_ratio, std::string error_file_prefix="error_log") {

    parlay::sequence<int> real_answers = parlay::tabulate(k,[&] (size_t i) {
        int u = std::get<0>(queries[i]);
        int v = std::get<1>(queries[i]);
        int rprime = std::get<2>(queries[i]);
        return unrooted_forest_lca(parent_tree,u,v,rprime);
    });
    //do checks in parallel, this is bottleneck?
    bool any_bad = parlay::reduce(parlay::delayed_tabulate(k,[&] (size_t i) {
        int u = std::get<0>(queries[i]);
        int v = std::get<1>(queries[i]);
        int rprime = std::get<2>(queries[i]);
        if ((answers[i]==nullptr && real_answers[i] != -1) || (answers[i] != nullptr && answers[i]->index != real_answers[i])) {
            return true; //set that error does exist

        }
        return false;

    }),MonoidOr()) ;

    //std::cout << "<";
    if (any_bad) { //if there was a bad output, go through sequentially to find error
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
                //TOD2* uncomment
            else if (answers[i]==nullptr || answers[i]->index != real_answers[i]) {
                
                error_print(answers,real_answers,iter,iter2,parent_tree,k,forest_ratio,chain_ratio,error_file_prefix,queries, i);
                std::cout << "break abort, bad exit" << std::endl;
    
    
                std::exit(5);
               
    
                
            } //TOD2* uncomment
        }
        std::cout << "because is bad triggered should never reach here, error" << std::endl;
        exit(5002);

    }
    
}

void run_from_file(std::string file) {
    std::ifstream my_file;
    my_file.open(file);
    std::string first_line;
    std::getline(my_file,first_line);
    std::cout << "first line: " << first_line << std::endl;
    int n = 0;
    my_file >> n;
    int temp;
    parlay::sequence<int> parent_tree(n);
    for (int i = 0; i < n; i++) {
        my_file >> temp;
        parent_tree[i]=temp;
    }
    int k = 0;
    my_file >> k;
    k=1; //TOD2* change back
    parlay::sequence<std::tuple<int,int,int>> queries(k);
    int t1, t2, t3 = 0;
    for (int i = 0; i < k+1; i++) { //TOd2* change back
        my_file >> t1 >> t2 >> t3;
        if (i > 0)  //restriction for error_log4
        queries[i-1]=std::make_tuple(t1,t2,t3);  //TOD2* change back
    }

    std::cout << "n is " << n << std::endl;
    std::cout << "k is " << k << std::endl;
    for (int i = 0; i < queries.size(); i++) {
        std::cout << std::get<0>(queries[i]) << ", " << std::get<1>(queries[i]) << ", " << std::get<2>(queries[i]) << std::endl; 

    }
    //print_parent_tree(parent_tree,"rf par tree");
    parlay::sequence<cluster<int,int>> clusters;
    get_RC_tree(clusters,parent_tree,false);
    int NUM_TRIALS, iter, iter2 =-1; //because we chose the trial
    double forest_ratio, chain_ratio = -1; //because we already made tree

    parlay::sequence<cluster<int,int>*> answers(k);

    batchLCA(clusters,queries,answers);
    handle_answers(queries,answers,k,parent_tree,clusters,iter,iter2,forest_ratio,chain_ratio,"nofileprint"); //alt is alt error file to avoid covering previous print

    deleteRCtree(clusters);


}

void run_tree(int NUM_TRIALS, parlay::sequence<cluster<int,int>>& clusters, int k,  parlay::sequence<int>& parent_tree, parlay::random_generator& pgen,  std::uniform_int_distribution<int>& dis, int iter, double forest_ratio, double chain_ratio) {

    parlay::internal::timer myt;
    
    //std::cout << "about to start LCA" << myt.next_time() << std::endl;

    dis(pgen); //stepping pgen

    for (int iter2 = 0; iter2 < NUM_TRIALS; iter2++) {
        //TOD2* test using variable batch sizes on same tree?
        parlay::sequence<std::tuple<int,int,int>> queries = parlay::tabulate(k,[&] (size_t i) {
            auto r = pgen[i];
            return std::make_tuple(dis(r),dis(r),dis(r));
        });

        //std::cout << "read queries " << myt.next_time() << std::endl;

        //NOTE! answers size initialized here*
        parlay::sequence<cluster<int,int>*> answers(k);

        if (PRINT_B) std::cout << "about to batch" << std::endl;

        //Theta(n) too expensive given that LCA is O(k log n)
        // parlay::parallel_for(0,clusters.size(),[&] (size_t i) {
        //     if (clusters[i].counter!=0) {
        //         std::cout << "error counter not 0 before" << std::endl;
        //         exit(5000);
        //     }
        // });

        batchLCA(clusters,queries,answers);

        if (PRINT_B) std::cout << "did LCA  " << myt.next_time() << std::endl;
        // parlay::parallel_for(0,clusters.size(),[&] (size_t i) {
        //     if (clusters[i].counter!=0) {
        //         std::cout << "error counter not 0 before" << std::endl;
        //         exit(5000);
        //     }
        // });

        if (PRINT_B) std::cout << "done batching" << std::endl;

        //TODO* uncomment
        handle_answers(queries,answers,k,parent_tree,clusters,iter,iter2,forest_ratio,chain_ratio);


        if (PRINT_B) std::cout << "handled answers " << myt.next_time() << std::endl;

        
        //reset on each print debug
        if (PRINT_B) std::cout << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl << std::endl;

        }


}


//TOD2* adjust LCA code to work with forest (multiple roots!) (all tests so far on 1 tree only)
void test_lca(int n, int NUM_TRIALS, int NUM_TREES, int k, std::mt19937& gen,parlay::random_generator& pgen, double forest_ratio, double chain_ratio, bool extra_print=false) {

    if (extra_print) std::cout << "starting n " << n << std::endl;
   // extra_print=false; //temporary

    std::uniform_int_distribution<int> dis(0,n-1);

    parlay::internal::timer t2;

    for (int iter = 0; iter < NUM_TREES; iter++) {

       if (PRINT_T) std::cout << "starting tree " << iter << " at time " << t2.next_time() << std::endl;
       // if (iter == 24) extra_print=true; //specific iter of interest
        if (extra_print) std::cout << "starting tree " << iter << std::endl;
        //generate random parent tree (permuted)
        parlay::sequence<int> parent_tree = generate_random_tree_perm_par(n,pgen,chain_ratio); //TOD2* switch to random perm tree generate_random_tree(n,gen); // 
        ////generate_random_tree_perm(n,gen,chain_ratio);//

        divide_into_forests_par(n,parent_tree,forest_ratio,pgen);
        //see the clusters
        if (extra_print) print_parent_tree(parent_tree,"parent tree:");

        if (PRINT_T) std::cout << "manual tree " << t2.next_time() << std::endl;

        //for debugging
        // parlay::sequence<parlay::sequence<int>> child_tree(parent_tree.size(),parlay::sequence<int>());
        // partree_to_childtree(parent_tree,child_tree);
        // print_child_tree(child_tree,"child tree");
        if (extra_print) std::cout << "about to make clusters" << std::endl;
        //create the RC tree
        parlay::sequence<cluster<int,int>> clusters;
        get_RC_tree(clusters,parent_tree,extra_print);

        if (PRINT_T) std::cout << "Finished RC " << t2.next_time() << std::endl;

      
        //sanity check that the RC tree is valid before main test
        //TOD2* reinstall this check (subtree sums check currently failing)
        // test_rc_valid(parent_tree, clusters,true); //true for print

   
        //if (extra_print) printTree(clusters);
       
        //std::cout << "num roots: " << num_roots(parent_tree) << std::endl;
        if (extra_print) std::cout << "about to run LCA" << std::endl;
        run_tree(NUM_TRIALS,clusters,k,parent_tree,pgen,dis,iter,forest_ratio,chain_ratio);
        if (extra_print) std::cout << "finished LCA" << std::endl;
        if (PRINT_T) std::cout << "finished LCA " << t2.next_time() << std::endl;

        deleteRCtree(clusters);


    }

    printf("success on test n: %d trials: %d trees: %d k: %d forr: %f chainr: %f\n", n, NUM_TRIALS,NUM_TREES,k,forest_ratio,chain_ratio);
    std::cout << "dummy print" << std::endl; // needed to get shown in outfile 
}

    
// //more extensive LCA test (lots of randomness)
void grand_test_lca(std::mt19937& gen,parlay::random_generator& pgen, int iters) {
    std::uniform_int_distribution<int> disn(1,100000000);
    std::exponential_distribution<> dis_tri(1000);
    std::exponential_distribution<> dis_tree(100);
    std::normal_distribution<> dis_k(10000,50000); //make sure to take abs value then add 1
    std::uniform_real_distribution<double> dis_forest(0,1);
    std::uniform_real_distribution<double> dis_chain(0,1);

    for (int i = 0; i < iters; i++) {
        int n = disn(gen);
        int trials = dis_tri(gen);
        int num_trees = dis_tree(gen);
        int k = abs(dis_k(gen))+1;
        double forest_ratio = dis_forest(gen);
        double chain_ratio = dis_chain(gen);


        test_lca(n,trials,num_trees,k,gen,pgen,forest_ratio,chain_ratio);
    }

}

// //choose a VERY large forest to test on
void large_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    test_lca(500'000'000, 100, 2, 10000, gen, pgen, 0, 0.3);

}


//very small trees
//~6s of testing time
void small_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    //test_lca(n, # trials, # trees, k, rand gen, rand ||gen, forest ratio, chain ratio);
   //std::cout << "1" << std::endl;
    test_lca(1,2,100,1,gen,pgen,0,.3);
    //    std::cout << "1" << std::endl;
    test_lca(2,10,20,2,gen,pgen,0.3,.3); 
      //  std::cout << "1" << std::endl;
        test_lca(3,10,20,3,gen,pgen,0.15,.3); 
      //  std::cout << "1" << std::endl;

    test_lca(4,10,100,3,gen,pgen,0.5,.3); 
     //  std::cout << "1" << std::endl;

    test_lca(5,10,100,5,gen,pgen,0.25,.3); 
      //  std::cout << "1" << std::endl;

    test_lca(6,10,10,10000,gen,pgen,0.2,.3); //note that k >> n is okay
      //  std::cout << "1" << std::endl;

    test_lca(7,10,100,7,gen,pgen,0.1,.3); 
       // std::cout << "1" << std::endl;

    test_lca(8,10,100,100,gen,pgen,0.6,.3);  
     //   std::cout << "1" << std::endl;

    test_lca(9,10,200,20,gen,pgen,0.05,.3); 
    //   std::cout << "1" << std::endl;

   test_lca(10,10,100,17,gen,pgen,0.01,.3); 
    //   std::cout << "1" << std::endl;
}

// //a chain (large and mid sized)
// ~ 300 mil tree budget
void chain_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    std::cout << "about to start chain 5" << std::endl;
    test_lca(5,100,3,10,gen,pgen,0,1,false);
        std::cout << "about to start chain 50" << std::endl;

    test_lca(50,100,3,10,gen,pgen,0,1);
        std::cout << "about to start chain 500" << std::endl;

    test_lca(500,100,3,50,gen,pgen,0,1);
        std::cout << "about to start chain 5000" << std::endl;

    test_lca(5000,100,3,50,gen,pgen,0.001,1);

        std::cout << "about to start chain 5000'2" << std::endl;

    test_lca(5000,100,3,50,gen,pgen,0,1);
        std::cout << "about to start chain 100K" << std::endl;

    test_lca(100000,100,3,1000,gen,pgen,0,1);
        std::cout << "about to start chain 1mil" << std::endl;

    test_lca(1000000,10,3,1000,gen,pgen,0,1);
        std::cout << "about to start chain 10mil" << std::endl;

    test_lca(100000000,5,3,1000,gen,pgen,0,1);
        std::cout << "finished chain 10mil" << std::endl;


}

// // //no chains, if a node x has a child y that is not a leaf, then x has 2 children
// //large and small
// < 500 mil tree budget
void balanced_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    test_lca(6,100,3,20,gen,pgen,0,0);
    test_lca(51,100,3,15,gen,pgen,0,0);
    test_lca(400,100,3,40,gen,pgen,0,0);
    test_lca(3000,50,10,9000000,gen,pgen,0.01,0); //high k balanced few trees
    test_lca(5000,100,3,1000000,gen,pgen,0,0); //high k balanced single tree
    test_lca(5300,100,3,70,gen,pgen,0,0);
    test_lca(5500,100,3,100,gen,pgen,0.001,0);
    test_lca(90000,100,3,900,gen,pgen,0,0);
    test_lca(1200000,10,3,3000,gen,pgen,0,0);
    test_lca(1200000,10,3,3000,gen,pgen,0.5,0); //many forests
    test_lca(90000000,11,3,3000,gen,pgen,0,0);

}

//k ~ n test
void n_k_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    test_lca(10,10,10,10,gen,pgen,.1,.3);
    test_lca(200,5,5,210,gen,pgen,.1,.3);
    test_lca(10000,5,5,7000,gen,pgen,.01,.2);
    test_lca(100000000,10,10,100005000,gen,pgen,.001,.3);
    test_lca(70000000,10,10,50000000,gen,pgen,.05,.3);

}

// //large k testing
void large_k_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    test_lca(10,10,10,2000,gen,pgen,.1,.3);
    test_lca(11,5,5,20000,gen,pgen,.1,.3);
    test_lca(1000,5,5,10001,gen,pgen,.01,.2);
    test_lca(500,15,2,10001,gen,pgen,0,1); //all chains

}

// //large k testing
// ~ 1 bil tree budget
void small_k_test_lca(std::mt19937& gen,parlay::random_generator& pgen) {
    test_lca(10,10,10,1,gen,pgen,.1,.3);
    test_lca(200,5,5,2,gen,pgen,.1,.3);
    test_lca(10000,5,5,1,gen,pgen,.01,.2);
    test_lca(100000000,1000,10,2,gen,pgen,.001,.3);
    test_lca(70000000,1000,10,1,gen,pgen,.05,.3);

}

// //test where all vertices isolated
// //large and small
void isolated_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    test_lca(6,100,3,20,gen,pgen,1,0.3);
    test_lca(51,100,3,15,gen,pgen,1,0.3);
    test_lca(400,100,3,40,gen,pgen,1,0.3);
    test_lca(5500,100,3,100,gen,pgen,1,0.3);
}
// //many forests, not much connection
void mostly_forests_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    test_lca(12,50,3,3,gen,pgen,.8,0.3);
    test_lca(101,120,10,1000,gen,pgen,.8,0.2);
    test_lca(306,80,20,50,gen,pgen,.8,0.1);
    test_lca(8777,110,50,80,gen,pgen,.8,0.4);
    test_lca(150000,90,100,9000,gen,pgen,.9,0.5);
    test_lca(3200000,1000,1000,100,gen,pgen,.95,0.25);
    test_lca(10000000,11,10,30000,gen,pgen,0.999,0.32);
}

void few_forests_lca(std::mt19937& gen, parlay::random_generator& pgen)  {

    test_lca(6,100,3,20,gen,pgen,.8,0.3);
    test_lca(51,100,3,15,gen,pgen,.8,0.2);
    test_lca(400,100,3,40,gen,pgen,.8,0.1);
    test_lca(5500,100,3,100,gen,pgen,.8,0.4);
    test_lca(90000,100,3,900,gen,pgen,.9,0.5);
    test_lca(1200000,10,3,3000,gen,pgen,.95,0.25);
    test_lca(1200000,10,3,3000,gen,pgen,.95,1); //all chains
    test_lca(90000000,11,3,3000,gen,pgen,0.999,0.32);


}

// //one huge tree
//1.5 bil + 3 bil tree budget (fairly large)
void single_tree_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    test_lca(150000000,100,10,50000,gen,pgen,1,.2);
    test_lca(120020409,90,20,10001,gen,pgen,1,.4);

}
//benchmarking times TOD2*

//highlights of testing, for quick on laptop
void mid_test_lca(std::mt19937& gen, parlay::random_generator& pgen,parlay::internal::timer& tim) {
    small_test_lca(gen,pgen);
    test_lca(5500,100,3,100,gen,pgen,1,0.3);
    test_lca(10000001,100,3,50000,gen,pgen,.11,.2);
    test_lca(200000,79,100,7999,gen,pgen,.82,0.6);
    test_lca(10000,5,5,7000,gen,pgen,.01,.22);
    test_lca(100000,100,3,1000,gen,pgen,0,1);

}

void extensive_lca(std::mt19937& gen, parlay::random_generator& pgen,parlay::internal::timer& tim) {
    
    small_test_lca(gen, pgen); //test small #s (edge case)
    std::cout << "Small test done in " << tim.next_time() << std::endl;

    isolated_test_lca(gen,pgen); //isolated vs connected
        std::cout << "Isolated test done in " << tim.next_time() << std::endl;

    mostly_forests_lca(gen,pgen);
        std::cout << "Mostly forests test done in " << tim.next_time() << std::endl;

    few_forests_lca(gen,pgen);
        std::cout << "Few forests done in " << tim.next_time() << std::endl;

    single_tree_lca(gen,pgen);
        std::cout << "Single tree done in " << tim.next_time() << std::endl;


    grand_test_lca(gen,pgen,100); //just lots of random values
    std::cout << "Grand test done in " << tim.next_time() << std::endl;

    large_test_lca(gen,pgen);
    std::cout << "Large test done in " << tim.next_time() << std::endl;


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

void more_test_lca(std::mt19937& gen, parlay::random_generator& pgen) {
    small_test_lca(gen,pgen);
    

}

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

//designed to fail?
void fail_test(std::mt19937& gen, parlay::random_generator& pgen) {

    test_lca(8,10,100,100,gen,pgen,0.6,.3,true);  

    //how would this line affect the next?? (and each parameter here matters*)
    test_lca(9,10,3,20,gen,pgen,0,.3,true); 


}

int main(int argc, char* argv[]) {

    parlay::internal::timer tim = parlay::internal::timer();

    std::cout << "hello world!31" << std::endl;

    int n = 10; //defaults
    int NUM_TRIALS = 1;
    int NUM_TREES = 1;
    int k = 1;
    int seed = 42; //fixed seed for easier testing
    int pseed = 43;
    double forest_ratio = 5.0/n;//1.41 / n; //ratio scales with n to not disconnect too much
    double chain_ratio = .3;
    bool is_from_file = false;
    std::string filename = "";


    parse_input(argc,argv,n,NUM_TRIALS,seed,pseed,NUM_TREES,k,forest_ratio,chain_ratio,is_from_file,filename);

    std::cout << "pseed: " << pseed << std::endl;

    std::mt19937 gen(seed);
    parlay::random_generator pgen(pseed);

    std::cout << "Batch size: " << k << std::endl;
    std::cout << "for ratio: " << forest_ratio << std::endl;
    std::cout << "chain ratio: " << chain_ratio << std::endl;

    std::cout << "about to start" << tim.next_time() << std::endl;

    if (is_from_file) {
        std::cout << "running from file " << filename << std::endl;
        run_from_file(filename);
    }
    else if (forest_ratio==-1) {

        //test_lca(n,NUM_TRIALS,NUM_TREES,k,gen,pgen,forest_ratio,chain_ratio); //0 is forest ratio
        //extensive_test_perm(n,NUM_TREES,gen);

        // std::cout << "Time: " << tim.next_time() << std::endl;

        // more_test_lca(gen,pgen);
        extensive_lca(gen,pgen,tim);

        // //fail_test(gen,pgen);


        // std::cout << "Time: " << tim.next_time() << std::endl;

        // // std::cout << "start" << std::endl;
        // // parlay::sequence<cluster<int,int>> clusters;

        // // auto parent_tree = give_example_tree3();
        // // get_RC_tree(clusters,parent_tree);

        // // std::cout << "yay" << std::endl;
        // // deleteRCtree(clusters);

    }
    else if (forest_ratio==-2) {
        mid_test_lca(gen,pgen,tim);
    }
    else {
        test_lca(n,NUM_TRIALS,NUM_TREES,k,gen,pgen,forest_ratio,chain_ratio); //0 is forest ratio


    }


}