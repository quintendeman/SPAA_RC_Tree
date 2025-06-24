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
#include "LCA_test.h"

//tree generation, functions
#include "random_trees.h"
#include "treePrimitives.h"

//ternarizer and tree generator
#include "ternarizer.h"
#include "RCdynamic.h"
#include "treeGen.h"

const bool PRINT_B = false; //progress prints
const bool PRINT_T = false; //timer prints

template<typename T, typename D>
void get_RC_tree(parlay::sequence<cluster<T,D>>& clusters,  parlay::sequence<T>& parent_tree, bool extra_print=false) {
    parlay::sequence<parlay::sequence<T>> G;
    G=convert_parents_to_graph(G,parent_tree);
    if (extra_print) std::cout << "made graph" << std::endl;
    create_base_clusters<T,T>(G, clusters, 3); //max deg=3; each vertex in the tree has maximum degree 3
    if (extra_print) std::cout << "made base clusters" << std::endl;
    create_RC_tree(clusters, static_cast<T>(parent_tree.size()), static_cast<T>(0), [] (int A, int B) {return A+B;}, false); //randomized=false -- use the deterministic contraction method //defretval=0 -- default cluster val (ex -inf, inf, 0) //print=false -- don't print contraction sizes
    if (extra_print) std::cout << "made tree " << std::endl;
}


template<typename T, typename D>
void error_print(parlay::sequence<cluster<T,D>*>& answers, parlay::sequence<T>& real_answers, int iter, int iter2, parlay::sequence<T>& parent_tree, int k, double forest_ratio, double chain_ratio, std::string error_file_prefix, parlay::sequence<std::tuple<T,T,T>>& queries, int i) {
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
template<typename T, typename D>
void handle_answers(parlay::sequence<std::tuple<T,T,T>>& queries, parlay::sequence<cluster<T,D>*>& answers, int k, parlay::sequence<T>& parent_tree, parlay::sequence<cluster<T,D>>& clusters, int iter, int iter2, double forest_ratio, double chain_ratio, std::string error_file_prefix="error_log") {

    parlay::sequence<T> real_answers = parlay::tabulate(k,[&] (size_t i) {
        T u = std::get<0>(queries[i]);
        T v = std::get<1>(queries[i]);
        T rprime = std::get<2>(queries[i]);
        return unrooted_forest_lca(parent_tree,u,v,rprime);
    });
    //do checks in parallel, this is bottleneck?
    bool any_bad = parlay::reduce(parlay::delayed_tabulate(k,[&] (size_t i) {
        T u = std::get<0>(queries[i]);
        T v = std::get<1>(queries[i]);
        T rprime = std::get<2>(queries[i]);
        if ((answers[i]==nullptr && real_answers[i] != -1) || (answers[i] != nullptr && answers[i]->index != real_answers[i])) {
            return true; //set that error does exist

        }
        return false;

    }),MonoidOr()) ;

    //std::cout << "<";
    if (any_bad) { //if there was a bad output, go through sequentially to find error
        for (int i  = 0; i < k; i++) {
            T u = std::get<0>(queries[i]);
            T v = std::get<1>(queries[i]);
            T rprime = std::get<2>(queries[i]);
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
template<typename T, typename D>
void run_from_file(std::string file) {
    std::ifstream my_file;
    my_file.open(file);
    std::string first_line;
    std::getline(my_file,first_line);
    std::cout << "first line: " << first_line << std::endl;
    int n = 0;
    my_file >> n;
    int temp;
    parlay::sequence<T> parent_tree(n);
    for (int i = 0; i < n; i++) {
        my_file >> temp;
        parent_tree[i]=temp;
    }
    int k = 0;
    my_file >> k;
    k=1; //TOD2* change back
    parlay::sequence<std::tuple<T,T,T>> queries(k);
    T t1, t2, t3 = 0;
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
    parlay::sequence<cluster<T,D>> clusters;
    get_RC_tree(clusters,parent_tree,false);
    int NUM_TRIALS, iter, iter2 =-1; //because we chose the trial
    double forest_ratio, chain_ratio = -1; //because we already made tree

    parlay::sequence<cluster<T,D>*> answers(k);

    batchLCA(clusters,queries,answers);
    handle_answers(queries,answers,k,parent_tree,clusters,iter,iter2,forest_ratio,chain_ratio,"nofileprint"); //alt is alt error file to avoid covering previous print

    deleteRCtree(clusters);


}
template<typename T, typename D>
void run_tree(int NUM_TRIALS, parlay::sequence<cluster<T,D>>& clusters, int k,  parlay::sequence<T>& parent_tree, parlay::random_generator& pgen,  std::uniform_int_distribution<T>& dis, int iter, double forest_ratio, double chain_ratio) {

    parlay::internal::timer myt;
    
    //std::cout << "about to start LCA" << myt.next_time() << std::endl;

    dis(pgen); //stepping pgen

    for (int iter2 = 0; iter2 < NUM_TRIALS; iter2++) {
        //TOD2* test using variable batch sizes on same tree?
        parlay::sequence<std::tuple<T,T,T>> queries = parlay::tabulate(k,[&] (size_t i) {
            auto r = pgen[i];
            return std::make_tuple(dis(r),dis(r),dis(r));
        });

        //std::cout << "read queries " << myt.next_time() << std::endl;

        //NOTE! answers size initialized here*
        parlay::sequence<cluster<T,D>*> answers(k);

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
template<typename T, typename D>
void test_lca(int n, int NUM_TRIALS, int NUM_TREES, int k, std::mt19937& gen,parlay::random_generator& pgen, double forest_ratio, double chain_ratio, bool extra_print) {


    if (extra_print) std::cout << "starting n " << n << std::endl;
   // extra_print=false; //temporary

    std::uniform_int_distribution<T> dis(0,n-1);

    parlay::internal::timer t2;

    for (int iter = 0; iter < NUM_TREES; iter++) {

       if (PRINT_T) std::cout << "starting tree " << iter << " at time " << t2.next_time() << std::endl;
       // if (iter == 24) extra_print=true; //specific iter of interest
        if (extra_print) std::cout << "starting tree " << iter << std::endl;
        //generate random parent tree (permuted)
        parlay::sequence<T> parent_tree = generate_random_tree_perm_par<T>(n,pgen,chain_ratio); //TOD2* switch to random perm tree generate_random_tree(n,gen); // 
        ////generate_random_tree_perm(n,gen,chain_ratio);//

        divide_into_forests_par<T>(n,parent_tree,forest_ratio,pgen);
        //see the clusters
        if (extra_print) print_parent_tree(parent_tree,"parent tree:");

        if (PRINT_T) std::cout << "manual tree " << t2.next_time() << std::endl;

        //for debugging
        // parlay::sequence<parlay::sequence<int>> child_tree(parent_tree.size(),parlay::sequence<int>());
        // partree_to_childtree(parent_tree,child_tree);
        // print_child_tree(child_tree,"child tree");
        if (extra_print) std::cout << "about to make clusters" << std::endl;
        //create the RC tree
        parlay::sequence<cluster<T,D>> clusters;
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

    
using wedge = std::tuple<long,long,double>;

parlay::sequence<long> tree_gen(long graph_size, parlay::sequence<cluster<long, double>>& clusters, double mean, double ln) {


    const double min_weight = 0.0;
    const double max_weight = 100.0f;
    auto distribution = exponential; // either exponential, geometric, constant or uniform

    int II = 0;

    // std::cout << "II " << II << std::endl;
      TreeGen<long, double> TG(graph_size, min_weight, max_weight, ln, mean, distribution, true, 0);

      TG.generateInitialEdges();

      auto retedges = TG.getAllEdges();

      parlay::sequence<std::pair<long,long>> ret_edges_formatted = parlay::map(retedges,[&] (std::tuple<long,long,double>& edge) {return std::make_pair(std::get<0>(edge),std::get<1>(edge));});

      auto initial_parent_tree = parentTree_from_treeGen(graph_size,ret_edges_formatted);
      //pseq(initial_parent_tree,"initial par tree");

      // auto retedges = generate_high_degree_graphs(graph_size);
    
      ternarizer<long, double> TR(graph_size, 0);

      auto ret_edge_modified = TR.add_edges(retedges);
      
      // TR.print_state();
      TR.verify_simple_tree();

      const long max_degree = 3;
      double defretval = 0.0;

      parlay::sequence<wedge> empty_edges_sequence;
      create_base_clusters(clusters, ret_edge_modified.second, max_degree, graph_size * extra_tern_node_factor); //note the .second needed
      create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;},false);

      parlay::sequence<std::pair<long,long>> all_edges = parlay::map(ret_edge_modified.second,[&] (std::tuple<long,long,double>& edge) {return std::make_pair(std::get<0>(edge),std::get<1>(edge));});

      
      return parentTree_from_treeGen(initial_parent_tree.size()+2*retedges.size(),all_edges);

}


long tree_gen_size(long graph_size, parlay::sequence<cluster<long, double>>& clusters, double mean, double ln, std::string dist) {

    const double min_weight = 0.0;
    const double max_weight = 100.0f;
    auto distribution=exponential;// either exponential, geometric, constant or uniform
    if (dist=="e") {
        distribution=exponential;

    }
    else if (dist=="u") {
        distribution=uniform;
    }
    else {
        std::cout << "error, wrong dist, abort " << std::endl;
        exit(10);
    }

    int II = 0;

    // std::cout << "II " << II << std::endl;
      TreeGen<long, double> TG(graph_size, min_weight, max_weight, ln, mean, distribution, true,0); //passing in seed of 0
      //std::cout << "seed used " << TG.seed << std::endl;

      TG.generateInitialEdges();

      auto retedges = TG.getAllEdges();

    
      ternarizer<long, double> TR(graph_size, 0);

      auto ret_edge_modified = TR.add_edges(retedges);
      
      // TR.print_state();
      TR.verify_simple_tree();

      const long max_degree = 3;
      double defretval = 0.0;

      parlay::sequence<wedge> empty_edges_sequence;
      create_base_clusters(clusters, ret_edge_modified.second, max_degree, graph_size * extra_tern_node_factor);
      create_RC_tree(clusters, graph_size*extra_tern_node_factor, defretval, [] (double A, double B) {return A+B;},false); //Q* added extra_tern_node_factor

      //return graph_size + 2*retedges.size(); old ternarizer size
      return graph_size;//graph_size * extra_tern_node_factor; //new ternarizer size

}


std::chrono::duration<double> get_single_runtime(parlay::random_generator& pgen,parlay::sequence<cluster<long, double>>& clusters, int k,std::uniform_int_distribution<long>& dis,parlay::sequence<long>& parent_tree,bool extra_testing=false) {

    dis(pgen); //resettle random generator
    parlay::sequence<std::tuple<long,long,long>> queries = parlay::tabulate(k,[&] (size_t i) {
        auto r = pgen[i];
        return std::make_tuple(dis(r),dis(r),dis(r));
    });

    //std::cout << "read queries " << myt.next_time() << std::endl;

    //NOTE! answers size initialized here*
    parlay::sequence<cluster<long,double>*> answers(k);

    auto static_creation_start = std::chrono::high_resolution_clock::now(); //

    batchLCA(clusters,queries,answers);

    auto static_creation_end = std::chrono::high_resolution_clock::now();
    //TOD2* add assertion here that queries check out -> within handle_answers
    if (extra_testing) handle_answers(queries,answers,k,parent_tree,clusters,0,0,0,0,"nofileprint"); 

    return static_creation_end-static_creation_start;

}

parlay::sequence<long> get_kvals(long max_k, int kscale) {
    parlay::sequence<long> kopts(kscale,1);
    
    for (int i = 1; i < kscale; i++) {
        if (i % 2 == 1) {
            kopts[i]=kopts[i-1] * 5;
        }
        else {
            kopts[i]=kopts[i-1] * 2;
        }

    }

    parlay::sequence<long> kvals = parlay::filter(kopts,[&] (long kcand) {return kcand <= max_k;}); 
    return kvals;

}
// //force a different tree each time, according to the TG seed change
// void bench_diff_trees(parlay::random_generator& pgen,long graph_size, long max_k, int trials_per, double mean, double ln,std::string dist_choice) {
//     std::cout << "bench diff trees " << std::endl;
//     int kscale=30;
//     long k = -1;
//     parlay::sequence<long> kvals = get_kvals(max_k,kscale);

//     const double min_weight = 0.0;
//     const double max_weight = 100.0f;
//     auto distribution=exponential;// either exponential, geometric, constant or uniform
//     if (dist_choice=="e") distribution=exponential;
//     else if (dist_choice=="u") distribution=uniform;
//     else {
//         std::cout << "error, wrong dist, abort " << std::endl;
//         exit(10);
//     }

//     TreeGen<long, double> TG(graph_size, min_weight, max_weight, ln, mean, distribution, true,15213); // note seed is 15213 (which also is default val)

//     parlay::sequence<long> parent_tree(1,1); //dummy to fill in necessary argument

//     double runtime=0;
//     for (int j = 0; j < kvals.size(); j++) {

//         k=kvals[j];
//         for (int iter = 0; iter < trials_per; iter++) {

//             parlay::sequence<cluster<long, double>> clusters; 

//             TG.generateInitialEdges(); //increments seed twice

//             auto retedges = TG.getAllEdges();


//             ternarizer<long, double> TR(graph_size, 0);

//             auto ret_edge_modified = TR.add_edges(retedges);

//             TR.verify_simple_tree();

//             const long max_degree = 3;
//             double defretval = 0.0;

//             parlay::sequence<wedge> empty_edges_sequence;
//             create_base_clusters(clusters, ret_edge_modified, max_degree, graph_size * extra_tern_node_factor);
//             create_RC_tree(clusters, graph_size, defretval, [] (double A, double B) {return A+B;},false);


//             long n = graph_size + 2*retedges.size(); //ternarized n

//             std::uniform_int_distribution<long> dis(0,n-1);
//             runtime= get_single_runtime(pgen,clusters,k,dis,parent_tree,false).count();
//             std::cout << parlay::internal::init_num_workers() << "," << n << "," << graph_size << "," << k << ", " << ln << "," << mean << "," << dist_choice <<  ", " << runtime << std::endl;

//             deleteRCtree(clusters);
//         }

        
//     }

// }



void bench(parlay::random_generator& pgen,long oldn, long max_k, int trials_per, double mean, double ln,std::string dist_choice) {
    int kscale=30;
   
    long k = -1;
    parlay::sequence<long> kopts(kscale,1);
    
    for (int i = 1; i < kscale; i++) {
        if (i % 2 == 1) {
            kopts[i]=kopts[i-1] * 5;
        }
        else {
            kopts[i]=kopts[i-1] * 2;
        }

    }

    parlay::sequence<long> kvals = parlay::filter(kopts,[&] (long kcand) {return kcand <= max_k;}); 
    parlay::sequence<cluster<long, double>> clusters; 

    parlay::sequence<long> parent_tree(1,1);
    //one tree for all testing
    tree_gen_size(oldn,clusters,mean,ln,dist_choice); 

    std::uniform_int_distribution<long> dis(0,oldn-1); //vary among the old value of n (only real vertices)

    double runtime=0;
    for (int j = 0; j < kvals.size(); j++) {
        k=kvals[j];
        for (int iter = 0; iter < trials_per; iter++) {
            runtime= get_single_runtime(pgen,clusters,k,dis,parent_tree,false).count();
            std::cout << parlay::internal::init_num_workers() << "," << oldn << "," << k << ", " << ln << "," << mean << "," << dist_choice <<  ", " << runtime << std::endl;
        }

        
    }
   

    deleteRCtree(clusters);



}


//pass in n,k via input parms
void bench_threads(parlay::random_generator& pgen, long oldn, long k, int NUM_TRIALS, double mean, double ln,std::string dist_choice) {

  
    auto distribution = exponential; // either exponential, geometric, constant or linear

    parlay::sequence<cluster<long, double>> clusters; 

    parlay::sequence<long> parent_tree(1,1); //not used, because handle answers set to false currently. If true, need to get parent tree from clusters
    //one tree for all testing
    //n is ternerized size
    tree_gen_size(oldn,clusters,mean,ln,dist_choice); 
   

    std::uniform_int_distribution<long> dis(0,oldn-1);

    double runtime=0;
    for (int iter = 0; iter < NUM_TRIALS; iter++) {
        runtime = get_single_runtime(pgen,clusters,k,dis,parent_tree,false).count();
        std::cout << parlay::internal::init_num_workers() << "," << oldn << "," << k << ", " << ln << "," << mean << ", " << dist_choice << ", " << runtime << std::endl;
        
    }
    deleteRCtree(clusters);


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

void test() {
    int n = 5;
    parlay::sequence<cluster<long, double>> clusters; 
    double mean = 2;
    double ln = 0.1;

    parlay::sequence<long> parent_tree = tree_gen(n,clusters,mean,ln); //one tree for all testing

    pseq(parent_tree,"parent tree");

    deleteRCtree(clusters);

}

int main(int argc, char* argv[]) {

    parlay::internal::timer tim = parlay::internal::timer();

    //std::cout << "hello world!31" << std::endl;

    // std::mt19937 my_gen(10);
    // std::uniform_int_distribution<int> hi(0,3);
    // for (int i = 0; i < 20; i++) std::cout << hi(my_gen) << " ";
    // std::cout << std::endl;

    int n = 10; //defaults
    int NUM_TRIALS = 1;
    int NUM_TREES = 1;
    int k = 1;
    int seed = 42; //fixed seed for easier testing
    int pseed = 43;
    double forest_ratio = 5.0/n;//1.41 / n; //ratio scales with n to not disconnect too much
    double chain_ratio = .3;
    bool is_from_file = false;
    double mean=8;
    double ln=.1;
    std::string filename = "";
    std::string dist_choice = "";


    parse_input(argc,argv,n,NUM_TRIALS,seed,pseed,NUM_TREES,k,forest_ratio,chain_ratio,is_from_file,filename,mean,ln,dist_choice);

    //std::cout << "pseed: " << pseed << std::endl;

    std::mt19937 gen(seed);
    parlay::random_generator pgen(pseed);

    // std::cout << "Batch size: " << k << std::endl;
    // std::cout << "for ratio: " << forest_ratio << std::endl;
    // std::cout << "chain ratio: " << chain_ratio << std::endl;

    // std::cout << "about to start" << tim.next_time() << std::endl;


    if (is_from_file) {
        std::cout << "running from file " << filename << std::endl;
        run_from_file<int,int>(filename);
    }

    //forr serves as flag for different run modes
    //run_from_file is true: load tree, query from file
    //-forr -1: extensive testing
    //-forr -2: some testing
    //-forr -4: bench on k
    //-forr -5: bench on threads
    //otherwise: run a small lca example (based on other input parameters)
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
        mid_test_lca<int,int>(gen,pgen);
    }
    
    else if (forest_ratio==-4) {
        bench(pgen,n,k,NUM_TRIALS,mean,ln,dist_choice);

    }
    else if (forest_ratio==-5) {
        //tim.next_time();
        //std::cout << "num_threads,n,ternsize,ln,mean,dist,time" << std::endl;
        bench_threads(pgen,n,k,NUM_TRIALS,mean,ln,dist_choice);
        //std::cout << "total runtime of tests: " << tim.next_time() << std::endl;
    }
    else if (forest_ratio==-6) {
        test();
    }
    // else if (forest_ratio==-7) {
    //     bench_diff_trees(pgen,n,k,NUM_TRIALS,mean,ln,dist_choice);
    // }
    else {
        test_lca<int,int>(n,NUM_TRIALS,NUM_TREES,k,gen,pgen,forest_ratio,chain_ratio); //0 is forest ratio


    }


}