//file for creating random trees (of bounded degree)

#ifndef RANTREE
#define RANTREE

#include <atomic>
#include <random>
#include <iostream>
#include <mutex>
#include <parlay/alloc.h>
#include "cluster.h"
#include "../examples/counting_sort.h"
#include <parlay/random.h>
#include "../examples/helper/graph_utils.h"
#include "parlay/primitives.h"
#include "treePrimitives.h"

std::mt19937 get_rand_gen(int seed) {
    if (seed == -1) {
        std::mt19937 gen(std::random_device{}());
        return gen;
    }
    else {
        std::mt19937 gen(seed);
        return gen;
    }
}


//another take on generate_tree_graph
//at each child, generate one or two children (naturally ternerized)
//root is 0
//TOD2* test on trees with greater than 2 children? 
template<typename T>
parlay::sequence<T> generate_random_tree(T num_elements, int seed=-1) {
    assert(num_elements > 0);
    parlay::sequence<T> parents = parlay::tabulate(num_elements,[&] (T v) {return (T) 0;});

    std::mt19937 gen = get_rand_gen(seed);

    std::uniform_real_distribution<double> dis(0, 1);

    T c = 1; //count # of elements already added to tree
    T par = 0; //the current parent, to which we add its children
    double p2 = .3; //probability of 2 children
    while (c < num_elements) {
        parents[c]=par;
        auto random_val = dis(gen);
        if (random_val < p2 && c+1 < num_elements) {
            parents[c+1]=par; 
            c += 2;
        }
        else {
            c += 1;
        }
        par += 1;

    }

    return parents;

}

//another take on generate_tree_graph
//at each child, generate one or two children (naturally ternerized)
//root is 0
//chain ratio is % of nodes with 1 child (as compared to 2 children) (with exception of endcap)
template<typename T>
parlay::sequence<T> generate_random_tree(T num_elements, std::mt19937& gen, double chain_ratio) {
    assert(num_elements > 0);
    parlay::sequence<T> parents = parlay::tabulate(num_elements,[&] (T v) {return (T) 0;});

    std::uniform_real_distribution<double> dis(0, 1);

    T c = 1; //count # of elements already added to tree
    T par = 0; //the current parent, to which we add its children
    while (c < num_elements) {
        parents[c]=par;
        auto random_val = dis(gen);
        if (random_val > chain_ratio && c+1 < num_elements) {
            parents[c+1]=par; 
            c += 2;
        }
        else {
            c += 1;
        }
        par += 1;

    }
    return parents;

}

//generate random tree of degree <= 3 in parallel
//at each child, generate one or two children 
//root is 0
//chain ratio is % of nodes with 1 child (as compared to 2 children) (with exception of endcap)
template<typename T>
parlay::sequence<T> generate_random_tree_par(T n, parlay::random_generator& pgen, double chain_ratio) {
    assert(n > 0);

    std::uniform_real_distribution<double> dis(0,1);
    dis(pgen); //resettle random generator
    parlay::sequence<double> random_bits(2*n+1,-1);

    //if the 1st n nodes all had 2 children, we would have 2*n+1 node tree
    parlay::sequence<T> extended = parlay::tabulate(2*n+1,[&] (T i) {return static_cast<T>((i-1) >> 1); });
    extended[0]=0; //0 is root
    //pseq(extended,"printing extended");
    parlay::sequence<bool> packbools = parlay::tabulate(2*n+1,[&] (T i) {
        if (i == 0 || i % 2 == 1) return true; //always keep at least one child
        auto r = pgen[i];
        double coin = dis(r);
        random_bits[i]=coin;
        if (coin > chain_ratio) return true;
        return false;
    });
    //pseq(packbools,"pack bools");
    //pseq(random_bits,"rand bits");
    //only take true positions
    auto filtered = parlay::pack(extended,packbools);
    //return the n-prefix of filtered
    return filtered.subseq(0,n);

}




//sequential way to create permutation. Given that parent tree creation is sequential, have permutation sequential too? 
template<typename T>
parlay::sequence<T> generate_random_perm_seq(T n, std::mt19937& gen) {
    parlay::sequence<T> perm(n,-1);
    parlay::sequence<T> cand = parlay::tabulate(n,[&] (T i) {return i;});
    std::uniform_int_distribution<int> dis(0, n-1); //note that we are slightly weighting lower #s to appear in perm
    T start_target = 0;
    T cur_target = start_target; //start with 0
    T choice_target;
    T choice_index;
    for (int counter = 0; counter < n; counter++) {
        
        choice_index = dis(gen) % cand.size();
        choice_target = cand[choice_index];
        //std::cout << cur_target << " " << choice_target << std::endl;
        perm[cur_target]=choice_target;
        cand.erase(cand.begin() + choice_index);
        if (cur_target != choice_target && choice_target != start_target) {
            cur_target=choice_target;
        }
        else {
            cur_target=cand[0];
            start_target = cand[0];
        }
    }

    // std::cout << "printing perm" << std::endl;
    // for (int i = 0; i < perm.size(); i++) {
    //     std::cout << i << ": " << perm[i] << std::endl;
    // }
    // std::cout << std::endl;
    return perm;

}


//randomly permute elements 0 through n-1, in parallel
//TOD2* use more delayed sequences instead of tabulates (here and throughout code)
//TOD2* neaten, slightly spaghetti
//use parlay::random_shuffle instead
//note: span not polylog due to concurrent writes? (use random shuffle instead)
template<typename T>
parlay::sequence<T> generate_random_perm_par(T n, parlay::random_generator& pgen) {
    parlay::sequence<T> perm(n,-1); //i maps to perm[i]
    parlay::sequence<T> cand = parlay::tabulate(n,[&] (T i) {return i;}); //set of #s that still lack a map (remaining domain)
    parlay::sequence<T> available_choices = parlay::tabulate(n,[&] (T i) {return i;}); //remaining range
    std::uniform_real_distribution<double> dis(0, 1);
    dis(pgen); //resettle generator

    //set counters to start at 0
    parlay::sequence<std::atomic<int>> counters = parlay::tabulate(cand.size(),[&] (size_t i) {return std::atomic<int>(0);});
   
    int base = 1; //base case size (1 means no special base case/low value case)
    //expected ~ O(log n) iters to terminate
    int counter = 0;
    while (cand.size() >= base) {
 
        //reset counters still in use
        parlay::parallel_for(0,cand.size(),[&] (size_t i) {
            counters[i]=0;
        });

     
        //pick random elt index within cand
        //if pick[i]=j, then cand[i] wants to map to choice[j]
        assert(available_choices.size() == cand.size());
        parlay::sequence<T> pick = parlay::tabulate(cand.size(),[&] (size_t i) {
            auto r = pgen[i];
            return static_cast<T>(floor(dis(r)*available_choices.size())); //random elt in 0 to cand.size()-1
        });

        //keep track of which picks go through
        //TOD2* change from atomic
        parlay::sequence<std::atomic<bool>> pick_remains(pick.size(),std::atomic<bool>(true));

        parlay::parallel_for(0,pick.size(),[&] (size_t i) {
            //if we were the first to pick a certain #, add to perm
            if (counters[pick[i]].fetch_add(1)==0) {
                perm[cand[i]]=available_choices[pick[i]];
                pick_remains[pick[i]]=false;

            }

        });

        //we continue processing elements whose perm has not been decided yet
        cand = parlay::filter(cand,[&] (T index) {
            return perm[index]==-1;
        });
        available_choices = parlay::pack(available_choices,pick_remains); 



    }

    //debugging
    auto results = histogram_by_key(perm);
    if (results.size() < n) {
        std::cout << "error, not all #s represented in perm" << std::endl;
        for (int i = 0; i < n; i++) {
            std::cout << perm[i] << " ";
        }
        std::cout << std::endl;
        exit(3001);
    }

    return perm;


}

//generate_random_tree, but with a permutation added to enhance randomness (catch bugs based on children relative ordering)
template<typename T>
parlay::sequence<T> generate_random_tree_perm(T num_elements, std::mt19937& gen, double chain_ratio) {
    parlay::sequence<T> parents = generate_random_tree(num_elements,gen,chain_ratio); //get parent tree, unpermuted

    parlay::sequence<T> perm = generate_random_perm_seq(num_elements,gen); //get permutation

    //print_parent_tree(parents,"par tree original");

    // std::cout << "printing perm" << std::endl;
    // for (int i = 0; i < perm.size(); i++) {
    //     std::cout << i << ": " << perm[i] << std::endl;
    // }
    // std::cout << std::endl;

    parlay::sequence<T> new_tree(num_elements,-1); 

    parlay::parallel_for(0,num_elements,[&] (size_t i) { //apply permutation
        new_tree[perm[i]]=perm[parents[i]];

    });

    //print_parent_tree(new_tree,"new tree");

    return new_tree; //return result of permutation

}


//generate_random_tree, but with a permutation added to enhance randomness (catch bugs based on children relative ordering)
template<typename T>
parlay::sequence<T> generate_random_tree_perm_par(T num_elements, parlay::random_generator& pgen, double chain_ratio) {
    parlay::sequence<T> parents = generate_random_tree_par(num_elements,pgen,chain_ratio); //get parent tree, unpermuted
    //TOD2* adjust
    std::uniform_int_distribution<size_t> dis; //range over all size_t vals
    size_t shuffle_seed = dis(pgen);
    auto vals = parlay::delayed_tabulate(num_elements,[&] (T i) {return i;});
    auto perm = random_shuffle(vals,parlay::random(shuffle_seed)); //get permutation

    // print_parent_tree(parents,"par tree original");

    // std::cout << "printing perm" << std::endl;
    // for (int i = 0; i < perm.size(); i++) {
    //     std::cout << i << ": " << perm[i] << std::endl;
    // }
    // std::cout << std::endl;

    parlay::sequence<T> new_tree(num_elements,-1); 

    parlay::parallel_for(0,num_elements,[&] (size_t i) { //apply permutation
        new_tree[perm[i]]=perm[parents[i]];

    });

    //print_parent_tree(new_tree,"new tree");

    return new_tree; //return result of permutation

}

//given a parent tree, permute it and return the new tree
template<typename T>
parlay::sequence<T> permute_tree(T num_elements, parlay::sequence<T>& parents, std::mt19937& gen) {

    parlay::sequence<T> perm = generate_random_perm_seq(num_elements,gen); //get permutation

    parlay::sequence<T> new_tree(num_elements,-1); 

    parlay::parallel_for(0,num_elements,[&] (size_t i) { //apply permutation
        new_tree[perm[i]]=perm[parents[i]];

    });
    return new_tree; //return result of permutation

}

/*
    Generate a simple, single rooted graph with each node having two children
    Then, randomly, change the parent of each node with a certain probability such that it picks something on the left of it

    Returns an array of parents such that the parents of index V would be parents[V]
*/
template <typename T>
parlay::sequence<T> generate_tree_graph(T num_elements)
{
    assert(num_elements > 0);

    parlay::sequence<T> dummy_initial_parents = parlay::tabulate(num_elements, [&] (T v) {return (T) 0;});

    std::mt19937 sgen(std::random_device{}());
    std::uniform_int_distribution<int> dis0(1,std::numeric_limits<int>::max()-1);
    parlay::random_generator gen(dis0(sgen)); //default seed is time(0)
    std::uniform_real_distribution<double> dis(0, 1);    

    parlay::parallel_for(0, num_elements, [&] (T v) {
        auto r = gen[v];
        double random_val = dis(r);

        static const double anywhere_left_weight = 0.5;
        static const double immediate_left_weight = 10.0;
        static const double root_weight = 0.01; /* warning, does not represet probability of a null cluster as degree capping may create more forests */

        static const double anywhere_prob = (anywhere_left_weight/(anywhere_left_weight+immediate_left_weight+root_weight));
        static const double root_prob = anywhere_prob + (root_weight/(anywhere_left_weight+immediate_left_weight+root_weight));

        if (random_val <= anywhere_prob && v > 0)
        {
            std::uniform_int_distribution<T> disint(0, v-1);

            dummy_initial_parents[v] = disint(r);
        }
        else if (random_val < root_prob)
        {
            dummy_initial_parents[v] = v;
        }
        else
            dummy_initial_parents[v] = v == 0 ? 0 : v - 1;
    });  
    
    return dummy_initial_parents;
}

#endif