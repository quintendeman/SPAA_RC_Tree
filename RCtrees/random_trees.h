//file for creating random trees (of bounded degree)

#ifndef RANTREE
#define RANTREE

#include <atomic>
#include <random>
#include <set>
#include <algorithm>
#include <iostream>
#include <mutex>
#include <parlay/alloc.h>
#include "cluster.h"
#include "../examples/counting_sort.h"
#include <parlay/random.h>

//given a tree written in terms of parents (index i has the edge i, tree[i]; tree[i] is the parent of i), reformat in terms of children (index i has a list (j_1,...j_m), where (i,j_x) are edges, and (j_x) are the children of i)
//in the original parents tree, tree, note that i=tree[i] iff i is the root.
//parlay::sequence<parlay::sequence<T>> child_tree(tree.size(),parlay::sequence<T>());
template<typename T>
parlay::sequence<parlay::sequence<T>> partree_to_childtree(parlay::sequence<T>& tree, parlay::sequence<parlay::sequence<T>>& child_tree) {
    
    for (int i = 0; i < tree.size(); i++) {
        if (i == tree[i]) continue;
        //TOD2 keep this sorted maybe? Note that sorts ascending currently.
        child_tree[tree[i]].push_back(i);
    }
    return child_tree;

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
//TOD2* layer with random permutation (for more rigorous testing, affected by parallel scheduling)
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
template<typename T>
parlay::sequence<T> generate_random_tree(T num_elements, std::mt19937& gen) {
    assert(num_elements > 0);
    parlay::sequence<T> parents = parlay::tabulate(num_elements,[&] (T v) {return (T) 0;});

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
// //randomly permute elements 0 through n-1, in parallel
// template<typename T>
// parlay::sequence<T> generate_random_perm(T n, parlay::random_generator& pgen) {
//     parlay::sequence<T> perm(n,-1);
//     parlay::sequence<T> cand = parlay::iota(n);
//     parlay::sequence<T> pick;
//     std::uniform_real_distribution<double> dis(0, 1);
//     parlay::sequence<std::atomic<int>> counters;
//     parlay::sequence<T> io = parlay::iota(n);

//     int base = 1; //base case size (1 means no special base case/low value case)
//     //expected ~ O(log n) iters to terminate
//     while (cand.size() >= base) {
//         io = parlay::iota(cand.size());
//         counters = parlay::tabulate(cand.size(),[&] (size_t i) {return std::atomic<int>(0);});
//         //pick random elt within cand
//         pick = parlay::tabulate(cand.size(),[&] (size_t i) {
//             return dis(pgen[i]) % cand.size();
//         });
//         parlay::parallel_for(0,pick.size(),[&] (size_t i) {
//             //if we were the first to pick a certain #, add to perm
//             if (counters[pick[i]].fetch_add(1)==0) {
//                 perm[cand[i]]=perm[cand[pick[i]]];

//             }

//         });
//         //we continue processing elements whose perm has not been decided yet
//         cand = parlay::filter(io,[&] (T index) {
//             return perm[index]==-1;
//         });



//     }

//     //debugging
//     auto results = histogram_by_key(perm);
//     if (results.size() < n) {
//         std::cout << "error, not all #s represented in perm" << std::endl;
//         for (int i = 0; i < n; i++) {
//             std::cout << perm[i] << " ";
//         }
//         std::cout << std::endl;
//         exit(3001);
//     }


// }


//generate_random_tree, but with a permutation added to enhance randomness (catch bugs based on children relative ordering)
template<typename T>
parlay::sequence<T> generate_random_tree_perm(T num_elements, std::mt19937& gen) {
    parlay::sequence<T> parents = generate_random_tree(num_elements,gen); //get parent tree, unpermuted

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



//find depth of tree (in child form)
//trimmed version of set_level_par
template<typename T>
int tree_depth(parlay::sequence<parlay::sequence<T>>& child_tree, T root1) {
    int depth = 0;

    parlay::sequence<T> stack; //stack from which we draw tasks
    stack.push_back(root1);

    parlay::sequence<parlay::sequence<T>> new_stack; 


    while (stack.size() > 0) {
        depth += 1;

        new_stack = parlay::sequence(stack.size(), parlay::sequence<T>());
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            auto s = stack[i];
            
            for (T j = 0; j < child_tree[s].size(); j++) new_stack[i].push_back(child_tree[s][j]);

        });

        //copy over new stack, but flattened for easier access
        stack = parlay::flatten(new_stack);
       
    }
    return depth;

}
template<typename T>
bool same_degrees(parlay::sequence<parlay::sequence<T>>& childtree1, parlay::sequence<parlay::sequence<T>>& childtree2) {
    T n = childtree1.size();

    auto rangn = parlay::iota(n);

      //confirm same # of leaves, degree 1 nodes, degree 2 nodes
    auto num_child_dist1 = parlay::histogram_by_key(parlay::map(rangn,[&] (int i) {return childtree1[i].size();}));
    auto num_child_dist2 = parlay::histogram_by_key(parlay::map(rangn,[&] (int i) {return childtree2[i].size();}));
    // for (int i = 0; i < num_child_dist1.size(); i++)
    // std::cout << num_child_dist1[i].first << " " << num_child_dist1[i].second << std::endl;
    if (num_child_dist1.size() != num_child_dist2.size()) {
        return false;
    }
    for (int i = 0; i < num_child_dist1.size(); i++) {
        int key = num_child_dist1[i].first;
        bool found_match = false;
        for (int j = 0; j < num_child_dist2.size(); j++) {
            if (key == num_child_dist2[j].first) {
                found_match=true;
                if (num_child_dist1[i].second != num_child_dist2[j].second) {
                    return false;

                }
            }
        }
        if (!found_match) {
            return false;
        }
    }
    return true;

}

//checking if 2 rooted trees are (possibly) equal by checking common properties
template<typename T>
bool possibly_equal(parlay::sequence<parlay::sequence<T>>& childtree1, parlay::sequence<parlay::sequence<T>>& childtree2, T root1, T root2) {

   
    //iso trees must have same # of nodes
    if (childtree1.size() != childtree2.size()) {
        return false;
    }
    //iso trees must have same depth
    if (tree_depth(childtree1,root1) != tree_depth(childtree2,root2)) {
        return false;
    }
    return same_degrees(childtree1,childtree2);

}

#endif