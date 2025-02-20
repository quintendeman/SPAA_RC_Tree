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


//note that we are NOT passing child_tree by reference, to get a new copy as the undirected graph
template<typename T>
parlay::sequence<parlay::sequence<T>> make_undirected_tree(parlay::sequence<parlay::sequence<T>> child_tree, parlay::sequence<T>& parent_tree) {
    for (int i = 0; i < child_tree.size(); i++) {
        child_tree[i].push_back(parent_tree[i]); //just push back the parent
    }
    return child_tree;

}

//given an unrooted tree, make it into a parent tree rooted at rprime
template<typename T>
void make_directed_tree(parlay::sequence<parlay::sequence<T>>& unrooted_tree, parlay::sequence<T>& parent_tree, T rprime) {

    std::deque<T> stack;
    parent_tree[rprime]=rprime;
    stack.push_back(rprime); 

    //int count = 0;
    while (stack.size() > 0) {
        // count += 1;
        // if (count > 50) exit(700); //for printing
        
        auto s = stack.back();
        stack.pop_back();

        // print_parent_tree(parent_tree,"intermediate directed tree");
        // std::cout << "s: " << s << std::endl;
        
        //the parent of each other edge in s is s
        for (int i = 0; i < unrooted_tree[s].size(); i++) {
            T child = unrooted_tree[s][i];
            //issue: The old root has itself as an edge
            if (child != parent_tree[s] &&  s != child) {//don't count the back edge
                parent_tree[child]=s;
                stack.push_back(child);
            }
        
        }
    }

}

//given a tree written in terms of parents (index i has the edge i, tree[i]; tree[i] is the parent of i), reformat in terms of children (index i has a list (j_1,...j_m), where (i,j_x) are edges, and (j_x) are the children of i)
//in the original parents tree, tree, note that i=tree[i] iff i is the root.
//parlay::sequence<parlay::sequence<T>> child_tree(tree.size(),parlay::sequence<T>());
template<typename T>
parlay::sequence<parlay::sequence<T>> partree_to_childtree(parlay::sequence<T>& tree, parlay::sequence<parlay::sequence<T>>& child_tree) {
    
    for (int i = 0; i < tree.size(); i++) {
        if (i == tree[i]) continue; //don't add the parent as a child
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


//find the root of a tree in a forest containing node u
//root = where the parent is itself (by our construction)
template<typename T>
T get_root(parlay::sequence<T>& tree, T u) {
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
    parlay::sequence<T> extended = parlay::tabulate(2*n+1,[&] (T i) {return (i-1) >> 1; });
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

//given parents graph (single tree), divide into ratio*n trees (for forest testing)
template<typename T>
void divide_into_forests(T n, parlay::sequence<T>& parents, double ratio,std::mt19937& gen) {

    std::uniform_real_distribution<double> dis(0,1);
    for (int i = 0; i < n; i++) {
        if (dis(gen) < ratio) {
            parents[i]=i; //cut off parent edge, now is root
        }
    }

}


//given parents graph (single tree), divide into ratio*n trees (for forest testing) in parallel
template<typename T>
void divide_into_forests_par(T n, parlay::sequence<T>& parents, double ratio,parlay::random_generator& pgen) {

    std::uniform_real_distribution<double> dis(0,1);
    dis(pgen);

    parlay::parallel_for(0,n,[&] (size_t i) {
        auto r = pgen[i];
        if (dis(r) < ratio) {
            parents[i]=i; //cut off parent edge, now is root
        }
        
    });

}

//count how many trees are in a forest
template<typename T>
T num_roots(parlay::sequence<T>& parents) {
    T counter = 0;
    for (int i = 0; i < parents.size(); i++) {
        if (parents[i]==i) counter += 1;
    }
    return counter;
}

//function in progress TOD2* finish or delete
// template<typename T>
// void print_forest_sizes(parlay::sequence<T>& parents) {
//     parlay::sequence<T> roots;
//     parlay::sequence<parlay::sequence<T>> child_tree;
//     partree_to_childtree(parents,child_tree);
//     for (int i = 0; i < parents.size(); i++) {
//         if (parents[i]==i) {
//             roots.push_back(i);
//         }
//     }
//     parlay::sequence<T> counters(roots.size(),1);

//     parlay::parallel_for(0,roots.size(),[&] (T i) {
//         for (int i = 0; i < )
        

//     })
// }


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



/*
    Generate a simple, single rooted graph with each node having two children
    Then, randomly, change the parent of each node with a certain probability such that it picks something on the left of it

    Returns an array of parents such that the parents of index V would be parents[V]
*/
template <typename T>
parlay::sequence<T> generate_tree_graph(T num_elements, int seed=time(0))
{
    assert(num_elements > 0);

    parlay::sequence<T> dummy_initial_parents = parlay::tabulate(num_elements, [&] (T v) {return (T) 0;});

    parlay::random_generator gen(seed); //default seed is time(0)
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
            std::uniform_int_distribution<> disint(0, v-1);

            dummy_initial_parents[v] = disint(gen);
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

/*
    Converts the parents array into a symmetric graph
*/
template <typename graph, typename T>
graph convert_parents_to_graph(graph G, parlay::sequence<T> parents)
{
    parlay::sequence<T> vertices = parlay::tabulate(parents.size(), [&] (T v) {return v;});

    G = parlay::map(vertices, [&] (T v) {
        if(parents[v] == v) // root
        {
            parlay::sequence<T> empty;
            return empty;
        }
        parlay::sequence<T> temp = parlay::tabulate(1, [&] (T i) {return i;});
        temp[0] = parents[v];
        return temp;
    });

    G = graph_utils<T>::symmetrize(G);

    return G;
}

/**
 * Ensures that the degree is capped for the parents vector
 * i.e. not too many nodes have the same parent
 * Uses locking! Fortunately, we don't care about the performance for the graph generation too much
*/
template <typename T>
parlay::sequence<T> degree_cap_parents(parlay::sequence<T> &parents, const T max_degree)
{
    parlay::sequence<std::atomic<T>> counts = parlay::tabulate(parents.size(), [] (size_t) {
       return std::atomic<T>(0); // Initialize each element with the value 0
    });


    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(v == parents[v])
            return;
        T parent_count = counts[parents[v]].fetch_add(1);
        if(parent_count < (max_degree - 1))
        {
            return;
        }
        else
            parents[v] = v;
    });

    return parlay::tabulate(parents.size(), [&] (T i) {
        return counts[i].load();
    });

}

template<typename T, typename D>
void degree_cap_add_edge(parlay::sequence<T> &parents, const T max_degree, parlay::sequence<std::tuple<T, T, D> >& tuples)
{
    parlay::sequence<std::atomic<T>> counts = parlay::tabulate(parents.size(), [] (size_t) {
       return std::atomic<T>(0); // Initialize each element with the value 0
    });

    std::vector<std::mutex> mutexes(counts.size());

    parlay::parallel_for(0, counts.size(), [&] (T chld) {
        if(parents[chld] != chld)
        {
            counts[parents[chld]].fetch_add(1);
        }
    });

    tuples = parlay::filter(tuples, [&] (auto tple) {
        auto& child = std::get<0>(tple);
        auto& parent = std::get<1>(tple);
        if(child == parent)
            return false;
        bool ret_val = false;

        mutexes[parent].lock();
        mutexes[child].lock();
        if(counts[child] < max_neighbours)
        {
            if(counts[parent] < max_neighbours - 2)
            {
                counts[parent]++;
                parents[child] = parent;
                ret_val = true;
            }
            else if(counts[parent] == max_neighbours - 1 && parents[parent] == parent)
            {
                counts[parent]++;
                parents[child] = parent;
                ret_val = true;        
            }
        }
        mutexes[parent].unlock();
        mutexes[child].unlock();
        return ret_val;
    });
    

    return;
}

#endif