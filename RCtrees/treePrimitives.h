//tree functions  not pertaining to randomness

//file for creating random trees (of bounded degree)

#ifndef TREE_PRIM_H
#define TREE_PRIM_H

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

template<typename T>
void make_directed_forest(parlay::sequence<parlay::sequence<T>>& unrooted_tree, parlay::sequence<T>& parent_tree) {
    for (int i = 0; i < parent_tree.size(); i++) {
        if (parent_tree[i]==-1) { //if unconnected portion, direct component with that node as root
            make_directed_tree(unrooted_tree,parent_tree,static_cast<T>(i));
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
         printf("(%d,%d) ", j , static_cast<int>(tree[j]));

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
    Converts the parents array into a symmetric graph
*/
template <typename graph, typename T>
graph convert_parents_to_graph(graph G, parlay::sequence<T>& parents)
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

//given a list of edges that form a (unbounded degree) tree, convert to parent tree form
//Sequential -- don't run often!
//TOD2* What data type for newsize? 
//slightly inefficient converting to long then back to T, should maintain in type T
template<typename T>
parlay::sequence<T> parentTree_from_treeGen(long newsize, parlay::sequence<std::pair<T,T>>& edges, bool extra_print=false) {


    //std::cout << "newsize: " << newsize << std::endl;
   

    if (extra_print) {
        std::cout << "printing edges" << edges.size() << std::endl;
        for (int i = 0; i < edges.size(); i++) {
            std::cout << edges[i].first << " " << edges[i].second << std::endl;
        }
    }   
    parlay::sequence<std::pair<T,T>> rev_edges = parlay::tabulate(edges.size(),[&] (size_t i) {return std::make_pair(edges[i].second,edges[i].first); } );
    if (extra_print) {
        std::cout << "printing rev edges " << rev_edges.size() << std::endl;
        for (int i = 0; i < rev_edges.size(); i++) {
            std::cout << rev_edges[i].first << " " << rev_edges[i].second << std::endl;
        }
    }
    edges.append(rev_edges);
    auto all_edges = parlay::remove_duplicates(edges);
    if (extra_print) {
        std::cout << "printing all edges " << all_edges.size() << std::endl;
        for (int i = 0; i < edges.size(); i++) {
            std::cout << all_edges[i].first << " " << all_edges[i].second << std::endl;
        }
    }
    auto edge_groups = parlay::group_by_key(all_edges);
    if (extra_print) {
        std::cout << "printing edge groups " << edge_groups.size() << std::endl;
        for (int i = 0; i < edge_groups.size(); i++) {
            std::cout << edge_groups[i].first << ":" << std::endl;
            pseq(edge_groups[i].second,"connects");
        }
    }
    //capture the range of ids that the (ternerized) edges have
    int edge_size = *parlay::max_element(parlay::map(edge_groups,[&] (auto i) {return i.first;}))+1;

    parlay::sequence<parlay::sequence<long>> unrooted_tree(newsize);
    parlay::parallel_for(0,edge_groups.size(),[&] (size_t i) {
        unrooted_tree[edge_groups[i].first]=parlay::tabulate(edge_groups[i].second.size(),[&] (size_t j) {return edge_groups[i].second[j];});
    });
    if (extra_print) {
        std::cout << "printing unrooted tree" << std::endl;
        for (int i = 0; i < unrooted_tree.size(); i++) {
            // I am commenting since this seemed to be giving an error but its hopyfully nothing since this is just a print
            // pseq(unrooted_tree[i],"neighbors of " +std::to_string(i));
        }
    }


    parlay::sequence<T> parent_tree(unrooted_tree.size(),-1);

     //TOD2* forest? make directed tree doesn't work on forests*
    make_directed_forest(unrooted_tree,parent_tree);

    if (extra_print) pseq(parent_tree,"parent tree treeGen");
    return parent_tree;

}

#endif