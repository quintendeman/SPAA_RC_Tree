//Code to handle LCA queries

#ifndef LCAHH //! cannot define function of this name, good to know
#define LCAHH

//header copied from RC.h
#include "RC.h"
#include "cluster.h"


// T - index type
// D - type of data edges/clusters are storing
// Given clusters U and V, find the LCA in the original tree of the representative vertices of U and V. Store the cluster for which this is the representative vertex in ans. The LCA is oriented by the root of the RC tree
template<typename T, typename D>
void LCA(cluster<T,D>* U, cluster<T,D>* V, cluster<T,D>*& ans) {
    std::cout << "Fill in later" << std::endl;

}

//Given an unrooted tree, the RC tree contraction enforces an orientation/root (in some fashion) (consider the root of the unrooted tree to be the root of the RC tree) 
//Thus, given the RC tree, get out this structure
//n = # of vertices
template<typename T, typename D>
parlay::sequence<T> clusters_to_parents(T n, parlay::sequence<cluster<T,D>>& clusters) {
    parlay::sequence<T> parents(n,-1);
    //map the cluster index (its id) to a value in 0 to n-1 (because cluster indices range to n+m > n)
    std::unordered_map<T,T> normalized_asg; 
    T hash_counter = 0; //id counter for hash assignment

    for (int i = 0; i < clusters.size(); i++) { 
            T id = clusters[i].index; //cluster id
          
            //if we're at a root
            if (clusters[i].parent == nullptr) {
                parents[id]=id;
            }
            //if we're not at a root
            else {
                T parent_id = clusters[i].parent->index;
                parents[id]=parent_id;  
            }
    }

    return parents;

}

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

//returns true if, *in the RC tree*, u's cluster is a descendant of v's cluster (v is ancestor)
template<typename T, typename D>
bool is_descendant(parlay::sequence<cluster<T,D>>& clusters, T u, T v) {
    T uold = u;
    T unew = clusters[u].parent.index;
    if (u == v) return true;
    while (uold != unew) {
        if (unew == v) return true;
        uold = unew;
        unew = clusters[unew].parent.index;
    }
    return false;

}

#endif //LCA