//Code to handle LCA queries

#ifndef LCAHH //! cannot define function of this name, good to know
#define LCAHH

//header copied from RC.h
#include "RC.h"
#include "cluster.h"
#include "include/parlay_hash/unordered_map.h" //for parlayhash concurrent hash map
#include "static_par_LCA.h"

void find_nodes_involved(parlay::sequence<cluster<T,D>>& clusters, int k, parlay::sequence<cluster<T,D>*> Uvertices, parlay::sequence<cluster<T,D>>* Vvertices, parlay::sequence<cluster<T,D>*> Rvertices, cluster<T,D>* root, parlay::sequence<bool>& marked_clusters,parlay::sequence<T>& involved_nodes) {

    parlay::sequence<T> stack;
    parlay::sequence<T> new_stack;


    parlay::parallel_for(0,k,[&] (size_t i) {
        // marked_clusters[Uvertices[i]]=true;
        // marked_clusters[Vvertices[i]]=true;
        // marked_clusters[Rvertices[i]]=true;
        new_stack[i].push_back(Uvertices[i]);
        new_stack[i].push_back(Vvertices[i]);
        new_stack[i].push_back(Rvertices[i]);
    });

    //TOD2* is remove duplicates too expensive? 
    //TOD2* does remove duplicates automatically output an array of reduced size? (or are there empty spaces?)
    stack=parlay::remove_duplicates(parlay::flatten(new_stack));
    involved_nodes.append(stack); //parallel insert the active nodes of this level 

    //bottom-up, propagate true marks
    while(stack.size() > 0) {
        new_stack=parlay::sequence(stack.size(),-1);
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            auto s = stack[i];
            marked_clusters[s]=true;

            new_stack[i]=clusters[s].parent->index; //unlike last bottom up, doesn't get rid of duplicates


        });

        stack=parlay::remove_duplicates(parlay::flatten(new_stack));
        involved_nodes.append(stack); //parallel insert the active nodes of this level 


    }


}

T instack_choose_boundary(parlay::sequence<T>& involved_nodes, parlay::sequence<T>& closest_boundary, parlay::sequence<cluster<T,D>>& clusters, cluster<T,D>* root, parlay::sequence<parlay::sequence<T>>& alt_child_tree, parlay::sequence<T>& alt_parent_tree, parlay::parlay_unordered_map<T,T>& index_map) {

    cluster<T,D>* cur = &clusters[involved_nodes[s]];
    cluster<T,D>* parent = &clusters[involved_nodes[alt_parent_tree[s]]];

    cluster<T,D>* left;
    cluster<T,D>* right;
    D def_val;
    parent->find_boundary_vertices(left,def_val,right,def_val,def_val);

    cluster<T,D>* child_left;
    cluster<T,D>* child_right;
    cur->find_boundary_vertices(child_left,def_val,child_right,def_val,def_val);

    //if this cluster is a binary cluster
    if (parent->state & (binary_cluster | base_edge) ) {
        T par_closest = closest_boundary[index_map[parent->index]];
        if (par_closest == index_map[child_left->index] || par_closest==index_map[child_right->index]) {
            return par_closest;
        }
        else { //pick the boundary NOT shared by the parent
            if (left->index==child_left->index || right->index == child_left->index) {
                return child_right->index;
            } 
            else {
                return child_left->index;
            }
        }

    }
    //TOD2* is nullary cluster one of the involved nodes?
    else if (parent->state & (nullary_cluster)) {
        std::cout << "did not account for this abort" << std::endl;
        exit(802);
    }
    else { //is unary cluster
        return left->index; //closer boundary vertex is boundary vertex of parent
    }


}
//in involved_nodes,
//for each binary cluster, find which boundary vertex is closer to the root
//for each unary cluster will just return the only boundary vertex
//TOD2* does the find_boundary_vertices call within this function look outside the subset of the tree we are considering, breaking the work bound?
//top-down computation
template<typename T, typename D>
void find_closest_boundary(parlay::sequence<T>& involved_nodes, parlay::sequence<T>& closest_boundary, parlay::sequence<cluster<T,D>>& clusters, cluster<T,D>* root, parlay::sequence<parlay::sequence<T>>& alt_child_tree, parlay::sequence<T>& alt_parent_tree, parlay::parlay_unordered_map<T,T>& index_map) {
    parlay::sequence<T> stack;
    stack.push_back(index_map[root->index]);

    parlay::sequence<parlay::sequence<T>> new_stack;

    while (stack.size() > 0) {
        new_stack = parlay::sequence<T>(stack.size(),parlay::sequence<T>());
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            auto s = stack[i];
            closest_boundary[s] = instack_choose_boundary(involved_nodes,closest_boundary,clusters,root,alt_child_tree,alt_parent_tree,index_map);
            for (T j : alt_child_tree[i]) {
                new_stack[i].push_back(j);
            }
            
        });
        stack=parlay::flatten(new_stack);

    }


}

void create_map_trees(parlay::sequence<T>& involved_nodes,parlay::sequence<cluster<T,D>>& clusters, cluster<T,D>* root, parlay::sequence<parlay::sequence<T>>& alt_child_tree, parlay::sequence<T>& alt_parent_tree, parlay::parlay_unordered_map<T,T>& index_map) {

    //map the value in involved_nodes to its index; we will then create the tree using the index
    //TOD2* is there a bulk insert command?
    parlay::parallel_for(0,involved_nodes.size(),[&] (size_t i) {
        index_map.insert(involved_nodes[i],i);
    });

  
    alt_child_tree = parlay::tabulate(involved_nodes.size(),[&] (size_t u) {
        T node = involved_nodes[u];
        parlay::sequence<T> my_children;
        for (int i = 0; i < clusters[u].children.size(); i++) {
            if (marked_clusters[clusters[u].children[i]->index]) {
                my_children.push_back(index_map.find(clusters[u].children[i]->index));
            }
        }
        return my_children;

    });

    alt_parent_tree = parlay::tabulate(involved_nodes.size(),[&] (size_t u) {
        T node = involved_nodes[u]
        if (node==root->index) { //parent of root is itself
            return index_map.find(root->index);
        }
        else {
            return index_map.find(clusters[node].parent->index);
        }
    });
}
// T - index type
// D - type of data edges/clusters are storing
// Given clusters U and V, find the LCA in the original tree of the representative vertices of U and V. Store the cluster for which this is the representative vertex in ans. The LCA is oriented by the root of the RC tree
//root = root of RC tree
//r = root around which we orient the LCA answer
template<typename T, typename D>
void LCA(parlay::sequence<cluster<T,D>>& clusters,  cluster<T,D>* root, parlay::sequence<cluster<T,D>*> Uvertices, parlay::sequence<cluster<T,D>>* Vvertices, parlay::sequence<cluster<T,D>*> Rvertices, cluster<T,D>*& ans) {
    std::cout << "Fill in later" << std::endl;

    if (Uvertices.size() != Vvertices.size() || Uvertices.size() != Rvertices.size()) {
        std::cout << "Error abort, LCA query not well formed, cluster input arrays of different length " << Uvertices.size() << " " << Vvertices.size() << " " << Rvertices.size() << std::endl;
        exit(801);
    }

    int nc = clusters.size(); //vertices + edges~
    int k = Uvertices.size(); //batch size

    //mark which clusters are involved in this calculation
    //in clusters, are the first |V| positions the first |V| base vertex clusters? (or the inplace replacements of these) TOD2*
    parlay::sequence<bool> marked_clusters(nc,false);
    parlay::sequence<T> involved_nodes; //all nodes involved in this calculation
    find_nodes_involved(clusters,k,Uvertices,Vvertices,Rvertices,root->index,marked_clusters,involved_nodes);
 
    //we make a child and parent tree completely self contained (so we can call the static LCA method); we use the index_map to go back and forth between
    parlay::parlay_unordered_map<T,T>(2*involved_nodes.size()) index_map;
    parlay::sequence<parlay::sequence<T>> alt_child_tree;
    parlay::sequence<T> alt_parent_tree;
    create_map_trees(involved_nodes,clusters,root,alt_child_tree,alt_parent_tree,index_map);

    //create static LCA structure
    parlay::sequence<int> head(2*k+1);
    parlay::sequence<LCAnode<T>> augmented_vertices = parlay::tabulate(k,[&] (size_t i) {return i;});
    preprocess_par(alt_parent_tree,alt_child_tree,root->index,augmented_vertices,head);

    auto la_table = preprocess_la(alt_parent_tree,augmented_vertices,root->index); //level ancestors structure

    parlay::sequence<T> closest_boundary(alt_parent_tree.size(),-1);
    find_closest_boundary(involved_nodes,closest_boundary,clusters,root,alt_child_tree,alt_parent_tree,index_map);


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