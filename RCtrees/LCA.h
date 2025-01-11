//Code to handle LCA queries

#ifndef LCAHH //! cannot define function of this name, good to know
#define LCAHH

//header copied from RC.h
#include "RC.h"
#include "cluster.h"
//TOD2* using OLD parlay hash, new one not working
#include "include/old_parlay_hash/unordered_map.h" //for parlayhash concurrent hash map
#include "static_par_LCA.h"
#include "random_trees.h"

#include<bitset>

//anything to do with an alt tree/LCA helper
template<typename T, typename D>
struct LCAhelper {
    parlay::sequence<T> involved_nodes;
    parlay::sequence<std::bitset<50>> bitsets; //TOD2* change from 50? 
    parlay::sequence<parlay::sequence<T>> la_table;
    parlay::sequence<LCAnode<T>> augmented_vertices;
    parlay::sequence<T> common_boundaries;
    parlay::sequence<T> alt_parent_tree;
    parlay::sequence<parlay::sequence<T>> alt_child_tree;
    parlay::sequence<int> head;
    parlay::sequence<T> closest_boundary;
    size_t max_level;
    cluster<T,D>* root;
    int sn;
    parlay::sequence<bool> marked_clusters;

    //default constructor
    // LCAhelper() {
    //     //this->index_map = parlay::parlay_unordered_map<T,T>(1); //this will get overwritten
    // }
    
};

//get the active RC nodes in this batch, using a bottom up sweep. 
//store the nodes in play in involved_nodes, and store a boolean yes/no is_involved in marked_clusters
//in a batch of size k, will mark O(k log(1+n/k))
//bottom-up
//FIXME TOD2* go back and reimplement this with atomics, instead of a possibly too expensive remove_duplicates call 
//TOD2* what if queries are on the log n different levels of the RC tree? So they traverse up at different speeds, preventing a proper union?
//also: in star, could be O(N) duplicates
//!!!!!!!
//FIXME
//!!!!!!!
//TOD2* find the roots here, instead of passing a single root in to begin with (because is forest, there may be many roots
//TOD2* root should not be stored in h, rather pass in the root separately? Or maybe have an augmented pointer from a vertex to its root? 
template<typename T, typename D>
void find_nodes_involved_old(parlay::sequence<cluster<T,D>>& clusters, int k, parlay::sequence<std::tuple<T,T>>& queries, LCAhelper<T,D>& h) {

    parlay::sequence<T> stack;
    parlay::sequence<T> new_stack(2*k,-1);

    parlay::parallel_for(0,k,[&] (size_t i) {
        //pointer math to avoid contention
        new_stack[2*i] = std::get<0>(queries[i]);
        new_stack[2*i+1] = std::get<1>(queries[i]);
    });

    // std::cout <<"debugging" <<std::endl;
    // for (int i = 0; i < queries.size(); i++) {
    //     std::cout << std::get<0>(queries[i]) << " " << std::get<1>(queries[i]) << std::endl;
    // }
    // for (int i = 0; i < new_stack.size(); i++) {
    //     std::cout << new_stack[i] << " ";
    // }
    // std::cout << std::endl;

    //TOD2* is remove duplicates too expensive? I think it's O(k) for array of length k, so okay (though not nice on constants)
    //TOD2* does remove duplicates automatically output an array of reduced size? (or are there empty spaces? empty spaces would kill efficiency)
    stack=parlay::remove_duplicates(parlay::filter(new_stack,[&] (T item) {return item != -1;})); 
    h.involved_nodes.append(stack); //parallel bulk inserts the active nodes of this level 

    // for (int i = 0; i < stack.size(); i++) {
    //     std::cout << stack[i] << " ";
    // }
    // std::cout << std::endl;

    //bottom-up, propagate true marks
    while(stack.size() > 0) {
        new_stack=parlay::sequence<T>(stack.size(),-1);
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            //a is an index in clusters
            auto a = stack[i];
            h.marked_clusters[a]=true;
            if (a != h.root -> index) { //if we aren't at the root
                new_stack[i]=clusters[a].parent->index; //unlike last bottom up, doesn't get rid of duplicates because we do that explicitly
            }   

        });
        //don't pass up -1s (this is one way to filter, another way is to keep a filled_nodes bool and pack by it)
        stack=parlay::remove_duplicates(parlay::filter(new_stack,[&] (T item) {return item != -1;})); 
        h.involved_nodes.append(stack); //parallel insert the active nodes of this level 
    }

    //TOD2* is this last remove duplicates call (on the whole of involved_nodes) too expensive?
    h.involved_nodes = parlay::remove_duplicates(h.involved_nodes);

    h.sn = h.involved_nodes.size(); //sn holds the # of nodes in the subset tree (the # of active internal nodes of the RC tree for this batch)
}

template<typename T, typename D>
void find_nodes_involved(parlay::sequence<cluster<T,D>>& clusters, int k, parlay::sequence<std::tuple<T,T>>& queries, LCAhelper<T,D>& h) {

    std::cout << "alive" << std::endl;

    parlay::sequence<T> stack;
    parlay::sequence<T> new_stack(2*k,-1);

    parlay::parallel_for(0,k,[&] (size_t i) {
        //pointer math to avoid contention
        new_stack[2*i] = std::get<0>(queries[i]);
        new_stack[2*i+1] = std::get<1>(queries[i]);
    });

    //this call of remove duplicates okay because is linear work in k, log depth (single call okay, repeated calls could be trouble)
    stack=parlay::remove_duplicates(parlay::filter(new_stack,[&] (T item) {return item != -1;})); 
    h.involved_nodes.append(stack); //parallel bulk inserts the active nodes of this level 

    //check whether node is being accessed first or not
    auto fetch_counters = parlay::tabulate(clusters.size(),[&] (size_t i) {return std::atomic<int>(1); });
    //-4 = unprocessed; -3 = in process; -2 = done, don't touch
    auto fetch_bool = parlay::tabulate(clusters.size(),[&] (size_t i) {return -4; });

    //mark that the queries have already been added to stack (ex. if an internal node and one of its descendants are both offered as queries, don't want to go up section from internal node to root twice)
    parlay::parallel_for(0,stack.size(),[&] (size_t i) {
        fetch_bool[stack[i]]=-2; 
    });

    //bottom-up, propagate true marks
    while(stack.size() > 0) {
            std::cout << "alive" << std::endl;

        new_stack=parlay::sequence<T>(stack.size(),-1);
        //set up race on adding parents
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            //a is an index in clusters
            auto a = stack[i];
            h.marked_clusters[a]=true;
            if (a != h.root -> index && fetch_bool[clusters[a].parent->index]==-4) { //if we aren't at the root and this is the first iteration the parent is accessed, then set up a race on the parent. 
                fetch_bool[clusters[a].parent->index]=-3;
            }   

        });

        //run the race on who gets to parent first, ensures parent only added once to stack
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            //a is an index in clusters
            auto a = stack[i];
            
            if (a != h.root -> index) {
                //must confirm not at root before accessing parent
                auto p = clusters[a].parent->index;
                    //only fetch add if fetch bool is true (relies on short circuiting of ands)
                    if (fetch_bool[p]==-3 && fetch_counters[p].fetch_add(2) == 1) { //if we aren't at the root and this is the first time this parent is accessed, then add the parent to the stack
                    new_stack[i]=p;
                    fetch_bool[p]=-2; //signal we are done with this node
                }
            }   

        });
        
        //don't pass up -1s (this is one way to filter, another way is to keep a filled_nodes bool and pack by it)
        stack=parlay::remove_duplicates(parlay::filter(new_stack,[&] (T item) {return item != -1;})); 
        h.involved_nodes.append(stack); //parallel insert the active nodes of this level 
    }

   
    h.sn = h.involved_nodes.size(); //sn holds the # of nodes in the subset tree (the # of active internal nodes of the RC tree for this batch)
        std::cout << "out" << std::endl;

}

//given a cluster whose index in involved_nodes is i_index, find the index (in alt_tree) of its boundary that is closer to the root, return it
//requires closest boundary to be computed for the parent of s (hence the top down computation with a stack)
template<typename T, typename D>
T instack_choose_boundary(T s, parlay::sequence<cluster<T,D>>& clusters, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {

    cluster<T,D>* cur = &clusters[h.involved_nodes[s]];
    cluster<T,D>* parent = &clusters[h.involved_nodes[h.alt_parent_tree[s]]];

  

    std::cout << "instack choose boundary for vertex " << cur->index << std::endl;
    // std::cout << "parent pointer " << parent << std::endl;
    std::cout << "\t parent index of " << cur->index << " is " << parent->index << std::endl;

    //if parent is root, then the root is the desired boudary
    if (parent->index == h.root->index) {
        std::cout << "\t choosing root" << std::endl;
        return index_map[h.root->index];

    }

    cluster<T,D>* left;
    cluster<T,D>* right;
    D def_val;
    // std::cout << left << " " << right << std::endl;

    parent->find_boundary_vertices(left,def_val,right,def_val,def_val);
    // std::cout << left << " " << right << std::endl;
    std::cout << "\t parent left boundary: " << left->index << std::endl;
    std::cout << "\t parent right boundary: " << right->index << std::endl;

    cluster<T,D>* child_left;
    cluster<T,D>* child_right;
    cur->find_boundary_vertices(child_left,def_val,child_right,def_val,def_val);

    

    std::cout << "\t child left boundary: " << child_left->index << std::endl;
    std::cout << "\t child right boundary: " << child_right->index << std::endl;

    std::cout << "\t parent state: " << parent->first_contracted_node->next->state << std::endl;

    //if the child (current vertex) is unary, then only one boundary to possibly return
    if (child_left == child_right) {
        return *index_map.Find(child_left->index);
    }

    //if the parent of this cluster is a binary cluster
    //TOD2* doublecheck or test? very possible is typo here
    if (parent->first_contracted_node->next->state & (binary_cluster | base_edge) ) {
        T par_closest = h.closest_boundary[*index_map.Find(parent->index)];
        //if the boundary of the parent that is closest to the root is one of the child's boundaries
        if (par_closest == *index_map.Find(child_left->index) || par_closest==*index_map.Find(child_right->index)) {
            std::cout << "\t choosing par closest" << std::endl;
            return par_closest;
        }
        else { //pick the boundary NOT shared by the parent
            if (left->index==child_left->index || right->index == child_left->index) {
                std::cout << "\t choosing right child boundary" << std::endl;

                return *index_map.Find(child_right->index);
            } 
            else if (left->index == child_right->index || right->index == child_right->index) {
                std::cout << "\t choosing left child boundary" << std::endl;

                return *index_map.Find(child_left->index);
            }
            else {
                std::cout << "error, parent shares no boundary with child, abort " << std::endl;
                exit(2001);
            }
        }

    }
    //TOD2* is nullary cluster one of the involved nodes?
    else if (parent->first_contracted_node->next->state & (nullary_cluster)) {
        std::cout << "did not account for this abort" << std::endl;
        exit(802);
    }
    else { //PARENT is unary cluster
        std::cout << "\t parent is unary "<< std::endl;

        if (child_left->index != child_right->index) { //then kid is binary cluster, in which case we want to return the parent's boundary 
            return *index_map.Find(left->index);
        }
        else { //kid is unary cluster, should have been caught earlier, abort
            std::cout << "Abort, kid is unary but this function should have returned on the earlier unary check" << std::endl;
            exit(803);

        }
    }

}


//in involved_nodes,
//for each binary cluster, find which boundary vertex is closer to the root
//for each unary cluster will just return the only boundary vertex
//TOD2* does the find_boundary_vertices call within this function look outside the subset of the tree we are considering, breaking the work bound? Think okay because constant max # of children.
//top-down computation
//if root, closest boundary is itself by convention
template<typename T, typename D>
void find_closest_boundary(parlay::sequence<cluster<T,D>>& clusters, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {

    //std::cout << "h3.11" << std::endl;

    h.closest_boundary = parlay::sequence<T>(h.sn,-1);
    parlay::sequence<T> stack;
    // std::cout << "root index is " << h.root -> index << std::endl;
    // std::cout << "alt index is " << *index_map.Find(h.root->index) << std::endl;
    stack.push_back(*index_map.Find(h.root->index));

    //std::cout << "h3.12" << std::endl;

    parlay::sequence<parlay::sequence<T>> new_stack;
    int count = 0;
    while (stack.size() > 0) {
        // std::cout << "printing stack iter " << count << std::endl;
        // for (int i = 0; i < stack.size(); i++) {
        //     std::cout << stack[i] << " ";
        // }
        // std::cout << std::endl;
        new_stack = parlay::sequence<parlay::sequence<T>>(stack.size(),parlay::sequence<T>());
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            //s is an index of the alt tree (not of RC tree/clusters)
            auto s = stack[i];
            //don't run a boundary check on the root, because root doesn't have boundaries
            if (h.involved_nodes[s]==h.root->index) {
                h.closest_boundary[s]=s;
            }
            else {
                h.closest_boundary[s] = instack_choose_boundary(s,clusters,h,index_map);
            }
            for (T j : h.alt_child_tree[s]) {
                new_stack[i].push_back(j);
            }     
        });
        stack=parlay::flatten(new_stack);
        if (count > 10) exit(1202);
        count += 1;
    }
}

//given an RC tree and a subset of the RC tree we are working on, extract this subset and put it into its own tree structure. 
//alt_child_tree holds this subset tree from the child's perspective: each node knows its list of children
//alt_parent_tree holds this (SAME) subset tree from the parent's perspective: each node only knows its parent
//index_map bridges the alt trees with the RC tree clusters: index_map[cluster index] = alt index. To go from alt tree to clusters, use involved_nodes: involved_nodes[alt index] = cluster index.
template<typename T, typename D>
void create_map_trees(parlay::sequence<cluster<T,D>>& clusters, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {
    //we make a child and parent tree completely self contained (so we can call the static LCA method); we use the index_map to go back and forth between

    //map the value in involved_nodes to its index; we will then create the tree using the index
    //TOD2* is there a bulk insert command?
    parlay::parallel_for(0,h.sn,[&] (T i) {
        index_map.Insert(h.involved_nodes[i],i);
    });

    //invariant: index_map and involved_nodes are inverses of each other
    bool should_error = false;
    parlay::parallel_for(0,h.sn,[&] (T i) {
        if (*index_map.Find(h.involved_nodes[i]) != i) {
            should_error=true;
            
        }
    });
    if (should_error) {
        std::cout << "aborting, index map broken" << std::endl;
        exit(1201);
    }

  
    h.alt_child_tree = parlay::tabulate(h.sn,[&] (size_t s) {
        T a = h.involved_nodes[s]; //a is the index in the RC tree/clusters list
        parlay::sequence<T> my_children;
        for (int j = 0; j < clusters[a].children.size(); j++) {
            //TOD2* is checking all children for membership like this too expensive?
            if (clusters[a].children[j] != nullptr && h.marked_clusters[clusters[a].children[j]->index]) {
                my_children.push_back(*index_map.Find(clusters[a].children[j]->index));
            }
        }
        return my_children;

    });

    h.alt_parent_tree = parlay::tabulate(h.sn,[&] (size_t u) {
        T a = h.involved_nodes[u];
        if (a==h.root->index) { 
            return *index_map.Find(h.root->index); //parent of root is itself
        }
        else {
            return *index_map.Find(clusters[a].parent->index);
        }
    });

   
}


//preprocessing for LCA. This section of the preprocessing creates the static structures (|| static LCA, level ancestors) that the algorithm uses as lemmas
template<typename T, typename D>
void static_preprocess(parlay::sequence<cluster<T,D>>& clusters, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {


    T alt_root = *index_map.Find(h.root->index);
    if (h.alt_parent_tree[alt_root] != alt_root) {
        std::cout << "error, root's parent is not itself " << std::endl;
        exit(2201);
    }
    //create static LCA structure
    h.head = parlay::sequence<int>(2*h.sn+1);
    h.augmented_vertices = parlay::sequence<LCAnode<T>>(h.sn);

    //give augmented vertices a (default) id. This id just used for printing.
    parlay::parallel_for(0,h.augmented_vertices.size(),[&] (size_t i) {
        h.augmented_vertices[i].id=i;

    });

    preprocess_par(h.alt_parent_tree,h.alt_child_tree,alt_root,h.augmented_vertices,h.head);

    // std::cout << "printing augment" << std::endl;
    // std::cout << "id,\t inl,\t preo,\t#ri0,\tlvl,\tascd,\t size" << std::endl;
    // for (int i =0 ; i < h.augmented_vertices.size(); i++) {
    //     h.augmented_vertices[i].print(); 
    // }
    //std::cout << "starting queries" << std::endl;


    //get the max level in the tree
    //TOD2* compute this once, outside (for constant factor efficiency)
    h.max_level = *parlay::max_element(parlay::map(h.augmented_vertices,[&] (LCAnode<T> node) { return node.level; }));

    h.la_table = preprocess_la(h.alt_parent_tree,h.augmented_vertices,alt_root); //level ancestors structure

}


//for each vertex write down if each of its ancestors is binary or unary. Do with a single top down sweep. Because we are only writing down a constant amount of information each time, this is O(# nodes) = O(k log(1+n/k)) work, instead of the O(k \log n) that naive level ancestors needs
template<typename T, typename D>
void get_ancestor_bitset(parlay::sequence<cluster<T,D>>& clusters, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {

    parlay::sequence<T> stack;

    int r = *index_map.Find(h.root->index);
    stack.push_back(r);
    //nullary is 0
    for (int iter = 0; iter < h.max_level+1; iter++) {
        parlay::sequence<parlay::sequence<T>> new_stack=parlay::sequence(stack.size(),parlay::sequence<T>());

        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            auto s = stack[i];

            h.bitsets[s]=h.bitsets[h.alt_parent_tree[s]]; //copy bitset from parent

            //if this cluster is unary
            if (clusters[h.involved_nodes[s]].first_contracted_node->next->state & (unary_cluster)) {
                h.bitsets[s].set(iter,1); //store a 1
            }
            for (T j : h.alt_child_tree[s]) {
                new_stack[i].push_back(j);
            }

        });

        stack=parlay::flatten(new_stack);

    }

}


//do additional preprocessing in preparation for queries. These structures give additional information about the RC tree (ancestors, boundaries, relation to root)
template<typename T, typename D>
void RC_preprocess(parlay::sequence<cluster<T,D>>& clusters, parlay::sequence<std::tuple<T,T>>& queries, int k, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {

    //std::cout << "h3.1" << std::endl;
    find_closest_boundary(clusters,h,index_map);
   //std::cout << "h3.2" << std::endl;

   
    //store ancestor information (binary or unary) in a bitset
    //in the i^th index in the bitset, the ancestor (of v) in level i is held
    //so in the 0^th position, the binary/unary value of the root is held
    //TOD2* can a variable sized (at compile) bitset be allocated like this? Or set to 50 to be safe? But very small probability that tree is too deep? (not with the deterministic contraction scheme?)
    //TOD2* assume all initialized as false (binary)?
    //TOD2* what to write for nullary cluster? currently writing 0
    h.bitsets = parlay::sequence<std::bitset<50>>(h.sn,0); 
    get_ancestor_bitset(clusters,h,index_map);

    std::cout << "printing all ancestor bitsets (cl. in.)" << std::endl;
    for (int i = 0 ; i < h.bitsets.size(); i++) {
        std::cout << "\t" << h.involved_nodes[i] << " " << h.bitsets[i].to_ulong() << std::endl;

    }

   // std::cout << "h3.3" << std::endl;

    //get common boundary with static LCA structure. Gives ALT TREE index
    h.common_boundaries = parlay::tabulate(k,[&] (size_t i) {
        T s = *index_map.Find(std::get<0>(queries[i])); //ALT TREE INDICES
        T t = *index_map.Find(std::get<1>(queries[i])); 
        return query(h.head,h.alt_parent_tree,h.augmented_vertices,s,t);

    });

    print_parent_tree(h.alt_parent_tree,"parent tree reprint");
    print_child_tree(h.alt_child_tree,"child tree reprint");

    std::cout << "common boundary debug" << std::endl;
    for (int i = 0; i < h.common_boundaries.size(); i++) {
        std::cout << "common boundary of " << *index_map.Find(std::get<0>(queries[i])) << " and " << *index_map.Find(std::get<1>(queries[i])) << " is " << h.common_boundaries[i] << std::endl;
    }

    //minor common boundary tests:
    bool throw_error = false;
    parlay::parallel_for(0,h.common_boundaries.size(),[&] (size_t i) {
        T s = *index_map.Find(std::get<0>(queries[i])); 
        if (query(h.head,h.alt_parent_tree,h.augmented_vertices,s,h.common_boundaries[i]) != h.common_boundaries[i])  {
            throw_error=true;
        }
        s = *index_map.Find(std::get<1>(queries[i])); 
        if (query(h.head,h.alt_parent_tree,h.augmented_vertices,s,h.common_boundaries[i]) != h.common_boundaries[i])  {
            throw_error=true;
        }

    });
    if (throw_error) {
        std::cout << "error abort, common boundary not ancestor of constituent" << std::endl;
        exit(2101);
        
    }

   // std::cout << "h3.4" << std::endl;


}


//given cluster list (clusters), a cluster (B), and a vertex (corresponds to index in clusters) u, find the vertex in the cluster path of B that is closest to u
//bitset is for u, helps us find the highest unary cluster
//TOD2* can I special set the bitset to not a fixed constant? Need to mess with at some point
//uses level ancestors as a subroutine
//takes in index u on CLUSTERS
//returns as index on CLUSTERS not on alt_tree
template<typename T, typename D>
T closest_on_cluster_path(parlay::sequence<cluster<T,D>>& clusters, cluster<T,D>* B, T u, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {
    std::cout << "running closest on cluster path for " << u  << " in " << B->index << std::endl;

    std::bitset<50> ubits=h.bitsets[*index_map.Find(u)]; //TOD2* correct fixed size for larger tree depths

    std::cout << "\tubits is " << ubits.to_ulong() << std::endl;

    int u_level = h.augmented_vertices[*index_map.Find(u)].level; //U's level in RC tree
    int b_level = h.augmented_vertices[*index_map.Find(B->index)].level;
    if (u_level < b_level) {
        std::cout << "error, B is ancestor of U in cluster path check" << std::endl;
        exit(1001);
    }

    std::cout << "\t" << u << " level: " << u_level << ", " << B->index << " level: " << b_level << std::endl;

    //zero out all but the first u_level+1 bits
    //TOD2* is casting to *unsigned* long an issue? vs signed long? (bitset only offers unsigned cast, I'd need to do an additional cast on top of it)
    std::bitset<50> rel_bits = ubits.to_ulong() & ((1 << (u_level+1))-1);

    //zero out any bits before bit b_level
    rel_bits = (rel_bits.to_ulong() / (1 << (b_level) )) * (1 << (b_level) );

    std::cout << "\trel bits is " << rel_bits << std::endl;
    
    if (rel_bits.to_ulong() == 0) { //if there are no unary clusters in this range, u itself must be on the cluster path TOD2* b itself must be the best node on cluster path?
        std::cout << "\tno unary clusters on path, return u itself" << std::endl;
        return u;
    }
    else { //there are at least 1 unary cluster, in which case use level ancestors on the highest one
        std::cout << "\tat least one unary cluster on path " << std::endl;
        int bit_level = r1(rel_bits.to_ulong()); //leftmost (highest) unary cluster is rightmost (because level 0 is root)
        std::cout << "\t bit level is " << bit_level << std::endl;
        int ancestor_comp_level = u_level - bit_level; //the i^th ancestor of u
        std::cout << "\t level in ancestor structure we are checking is " << ancestor_comp_level << std::endl;
        T highest_unary = query_la(h.la_table,*index_map.Find(u),ancestor_comp_level);
        std::cout << "\tthat unary is " << highest_unary << std::endl;
        cluster<T,D>* left;
        D def_val;
        //the vertex on the cluster path is the (only) boundary vertex of the highest unary
        clusters[h.involved_nodes[highest_unary]].find_boundary_vertices(left,def_val,left,def_val,def_val);

        std::cout << "\tboundary of that unary is " << left->index << std::endl;
        return left->index;

    }   
}


//given a vertex b and a vertex f inside cluster B, find the cluster that is an immediate child of B (in the RC tree) that has f as an ancestor
//return as a CLUSTER index
template<typename T, typename D>
T highest_ancestor_in_cluster(T b, T f, parlay::sequence<cluster<T,D>>& clusters, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {
    if (b==f) {
        return f; //if common boundary is the cluster, then return itself?
    }
    //std::cout << "high ancestor " << "b: " << b << ", f: " << f << std::endl;
    for (int i = 0; i < clusters[b].children.size(); i++) {
        cluster<T,D>* child = clusters[b].children[i];
        //don't look at nullptr, nor clusters that are not part of the relevant subset
        //TOD2* causes work efficient issues? Think okay
        if (child==nullptr || h.marked_clusters[child->index]==false) continue; 

        //std::cout << "looking at child: " << child->index << std::endl;
        T chi_alt = *index_map.Find(child->index);
        //if this child is the ancestor of f
        if (query(h.head,h.alt_parent_tree,h.augmented_vertices,*index_map.Find(f),chi_alt) == chi_alt) {
            return child->index;
        }
    }
    //should never reach here because if f is a child of b, then one of B's children must contain f
    std::cout << "should never have reached here in highest ancestor, error" << std::endl;
    std::cout << "b is: " << b << ", " << "f is " << f << std::endl;
    exit(1003);
}



//given cluster represented by vertex a, vertex b inside cluster A (could be very deep inside or shallow), check if b is between f and the root
//requirement: f is a child of b in the RC tree
//returns true if b is between f and the root
template<typename T, typename D>
bool is_cluster_between(T b, T f, parlay::sequence<cluster<T,D>>& clusters, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {

    if (h.augmented_vertices[index_map[b]].level > h.augmented_vertices[index_map[f]].level) {
        std::cout << "error, f higher in RC tree than b " << std::endl;
        exit(1002);
    }
    //get the highest cluster containing f; note that X has b has a boundary vertex
    T x = highest_ancestor_in_cluster(b,f,clusters,h,index_map);
    T x_boundary_rd = h.involved_nodes[h.closest_boundary[*index_map.Find(x)]]; //boundary of x in direction of root

    //get the boundaries of B
    cluster<T,D>* left;
    D def_val;
    cluster<T,D>* right;
    clusters[b].find_boundary_vertices(left,def_val,right,def_val,def_val);

    std::cout << "Asking the question is cluster " << b << " between " << f << " and root?" << std::endl;
    std::cout << "\thighest ancestor of " << f << "descended from  " << b << " is " << x << std::endl;
    std::cout << "\tthe boundary of " << x << " in direction of root is " << x_boundary_rd << std::endl;
    std::cout << "\tthe boundary of " << b << " are " << left->index << " and " << right->index << std::endl;

    std::cout << "\treturning " << !(x_boundary_rd == left->index || x_boundary_rd == right->index) << std::endl;

    //if x's boundary facing the root is one of B's boundaries, then x is between so false
    if (x_boundary_rd == left->index || x_boundary_rd == right->index) return false;
    //if x's boundary that faces r does not face outward, then it faces b, so return true
    return true;


}



//do subcase of casework for LCA calculation, when c==u or c==v
template<typename T, typename D>
T case_c_is_uv(parlay::sequence<cluster<T,D>>& clusters,parlay::sequence<std::tuple<T,T>>& queries, parlay::sequence<T>& answers, int k, LCAhelper<T,D>& h, T c, T u, T v,parlay::parlay_unordered_map<T,T>& index_map) {

    std::cout << "in root case/c is uv case" << std::endl;

    T r = h.root->index; 

    if (u==r) { //if u root, u is LCA
        return u;
    }
    else if (v == r) { //if v root, v is LCA
        return v;
    }
    //if the common boundary is original root of the RC tree, common boundary is LCA
    else if (c==r) {
        return c;
    }

    //here, the chosen vertex is the vertex that is NOT equal to c
    T other_vertex = u; //which of u or v c is equal to
    if (c==u) other_vertex=v;
    //is x well defined, for c==other_vertex? (piont of this call?)
    T x = highest_ancestor_in_cluster(c,other_vertex,clusters,h,index_map);

    if (clusters[x].first_contracted_node->next->state & unary_cluster) { //if X is unary
        return c;
    }
    else { //if X is binary
        if (is_cluster_between(c,other_vertex,clusters,h,index_map)) {
            return c; //is c not other_vertex (because c in between, c is LCA)
        }
        else {
            return closest_on_cluster_path(clusters,&clusters[x],other_vertex,h,index_map);

        }

    }

}

//do subcase of LCA query where the common boundary is not u nor v, and none of c,u,v are the root
template<typename T, typename D>
T case_c_not_uv(parlay::sequence<cluster<T,D>>& clusters,parlay::sequence<std::tuple<T,T>>& queries, parlay::sequence<T>& answers, int k, LCAhelper<T,D>& h, T c, T u, T v,parlay::parlay_unordered_map<T,T>& index_map) {
  //  std::cout << "h4.11" << std::endl;
    T x = highest_ancestor_in_cluster(c,u,clusters,h,index_map);
   //     std::cout << "h4.12" << std::endl;

    T y = highest_ancestor_in_cluster(c,v,clusters,h,index_map);

    //       std::cout << "h4.125" << std::endl;

    //TOD2* could plug in x/y here, should not change result
    bool u_bet = is_cluster_between(c,u,clusters,h,index_map);
    //    std::cout << "h4.13" << std::endl;

    bool v_bet = is_cluster_between(c,v,clusters,h,index_map);
    //    std::cout << "h4.14" << std::endl;

    std::cout << "Debugging info for case c not uv " << std::endl;
    std::cout << "\tu's highest ancestor in c is " << x << std::endl;
    std::cout << "\tv's high ancestor in c is " << y << std::endl;
    std::cout << "\tc is between u and root is " << u_bet << std::endl;
    std::cout << "\tc is between v and root is " << v_bet << std::endl;

    if (u_bet && v_bet) return c;
    else if (u_bet) return closest_on_cluster_path(clusters,&clusters[y],v,h,index_map);
    else if (v_bet) return closest_on_cluster_path(clusters,&clusters[x],u,h,index_map);

    std::cout << "should never reach here case c not uv" << std::endl;
    exit(1003);

}


//after preprocessing, answer each query in parallel. Separate into cases based on c and r
template<typename T,typename D>
void batch_fixed_LCA_casework(parlay::sequence<cluster<T,D>>& clusters,parlay::sequence<std::tuple<T,T>>& queries, parlay::sequence<T>& answers, int k, LCAhelper<T,D>& h,parlay::parlay_unordered_map<T,T>& index_map) {
    T r = h.root->index;
    parlay::parallel_for(0,k,[&] (size_t i) {
        T u = std::get<0>(queries[i]); //queries, in CLUSTER index
        T v = std::get<1>(queries[i]); //CLUSTER index
        T c = h.involved_nodes[h.common_boundaries[i]]; //common boundary CLUSTER index

        printf("(u,v,c): (%d %d %d)\n",u,v,c);
        //first cases, relations between u,v,c,r
        if (u==r || v == r || c ==r || c == u || c == v) //if u root, u is LCA
            answers[i]=case_c_is_uv(clusters,queries,answers,k,h,c,u,v,index_map);
        //if the common boundary is not u nor v, and the root is not u,v,c
        else 
            answers[i]=case_c_not_uv(clusters,queries,answers,k,h,c,u,v,index_map);
    });

}

// T - index type
// D - type of data edges/clusters are storing
//Given a sequence of queries in the form of (u_i,v_i), where u_i,v_i are vertices (not clusters), find the LCA of u_i and v_i in the tree rooted at **the RC tree root** and store the answer in the sequence answers. 
//root = original root of RC tree
template<typename T, typename D>
void batch_fixed_LCA(parlay::sequence<cluster<T,D>>& clusters,  cluster<T,D>* root, parlay::sequence<std::tuple<T,T>>& queries, parlay::sequence<T>& answers) {

    //std::cout << "h1" << std::endl;

    //std::cout << "new call of batch fixed LCA" << std::endl;

    int nc = clusters.size(); //vertices + edges~
    int k = queries.size(); //batch size

    LCAhelper<T,D> h;//= LCAhelper<T,D>();
    h.root=root;

    //mark which clusters are involved in this calculation
    //in clusters, are the first |V| positions the first |V| base vertex clusters? (or the inplace replacements of these) TOD2*
    //marked_clusters follows the indexing in clusters
    h.marked_clusters = parlay::sequence<bool>(nc,false);
    //find all nodes involved in this calculation. the indexing in involved_nodes is arbitrary (but will become important, beacuse we will map to it with index_map)
    find_nodes_involved(clusters,k,queries,h);

    parlay::parlay_unordered_map<T,T> index_map = parlay::parlay_unordered_map<T,T>(2*h.sn); //will get overwritten

    std::cout << "printing involved nodes in form (alt index,original index)" <<std::endl;
    for (int i = 0; i < h.involved_nodes.size(); i++) {
        std::cout << "\t" << i << " " << h.involved_nodes[i] << std::endl;
    }

    //create the subset of the tree we are working on (from involved_nodes), the hashmap to get from the RCtree to this alt tree and back
    create_map_trees(clusters,h,index_map);

    //std::cout << "h2" << std::endl;
    //prepare the static data structure (unrelated to RC)
    static_preprocess(clusters,h,index_map);

    //std::cout << "h3" << std::endl;


    //print_child_tree(h.alt_child_tree,"printing RC subset child tree");
    print_parent_tree(h.alt_parent_tree,"printing RC subset parent tree");
    // for (int i = 0; i < h.involved_nodes.size(); i++) {
    //     std::cout << i << " " << h.involved_nodes[i] << " " << *index_map.Find(h.involved_nodes[i]) << std::endl;
    // }



    //prepare the RC-specific augmented info
    RC_preprocess(clusters,queries,k,h,index_map);

    std::cout << "printing out closest boundary (in cluster index) : " << std::endl;
    for (int i = 0; i < h.closest_boundary.size(); i++) {
        std::cout << "(" << h.involved_nodes[i] << ", " << h.involved_nodes[h.closest_boundary[i]] <<") ";
    }
    std::cout << std::endl;

    //std::cout << "h4" << std::endl;



    batch_fixed_LCA_casework(clusters,queries,answers,k,h,index_map);

    //std::cout << "h5" << std::endl << std::endl << std::endl;


}


//given queries of form (u_i,v_i,r_i), find the LCA when the tree is rooted at r_i
//first, answer queries based on the original root of the RC tree, then combine with logic_lca
template<typename T, typename D>
void batchLCA(parlay::sequence<cluster<T,D>>& clusters,  cluster<T,D>* root, parlay::sequence<std::tuple<T,T,T>>& queries, parlay::sequence<cluster<T,D>*>& answers) {


    int k = queries.size();

    //do LCA queries for (u_i,v_i) (u_i,r_i) (v_i,r_i)
    parlay::sequence<T> answers_uv(k,-1);
    parlay::sequence<T> answers_ur(k,-1);
    parlay::sequence<T> answers_vr(k,-1);
    parlay::sequence<std::tuple<T,T>> queries_uv=parlay::tabulate(k,[&] (size_t i) {
        return std::make_tuple(std::get<0>(queries[i]),std::get<1>(queries[i]));
    });
    parlay::sequence<std::tuple<T,T>> queries_ur=parlay::tabulate(k,[&] (size_t i) {
        return std::make_tuple(std::get<0>(queries[i]),std::get<2>(queries[i]));
    });
    parlay::sequence<std::tuple<T,T>> queries_vr=parlay::tabulate(k,[&] (size_t i) {
        return std::make_tuple(std::get<1>(queries[i]),std::get<2>(queries[i]));
    });

    std::cout << "query print" << std::endl;
    for (int i = 0; i < answers_uv.size(); i++) {
        std::cout << i << ": " << std::get<0>(queries[i]) << " " << std::get<1>(queries[i]) << " " << std::get<2>(queries[i]) << " " << answers_uv[i] << " " << answers_ur[i] << " " << answers_vr[i] << std::endl;
    }
    std::cout << "Starting 1st query call u v" << std::endl;
    batch_fixed_LCA(clusters,root,queries_uv,answers_uv);
    std::cout << std::endl << std::endl << std::endl;

    std::cout << "Starting 2nd query call u r" << std::endl;
    batch_fixed_LCA(clusters,root,queries_ur,answers_ur);
    std::cout << std::endl << std::endl << std::endl;


    std::cout << "Starting 3rd query call v r" << std::endl;
    batch_fixed_LCA(clusters,root,queries_vr,answers_vr);
    std::cout << std::endl << std::endl << std::endl;

    std::cout << "answer print" << std::endl;
    for (int i = 0; i < answers_uv.size(); i++) {
        std::cout << i << ": " << std::get<0>(queries[i]) << " " << std::get<1>(queries[i]) << " " << std::get<2>(queries[i]) << " " << answers_uv[i] << " " << answers_ur[i] << " " << answers_vr[i] << std::endl;
    }

    answers=parlay::tabulate(k,[&] (size_t i) {
        return &clusters[logic_lca(answers_uv[i],answers_ur[i],answers_vr[i])]; //this logic in Vanilla LCA
    });

}


#endif //LCA