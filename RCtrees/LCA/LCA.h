//Code to handle LCA queries
//This function takes in input directly, converts to single tree, fixed root version, then calls the LCA code in fixedLCA.h

#ifndef LCAHH //! cannot define function of this name, good to know
#define LCAHH

//header copied from RC.h
#include "RC.h"
#include "cluster.h"
//TOD2* using OLD parlay hash, new one not working
#include "include/old_parlay_hash/unordered_map.h" //for parlayhash concurrent hash map
#include "static_par_LCA.h"
#include "static_seqn_LCA.h"
#include "random_trees.h"
#include "fixedLCA.h"

#include<bitset>


//given a set of roots and involved nodes, assign each involved node to a root
template<typename T, typename D>
void populate_root_ids(parlay::sequence<cluster<T,D>>& clusters, parlay::sequence<bool>& marked_clusters,parlay::sequence<T>& involved_nodes, parlay::sequence<T>& roots,parlay::sequence<T>& root_ids,parlay::parlay_unordered_map<T,T>& index_map) {

    //can run roots independently because disconnected
    parlay::parallel_for(0,roots.size(),[&] (size_t i) {
        parlay::sequence<T> stack;
        parlay::sequence<parlay::sequence<T>> new_stack;
        stack.push_back(roots[i]);
        while (stack.size() > 0) {
            new_stack = parlay::sequence<parlay::sequence<T>>(stack.size(),parlay::sequence<T>());
            parlay::parallel_for(0,stack.size(),[&] (size_t j) {
                T a = stack[j]; //current tree node index we are looking at
                root_ids[*index_map.Find(a)]=i;
                //TOD2* is checking all children within bound (given many of the children could not be marked)? Requires bounded # of children here
                //Add children to stack, to propagate down tree
                for (int l = 0; l < clusters[a].children.size(); l++) {
                    if (clusters[a].children[l] != nullptr && marked_clusters[clusters[a].children[l]->index]) {
                        new_stack[j].push_back(*index_map.Find(clusters[a].children[l]->index));
                    }
                }
            });
            stack = parlay::flatten(new_stack);
        }
    });
}



//forest version of getting involved nodes, roots
template<typename T, typename D>
void forest_find_nodes_involved(parlay::sequence<cluster<T,D>>& clusters, int k, parlay::sequence<std::tuple<T,T,T>>& queries, parlay::sequence<bool>& marked_clusters,parlay::sequence<T>& involved_nodes, parlay::sequence<T>& roots) {

    parlay::sequence<T> stack;
    parlay::sequence<T> new_stack(3*k,-1);

    parlay::parallel_for(0,k,[&] (size_t i) {
        new_stack[3*i] = std::get<0>(queries[i]);
        new_stack[3*i+1] = std::get<1>(queries[i]);
        new_stack[3*i+2] = std::get<2>(queries[i]);

    });

    stack=parlay::remove_duplicates(parlay::filter(new_stack,[&] (T item) {return item != -1;})); 
    involved_nodes.append(stack);
    roots.append(parlay::filter(stack,[&] (T item) {return isNullary(&clusters[item]); })); //get initial roots


    auto fetch_counters = parlay::tabulate(clusters.size(),[&] (size_t i) {return std::atomic<int>(1); });
    //-4 = unprocessed; -3 = in process; -2 = done, don't touch
    auto fetch_bool = parlay::tabulate(clusters.size(),[&] (size_t i) {return -4; });

    parlay::parallel_for(0,stack.size(),[&] (size_t i) {
        fetch_bool[stack[i]]=-2; 
    });

    //bottom-up, propagate true marks
    while(stack.size() > 0) {

        new_stack=parlay::sequence<T>(stack.size(),-1);

        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            //a is an index in clusters
            auto a = stack[i];
            h.marked_clusters[a]=true;
            if (!isNullary(&clusters[a]) && fetch_bool[clusters[a].parent->index]==-4) { //if we aren't at the root and this is the first iteration the parent is accessed, then set up a race on the parent. 
                fetch_bool[clusters[a].parent->index]=-3;
            }   

        });

        //run the race on who gets to parent first, ensures parent only added once to stack
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            //a is an index in clusters
            auto a = stack[i];
            
            if (!isNullary(&clusters[a])) { 
                //must confirm not at root before accessing parent
                auto p = clusters[a].parent->index;
                    if (fetch_bool[p]==-3 && fetch_counters[p].fetch_add(2) == 1) { //if we aren't at the root and this is the first time this parent is accessed, then add the parent to the stack
                    new_stack[i]=p;
                    fetch_bool[p]=-2; //signal we are done with this node
                }
            }   

        });
        
        //don't pass up -1s (this is one way to filter, another way is to keep a filled_nodes bool and pack by it)
        stack=parlay::remove_duplicates(parlay::filter(new_stack,[&] (T item) {return item != -1;})); 
        involved_nodes.append(stack); //parallel insert the active nodes of this level 
        roots.append(parlay::filter(stack,[&] (T item) {return isNullary(&clusters[item]); })); //add any roots found

    }

}


//given queries of form (u_i,v_i,r_i), find the LCA when the tree is rooted at r_i
//first, answer queries based on the original root of the RC tree, then combine with logic_lca
template<typename T, typename D>
void batchLCA_singletree(parlay::sequence<cluster<T,D>>& clusters,  cluster<T,D>* root, parlay::sequence<std::tuple<T,T,T>>& queries, parlay::sequence<cluster<T,D>*>& answers) {


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

//do batch LCA on a forest (not just single tree)
//return nullptr if the query vertices are in separate trees (meaning that they can't have a common ancestor)
template<typename T, typename D>
void batchLCA(parlay::sequence<cluster<T,D>>& clusters,  parlay::sequence<std::tuple<T,T,T>>& queries, parlay::sequence<cluster<T,D>*>& answers) {

    int nc = clusters.size(); //vertices + edges~
    int k = queries.size(); //batch size

    //find which nodes are involved in the calculations
    parlay::sequence<T> involved_nodes;
    parlay::sequence<bool> marked_clusters(nc,false); 
    parlay::sequence<T> roots;
    forest_find_nodes_involved(clusters,k,queries,marked_clusters,involved_nodes,roots);
    int sn = involved_nodes.size();

    //index map maps between involved_nodes index and cluster index
    parlay::parlay_unordered_map<T,T> index_map = parlay::parlay_unordered_map<T,T>(2*sn); 

     parlay::parallel_for(0,sn,[&] (T i) {
        index_map.Insert(h.involved_nodes[i],i);
    });

    //store which root each involved_node belongs to
    parlay::sequence<T> root_ids(h.involved_nodes.size(),-1);

    populate_root_ids(clusters,marked_clusters,involved_nodes,roots,root_ids,index_map);

    parlay::sequence<T> rangk = parlay::tabulate(k,[&] (T i) {return i;});

    //TOD2* within a lambda function like this, better to write arugment type as reference? std::tuple<T,T,T>& query? Or does [&] already handle that?

    //find which queries all contain members in the same tree
    auto answerable_query_indices = parlay::filter(rangk,[&] (T i) {
        return root_ids[std::get<0>(queries[i])] == root_ids[std::get<1>(queries[i])] && root_ids[std::get<0>(queries[i])] == root_ids[std::get<2>(queries[i])];
    });
    //find which queries do not contain all members in same tree
    auto unanswerable_query_indices = parlay::filter(rangk,[&] (T i) {
        return root_ids[std::get<0>(queries[i])] != root_ids[std::get<1>(queries[i])] || root_ids[std::get<0>(queries[i])] != root_ids[std::get<2>(queries[i])];
    });
    //if query not connected, then answer is nullptr
    parlay::parallel_for(0,unanswerable_query_indices.size(),[&] (size_t i) {
        answers[i]=nullptr;
    });
    //creates sequence of key-sequence pairs of form
    //(ROOT, QUERY INDICES BELONGING TO THAT ROOT)
    auto queries_by_tree = parlay::group_by_key(parlay::map(answerable_query_indices,[&] (T i) {
        return std::make_pair(root_ids[std::get<0>(queries[i])],i);
    }));
    //in parallel, look at the queries for each root
    parlay::parallel_for(0,queries_by_tree.size(),[&] (size_t i) {
        //list of queries associated with root queries_by_tree[i].first
        query_list = parlay::map(queries_by_tree[i].second,[&] (T j) {return queries[j]});
        parlay::sequence<cluster<T,D>*> answer_list;

        //run LCA on these queries
        batchLCA_singletree(clusters,&clusters[queries_by_tree[i].first],query_list,answer_list);
        
        //store answers back in appropriate index
        parlay::parallel_for(0,answer_list.size(),[&] (size_t j) {
            answers[queries_by_tree[i].second[j]]=answer_list[j];
        });
    });


}


#endif //LCA