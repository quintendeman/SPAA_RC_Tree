//Code to handle LCA queries
//This function takes in input directly, converts to single tree, fixed root version, then calls the LCA code in fixedLCA.h

#ifndef LCAHH //! cannot define function of this name, good to know
#define LCAHH


//header copied from RC.h
#include "RC.h"
#include "cluster.h"
#include "parhash_include/parlay_hash/unordered_map.h" //for parlayhash concurrent hash map
#include "fixedLCA.h"

#include<bitset>

//whether or not to print debug info
const bool PRINT_L = false;

//given a set of roots and involved nodes, assign each involved node to a root
template<typename T, typename D>
void populate_root_ids(parlay::sequence<cluster<T,D>>& clusters, parlay::sequence<T>& involved_nodes, parlay::sequence<T>& roots,parlay::sequence<T>& root_ids,parlay::parlay_unordered_map<T,T>& index_map) {

    //can run roots independently because disconnected
    parlay::parallel_for(0,roots.size(),[&] (size_t roots_index) {
        parlay::sequence<T> stack; //stack holds CLUSTER index
        parlay::sequence<parlay::sequence<T>> new_stack;
        stack.push_back(roots[roots_index]); //index in clusters of each root
        while (stack.size() > 0) {
            new_stack = parlay::sequence<parlay::sequence<T>>(stack.size(),parlay::sequence<T>());
            parlay::parallel_for(0,stack.size(),[&] (size_t j) {
                T a = stack[j]; //current RC tree node index we are looking at
                //root ids holds ALT index
                root_ids[*index_map.Find(a)]=*index_map.Find(roots[roots_index]);
               
                //Add children to stack, to propagate down tree
                for (int l = 0; l < clusters[a].children.size(); l++) {
                    if (clusters[a].children[l] != nullptr && index_map.Find(clusters[a].children[l]->index) != std::nullopt) {
                        new_stack[j].push_back(clusters[a].children[l]->index);
                    }
                }
            });
            stack = parlay::flatten(new_stack);
        }
    });
}



template<typename T, typename D>
void forest_find_nodes_involved(parlay::sequence<cluster<T,D>>& clusters, int k, parlay::sequence<std::tuple<T,T,T>>& queries, parlay::sequence<T>& involved_nodes, parlay::sequence<T>& roots) {

    parlay::sequence<T> stack;
    parlay::sequence<T> new_stack(3*k,-1);

    parlay::parallel_for(0,k,[&] (size_t i) {
        //pointer math to avoid contention
        new_stack[3*i] = std::get<0>(queries[i]);
        new_stack[3*i+1] = std::get<1>(queries[i]);
        new_stack[3*i+2] = std::get<2>(queries[i]);
    });

    parlay::sequence<bool> add_to_initial_stack(new_stack.size(),false);
    parlay::parallel_for(0,new_stack.size(),[&] (size_t i) {
        //if we grab this elt first
        if (clusters[new_stack[i]].counter.fetch_add(1)==0) {
            add_to_initial_stack[i]=true;
        }

    });
    stack = parlay::pack(new_stack,add_to_initial_stack); 
    involved_nodes.append(stack);
    roots.append(parlay::filter(stack,[&] (T item) {return isNullary(&clusters[item]); }));

    //bottom-up, propagate true marks
    while(stack.size() > 0) {
        new_stack=parlay::sequence<T>(stack.size(),-1);
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            //a is an index in clusters
            auto a = stack[i];
            //load needed because parent is atomic
            if (!isNullary(&clusters[a]) && clusters[a].parent.load()->counter.fetch_add(1)==0) { //if we aren't at the root and we are the first child here
                new_stack[i]=clusters[a].parent.load()->index; //look @ in the next iteration
            }   

        });
        //don't pass up -1s (this is one way to filter, another way is to keep a filled_nodes bool and pack by it)
        stack=parlay::filter(new_stack,[&] (T item) {return item != -1;}); 
        involved_nodes.append(stack); //parallel insert the active nodes of this level 
        roots.append(parlay::filter(stack,[&] (T item) {return isNullary(&clusters[item]); }));
    }

    //reset all counters touched*
    parlay::parallel_for(0,involved_nodes.size(),[&] (size_t i) {
        clusters[involved_nodes[i]].counter=0;
    });
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
    if (PRINT_L) {
        std::cout << "query print" << std::endl;
        for (int i = 0; i < answers_uv.size(); i++) {
            std::cout << i << ": " << std::get<0>(queries[i]) << " " << std::get<1>(queries[i]) << " " << std::get<2>(queries[i]) << " " << answers_uv[i] << " " << answers_ur[i] << " " << answers_vr[i] << std::endl;
        }
    }
    if (PRINT_L) {
        std::cout << "Starting 1st query call u v" << std::endl;
        std::cout << "root of choice " << root->index << std::endl;
    }
    batch_fixed_LCA(clusters,root,queries_uv,answers_uv);
    if (PRINT_L) {
        std::cout << std::endl << std::endl << std::endl;
        std::cout << "Starting 2nd query call u r" << std::endl;
    }
    batch_fixed_LCA(clusters,root,queries_ur,answers_ur);
    if (PRINT_L) {
        std::cout << std::endl << std::endl << std::endl;
        std::cout << "Starting 3rd query call v r" << std::endl;
    }
    batch_fixed_LCA(clusters,root,queries_vr,answers_vr);
    if (PRINT_L) {
        std::cout << std::endl << std::endl << std::endl;
        std::cout << "answer print uv ur' vr'" << std::endl;
        for (int i = 0; i < answers_uv.size(); i++) {
            std::cout << i << ": " << std::get<0>(queries[i]) << " " << std::get<1>(queries[i]) << " " << std::get<2>(queries[i]) << " " << answers_uv[i] << " " << answers_ur[i] << " " << answers_vr[i] << std::endl;
        }
    }

    answers=parlay::tabulate(k,[&] (size_t i) {
        return &clusters[logic_lca(answers_uv[i],answers_ur[i],answers_vr[i])]; //this logic in Vanilla LCA
    });
    parlay::parallel_for(k,[&] (size_t i) {
        if (answers[i]->index != logic_lca(answers_uv[i],answers_ur[i],answers_vr[i])) {
            std::cout << "abort, answer index not the lca given" << std::endl;
            exit(7002);
        }
    });
    if (PRINT_L) std::cout << "exit" << std::endl;

}

//do batch LCA on a forest (not just single tree)
//return nullptr if the query vertices are in separate trees (meaning that they can't have a common ancestor)
template<typename T, typename D>
void batchLCA(parlay::sequence<cluster<T,D>>& clusters,  parlay::sequence<std::tuple<T,T,T>>& queries, parlay::sequence<cluster<T,D>*>& answers) {
    parlay::internal::timer t3;
    t3.start();

    int nc = clusters.size(); //vertices + edges~
    int k = queries.size(); //batch size

    if (nc > std::numeric_limits<T>::max()/2-10) {
        std::cout << "too many elements requested for T type, aborting" << std::endl;
        exit(10001);
    }

    //find which nodes are involved in the calculations
    parlay::sequence<T> involved_nodes;
    parlay::sequence<T> roots; //stores roots, by RC tree index
    forest_find_nodes_involved(clusters,k,queries,involved_nodes,roots);
    int sn = involved_nodes.size();

    if (PRINT_L) {
        std::cout << "printing queries " << std::endl;
        for (int i = 0; i < queries.size(); i++) {
            std::cout << std::get<0>(queries[i]) << " " << std::get<1>(queries[i]) << " " << std::get<2>(queries[i]) << std::endl;
        }

        std::cout << "printing involved nodes" << std::endl;
        for (int i = 0; i < involved_nodes.size(); i++) {
            std::cout << involved_nodes[i] << " ";
        }
        std::cout << std::endl;
    }
    //index map maps between involved_nodes index and cluster index
    parlay::parlay_unordered_map<T,T> index_map = parlay::parlay_unordered_map<T,T>(2*sn); 

     parlay::parallel_for(0,sn,[&] (T i) {
        index_map.Insert(involved_nodes[i],i);
    }); 

    if (PRINT_L) {
        std::cout << "roots" << std::endl;
        for (int i = 0 ; i < roots.size(); i++) {
            std::cout << roots[i] << std::endl;
        }
    }

    //store which root each involved_node belongs to
    //store in the language of the alt tree (because root_ids is indexed by the alt tree)
    parlay::sequence<T> root_ids(involved_nodes.size(),-1);

    populate_root_ids(clusters,involved_nodes,roots,root_ids,index_map);

    if (PRINT_L) {
        std::cout << "root ids: " << std::endl;
        for (int i = 0; i < root_ids.size(); i++) {
            std::cout << root_ids[i] << " " ;
        }
        std::cout << std::endl;
    }

    parlay::sequence<int> rangk = parlay::tabulate(k,[&] (int i) {return i;});

    if (PRINT_L) {
        std::cout << "printing index map for queries" << std::endl;
        for (int i = 0; i < queries.size(); i++) {
            std::cout << std::get<0>(queries[i]) << ": " << std::endl;
            std::cout << *index_map.Find(std::get<0>(queries[i])) << std::endl;
            std::cout << std::get<1>(queries[i]) << ": " << *index_map.Find(std::get<1>(queries[i])) << std::endl;
            std::cout << std::get<2>(queries[i]) << ": " << *index_map.Find(std::get<2>(queries[i])) << std::endl;

        }
    }


    //find which queries all contain members in the same tree
    auto answerable_query_indices = parlay::filter(rangk,[&] (int i) {
        return root_ids[*index_map.Find(std::get<0>(queries[i]))] == root_ids[*index_map.Find(std::get<1>(queries[i]))] && root_ids[*index_map.Find(std::get<0>(queries[i]))] == root_ids[*index_map.Find(std::get<2>(queries[i]))];
    });

    //find which queries do not contain all members in same tree
    auto unanswerable_query_indices = parlay::filter(rangk,[&] (int i) {
         return root_ids[*index_map.Find(std::get<0>(queries[i]))] != root_ids[*index_map.Find(std::get<1>(queries[i]))] || root_ids[*index_map.Find(std::get<0>(queries[i]))] != root_ids[*index_map.Find(std::get<2>(queries[i]))];
    });
    if (PRINT_L) {
        for (int i = 0; i < answerable_query_indices.size(); i++) {
            std::cout << "answerable " << i << " " << answerable_query_indices[i] << std::endl;
        }
        for (int i = 0; i < unanswerable_query_indices.size(); i++) {
            std::cout << "unanswerable " << i << " " << unanswerable_query_indices[i] << std::endl;
        }
    }
  
    //if query not connected, then answer is nullptr
    parlay::parallel_for(0,unanswerable_query_indices.size(),[&] (size_t i) {
        answers[i]=nullptr;
    });

    //creates sequence of key-sequence pairs of form
    //(ROOT, QUERY INDICES BELONGING TO THAT ROOT)
    auto queries_by_tree = parlay::group_by_key(parlay::map(answerable_query_indices,[&] (int i) {
        return std::make_pair(root_ids[*index_map.Find(std::get<0>(queries[i]))],i);
    }));
    if (PRINT_L) {
        for (int i = 0; i < queries_by_tree.size(); i++) {
            std::cout << i << " " << queries_by_tree[i].first << ": ";
            for (int j = 0; j < queries_by_tree[i].second.size();j++) {
                std::cout << queries_by_tree[i].second[j] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "about to single tree LCA" << std::endl;
    }
    
    std::string ps = "processed global info time " + std::to_string(t3.next_time()) + "\n";
    //std::cout << ps;

    //in parallel, look at the queries for each root 
    parlay::parallel_for(0,queries_by_tree.size(),[&] (size_t i) {

        //list of queries associated with root queries_by_tree[i].first
        parlay::sequence<std::tuple<T,T,T>> query_list = parlay::map(queries_by_tree[i].second,[&] (int j) {return queries[j];});
        parlay::sequence<cluster<T,D>*> answer_list;

        //run LCA on these queries
        batchLCA_singletree(clusters,&clusters[involved_nodes[queries_by_tree[i].first]],query_list,answer_list); 
        
        //store answers back in appropriate index
        parlay::parallel_for(0,answer_list.size(),[&] (size_t j) {
            answers[queries_by_tree[i].second[j]]=answer_list[j];
        });
    });
    


}


#endif //LCA