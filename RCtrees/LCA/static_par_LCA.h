
//encode the LCA finding method from Schieber & Vishkin (1988)
//this version is trying to do in parallel, unlike the sequential static_seqn_LCA.h
//note that the query is the same, only the preprocessing is in parallel

//note: scrolled through Tom Tseng Euler Tour github + code, cite*? TOD2*
//https://github.com/tomtseng/parallel-euler-tour-tree/tree/master/src/dynamic_trees/euler_tour_tree

#ifndef STATIC_PAR_LCA
#define STATIC_PAR_LCA

#include<math.h>
#include "static_seqn_LCA.h"
#include "VanillaLCA.h"

//traverse the tree in order and mark the preorder field in LCAnode
//top-down calculation (parallel scan)
//requires size of subtrees to be already calculated
//tree in children format
//child_tree : tree, in child format
//root : root of tree
//augmented vertices : list of the additional info stored with vertices (preorder, level, etc.)
template<typename T>
void set_preorder_par(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {

    parlay::sequence<T> stack; //stack from which we draw tasks
    parlay::sequence<parlay::sequence<T>> new_stack; 

    //set initial stack val
    stack.push_back(root);
    augmented_vertices[root].preorder=1;

    while (stack.size() > 0) {

        // print_augmented_vertices(augmented_vertices); //DEBUGGING
        // std::cout << "printing stack (set preorder ) " << std::endl;
        // for (T i : stack) {
        //     std::cout << i << " ";
        // }
        // std::cout << std::endl;

        new_stack = parlay::sequence(stack.size(), parlay::sequence<T>());

        parlay::parallel_for(0,stack.size(),[&] (size_t i) {

            auto s = stack[i];
            int running_count = augmented_vertices[s].preorder+1;
            for (T j = 0; j < child_tree[s].size(); j++) {
                augmented_vertices[child_tree[s][j]].preorder = running_count;
                running_count += augmented_vertices[child_tree[s][j]].size;
                //have a stack of stacks to avoid races/concurrent inserts
                new_stack[i].push_back(child_tree[s][j]);

            }

        });

        //copy over new stack, but flattened for easier access
        stack = parlay::flatten(new_stack);
       
    }

}


//traverse the tree in order and mark the preorder field in LCAnode
//top-down calculation (parallel scan)
//requires size of subtrees to be already calculated
//tree in children format
//child_tree : tree, in child format
//root : root of tree
//augmented vertices : list of the additional info stored with vertices (preorder, level, etc.)
template<typename T>
void set_preorder_par2(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {

    parlay::sequence<T> stack; //stack from which we draw tasks
    parlay::sequence<parlay::sequence<T>> new_stack; 

    //set initial stack val
    stack.push_back(root);
    augmented_vertices[root].preorder=1;

    while (stack.size() > 0) {

        //print_augmented_vertices(augmented_vertices); //DEBUGGING

        new_stack = parlay::sequence(stack.size(), parlay::sequence<T>());

        parlay::parallel_for(0,stack.size(),[&] (size_t i) {

            auto s = stack[i];
            int running_count = augmented_vertices[s].preorder+1;
            for (T j = 0; j < child_tree[s].size(); j++) {
                augmented_vertices[child_tree[s][j]].preorder = running_count;
                running_count += augmented_vertices[child_tree[s][j]].size;
                //have a stack of stacks to avoid races/concurrent inserts
                new_stack[i].push_back(child_tree[s][j]);

            }

        });

        //copy over new stack, but flattened for easier access
        stack = parlay::flatten(new_stack);
       
    }

}





//set the level (depth in tree) of each vertex
//example of top-down computation
//TOD2* is this how you do a top down computation?
template<typename T>
void set_level_par(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {
    //pair is (index id, depth)

    parlay::sequence<std::pair<T,T>> stack; //stack from which we draw tasks
    stack.push_back(std::make_pair(root,0));

    parlay::sequence<parlay::sequence<std::pair<T,T>>> new_stack; 

    while (stack.size() > 0) {

        //print_augmented_vertices(augmented_vertices); //DEBUGGING

        new_stack = parlay::sequence(stack.size(), parlay::sequence<std::pair<T,T>>());
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {

            auto s = stack[i];
            
            augmented_vertices[s.first].level=s.second; //set level
            for (T j = 0; j < child_tree[s.first].size(); j++) {
                //have a stack of stacks to avoid races/concurrent inserts
                new_stack[i].push_back(std::make_pair(child_tree[s.first][j], s.second+1));

            }

        });

        //copy over new stack, but flattened for easier access
        stack = parlay::flatten(new_stack);
       
    }

}


//given a set of vertices, separate them into a sequence of sequences of levels (ith level in ith index)
//use a top-down computation; this is like setting the level, except that we are recording the stack and not doing anything to the nodes of each level
//this is needed for a proper bottom-up computation
//TOD2* return by reference not copy?
template<typename T>
parlay::sequence<parlay::sequence<T>> level_batch(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices, int max_level) {

    parlay::sequence<parlay::sequence<T>> level_record(max_level+1,parlay::sequence<T>());

    parlay::sequence<T> stack; //stack from which we draw tasks
    stack.push_back(root);

    parlay::sequence<parlay::sequence<T>> new_stack; 

    for (int i = 0; i < max_level+1; i++) {

        level_record[i]=parlay::tabulate(stack.size(),[&] (size_t j) {return stack[j];});

        new_stack = parlay::sequence(stack.size(), parlay::sequence<T>());
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {

            auto s = stack[i];
            
            for (T j = 0; j < child_tree[s].size(); j++) {
                //have a stack of stacks to avoid races/concurrent inserts
                new_stack[i].push_back(child_tree[s][j]);

            }

        });

        //copy over new stack, but flattened for easier access
        stack = parlay::flatten(new_stack);
       
    }
    return level_record;

}


//set the subtree size of each vertex
//example of bottom up computation
//requires level tagging to be complete (set_level_par must have already been called)
//TOD2* note to self: having a level record a much cleaner way of doing sweeps! 
template<typename T>
void set_size_par(parlay::sequence<parlay::sequence<T>>& child_tree, parlay::sequence<T>& parent_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices, int max_level) {

  
    //set default size value
    parlay::parallel_for(0,augmented_vertices.size(),[&] (size_t i) {
        augmented_vertices[i].size=1; 
    });

    //do top-down sweep to get record of which nodes on which levels in sequence of sequence form
    auto level_record = level_batch(child_tree,root,augmented_vertices,max_level);

    // std::cout << "printing level record" << std::endl;
    // for (int i = 0; i < level_record.size(); i++) {
    //     std::cout << i << ": ";
    //     for (T j : level_record[i]) {
    //         std::cout << j << " ";
    //     }
    //     std::cout << std::endl;
    // }

    //std::cout << "Debugging in set size" << std::endl;
    for (int iter = max_level; iter >= 0; iter--) {

        //print_augmented_vertices(augmented_vertices); //DEBUGGING
       

        //invariant: in the iteration iter, all of the vertices in the stack have level iter
        parlay::parallel_for(0,level_record[iter].size(),[&] (size_t i) {
            auto s = level_record[iter][i];
            for (T j = 0; j < child_tree[s].size(); j++) {
                augmented_vertices[s].size += augmented_vertices[child_tree[s][j]].size;
            }
            
        });
        
    }

}


//note that inlabel does not need a tree recursive computation, because of interval property of preorder
template<typename T>
void set_inlabel_par(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& av) {
    long long one=1;
    parlay::parallel_for(0,av.size(),[&] (size_t node) {
        int i = floor(log2((av[node].preorder-1) ^ (av[node].preorder+av[node].size-1)));
        av[node].inlabel=(one<<i)*floor((av[node].preorder+av[node].size-1)/(one<<i));

    });

}



//set the ascendant value for each vertex
//top down computation
template<typename T>
void set_ascendant_par(parlay::sequence<T>& parent_tree, parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {
    long long one =1;

    T l = floor(log2(parent_tree.size())); //l (from paper) -- level depth of full binary tree with n nodes

    parlay::sequence<T> stack;
    parlay::sequence<parlay::sequence<T>> new_stack;
    stack.push_back(root);

    while (stack.size() > 0) {

        new_stack = parlay::sequence(stack.size(), parlay::sequence<T>());

        parlay::parallel_for(0,stack.size(),[&] (size_t j) {
            auto s = stack[j];

            if (s == root) {
            augmented_vertices[s].ascendant = one << l; //2^l
            }
            else {
                //same inlabel, just pass down same ascendant value
                if (augmented_vertices[s].inlabel == augmented_vertices[parent_tree[s]].inlabel) {
                    augmented_vertices[s].ascendant = augmented_vertices[parent_tree[s]].ascendant; 
                }
                //different inlabel, add increment to ascendant value of s
                else {
                
                    T i = r1(augmented_vertices[s].inlabel);

                    augmented_vertices[s].ascendant = augmented_vertices[parent_tree[s]].ascendant | (one << i); //or instead of +


                }

            }

            for (T k = 0; k < child_tree[s].size(); k++) {
                new_stack[j].push_back(child_tree[s][k]);
            }


        });

        stack=parlay::flatten(new_stack);
        
    }


}

//set the head. requires head to be previously allocated (to size n)
template<typename T>
void set_head_par(parlay::sequence<T>& head, parlay::sequence<LCAnode<T>>& augmented_vertices, parlay::sequence<T>& parent_tree) {
    //can loop over in parallel because each node marks at most one (distinct) head location
    parlay::parallel_for(0,parent_tree.size(),[&] (size_t i) {
          //if the inlabel of vertex i is NOT the same as its parent, this means that vertex i is the head of a path (path of same inlabels partition the tree) -- mark this in head
        //root is also a HEAD by def
        if (augmented_vertices[i].inlabel != augmented_vertices[parent_tree[i]].inlabel || i==parent_tree[i]) {
            head[augmented_vertices[i].inlabel] = i;

        }


    });
    
}


template<typename T>
void preprocess_par(parlay::sequence<T>& parent_tree, parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices, parlay::sequence<T>& head) {


    set_level_par(child_tree,root,augmented_vertices);

    //get the max level in the tree
    //TOD2* compute this once, outside (for constant factor efficiency)
    int max_level = *parlay::max_element(parlay::map(augmented_vertices,[&] (LCAnode<T> node) {
        return node.level;
    }));

    //std::cout << "max level is " << max_level << std::endl;

    set_size_par(child_tree,parent_tree,root,augmented_vertices,max_level);
    set_preorder_par(child_tree,root,augmented_vertices);

    set_inlabel_par(child_tree,root,augmented_vertices);
    set_ascendant_par(parent_tree,child_tree,root,augmented_vertices);
    set_head_par(head,augmented_vertices,parent_tree);


} 




#endif