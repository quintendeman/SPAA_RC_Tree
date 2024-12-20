
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


//set the subtree size of each vertex
//example of bottom up computation
template<typename T>
void set_size_par(parlay::sequence<parlay::sequence<T>>& child_tree, parlay::sequence<T>& parent_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {

    //rangn
    parlay::sequence<T> rangn = parlay::tabulate(child_tree.size(),[&] (T i) {return i;});

    //get out the leaves, we start here
    parlay::sequence<T> leaves = parlay::filter(rangn,[&] (size_t i) {
        return child_tree[i].size()==0;
    });

    //set default size value
    parlay::parallel_for(0,augmented_vertices.size(),[&] (size_t i) {
        augmented_vertices[i].size=1; 
    });

    parlay::sequence<T> stack = parlay::tabulate(leaves.size(),[&] (size_t i) {return leaves[i];});
    parlay::sequence<T> new_stack(stack.size(),-1);
    parlay::sequence<bool> filled_indices(stack.size(),false);


    while (stack.size() > 0) {
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            auto s = stack[i];
            //add child sizes
            for (T j = 0; j < child_tree[s].size(); j++) {
                augmented_vertices[s].size += augmented_vertices[child_tree[s][j]].size;
            }

            //if this node is the 0th child of its parent, let this node add its parent to the stack, otherwise not
            if (child_tree[parent_tree[s]][0]==s) {
                new_stack[i]=parent_tree[s];
                filled_indices[i]=true;

            }
        });
        
        stack = parlay::pack(new_stack,filled_indices); //the stack is smaller than before, because there are probably less nodes in level l-1 than level l
        new_stack=parlay::tabulate(stack.size(),[&] (size_t i){return -1;});
        filled_indices=parlay::tabulate(stack.size(),[&] (size_t i) {return false;});
        
    }

}


//set the subtree size of each vertex
//example of bottom up computation
//tag stack with child to avoid overwrites (from leaves being at different tree levels)
//requires level tagging to be complete (set_level_par must have already been called)
template<typename T>
void set_size_alt_par(parlay::sequence<parlay::sequence<T>>& child_tree, parlay::sequence<T>& parent_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {

    //get the max level in the tree
    //TOD2* compute this once, outside (for constant factor efficiency)
    int max_level = *parlay::max_element(parlay::map(augmented_vertices,[&] (LCANode<T> node) {
        return node.level;
    }));

    //rangn
    parlay::sequence<T> rangn = parlay::tabulate(child_tree.size(),[&] (T i) {return i;});


    //set default size value
    parlay::parallel_for(0,augmented_vertices.size(),[&] (size_t i) {
        augmented_vertices[i].size=1; 
    });
    //mark parts of stack that are still relevant
    parlay::sequence<bool> filled_indices(stack.size(),false);


    parlay::sequence<T> bottom_level = parlay::filter(rangn,[&] (size_t i) {
        return augmented_vertices[i].level==max_level;
    });
    parlay::sequence<std::pair<T,T>> stack = parlay::tabulate(bottom_level.size(),[&] (size_t i) {return std::make_pair(bottom_level[i],-1); });

    parlay::sequence<std::pair<T,T>> new_stack(stack.size(),std::make_pair(-1,-1));

    while (stack.size() > 0) {
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            auto s = stack[i].first;
            auto chi = stack[i].second;
            if (chi != -1) {
                //note that because we traverse from deepest to shallow, the child will be complete when it is added to the size of the parent
                augmented_vertices[s].size += augmented_vertices[chi].size;

            }
    
            //if we have not reached the root yet, add this edge to the stack
            if (parent_tree[s] != s) {
                new_stack[i]=std::make_pair(parent_tree[s],s);
                filled_indices[i]=true;

            }
            
        });
        
        stack = parlay::pack(new_stack,filled_indices); //the stack is smaller than before, because shallower levels hit the root first
        
        //reset the filler variables
        new_stack=parlay::tabulate(stack.size(),[&] (size_t i){return std::make_pair(-1,-1);});
        filled_indices=parlay::tabulate(stack.size(),[&] (size_t i) {return false;});
        
    }

}


//set the subtree size of each vertex
//example of bottom up computation
//requires level tagging to be complete (set_level_par must have already been called)
template<typename T>
void set_size_improved_par(parlay::sequence<parlay::sequence<T>>& child_tree, parlay::sequence<T>& parent_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices, int max_level) {

  
    //rangn
    parlay::sequence<T> rangn = parlay::tabulate(child_tree.size(),[&] (T i) {return i;});


    //set default size value
    parlay::parallel_for(0,augmented_vertices.size(),[&] (size_t i) {
        augmented_vertices[i].size=1; 
    });
    //mark parts of stack that are still relevant
    parlay::sequence<bool> filled_indices(stack.size(),false);


    parlay::sequence<T> bottom_level = parlay::filter(rangn,[&] (size_t i) {
        return augmented_vertices[i].level==max_level;
    });
    parlay::sequence<T> stack = parlay::tabulate(bottom_level.size(),[&] (size_t i) {return bottom_level[i]; });

    parlay::sequence<T> new_stack(stack.size(),-1);

    for (int iter = max_level; iter > 0; iter--) {
        //invariant: in the iteration iter, all of the vertices in the stack have level iter
        parlay::parallel_for(0,stack.size(),[&] (size_t i) {
            auto s = stack[i];
            for (T j = 0; j < child_tree[s].size(); j++) {
                augmented_vertices[s].size += augmented_vertices[child_tree[s][j]].size;
            }

            //if this node is the 0th child of its parent, let this node add its parent to the stack, otherwise not
            //TOD2* critically depends on order that children are added to child tree - if the backedge to the parent is added first, then this check does not work
            if (child_tree[parent_tree[s]][0]==s) {
                new_stack[i]=parent_tree[s];
                filled_indices[i]=true;

            }

            //if we have not reached the root yet, add this edge to the stack
            if (parent_tree[s] == s) {
                std::cout << "error, reached root early" << std::endl;
                exit(803);
            }
    
            new_stack[i]=parent_tree[s];
            filled_indices[i]=true;

            
        });
        
        stack = parlay::pack(new_stack,filled_indices); //the stack is smaller than before, because shallower levels hit the root first
        
        //reset the filler variables
        new_stack=parlay::tabulate(stack.size(),[&] (size_t i){return -1;});
        filled_indices=parlay::tabulate(stack.size(),[&] (size_t i) {return false;});
        
    }

}


//note that inlabel does not need a tree recursive computation, because of interval property of preorder
template<typename T>
void set_inlabel_par(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& av) {

    parlay::parallel_for(0,av.size(),[&] (size_t node) {
        int i = floor(log2((av[node].preorder-1) ^ (av[node].preorder+av[node].size-1)));
        av[node].inlabel=(1<<i)*floor((av[node].preorder+av[node].size-1)/(1<<i));

    });

}



//set the ascendant value for each vertex
//top down computation
template<typename T>
void set_ascendant_par(parlay::sequence<T>& parent_tree, parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {

    T l = floor(log2(parent_tree.size())); //l (from paper) -- level depth of full binary tree with n nodes

    parlay::sequence<T> stack;
    parlay::sequence<parlay::sequence<T>> new_stack;
    stack.push_back(root);

    while (stack.size() > 0) {

        new_stack = parlay::sequence(stack.size(), parlay::sequence<T>());

        parlay::parallel_for(0,stack.size(),[&] (size_t j) {
            auto s = stack[j];

            if (s == root) {
            augmented_vertices[s].ascendant = 1 << l; //2^l
            }
            else {
                //same inlabel, just pass down same ascendant value
                if (augmented_vertices[s].inlabel == augmented_vertices[parent_tree[s]].inlabel) {
                    augmented_vertices[s].ascendant = augmented_vertices[parent_tree[s]].ascendant; 
                }
                //different inlabel, add increment to ascendant value of s
                else {
                
                    T i = r1(augmented_vertices[s].inlabel);

                    augmented_vertices[s].ascendant = augmented_vertices[parent_tree[s]].ascendant | (1 << i); //or instead of +


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
    int max_level = *parlay::max_element(parlay::map(augmented_vertices,[&] (LCANode<T> node) {
        return node.level;
    }));

    set_size_improved_par(child_tree,parent_tree,root,augmented_vertices,max_level);
    set_preorder_par(child_tree,root,augmented_vertices);

    set_inlabel_par(child_tree,root,augmented_vertices);
    set_ascendant_par(parent_tree,child_tree,root,augmented_vertices);
    set_head_par(head,augmented_vertices,parent_tree);


} 




#endif