//encode the LCA finding method from Schieber & Vishkin (1988)
//this version does naively (not trying to be efficient or ||)
#ifndef SSLCA
#define SSLCA

#include<math.h>

//T is indexing type
template<typename T>
struct LCAnode {
    T id; //id of node (ex. array index)
    T inlabel; 
    T preorder; //position visited in preorder traversal (lchild yourself rchild)
    T num_rightmost_zeroes; //of preorder
    T level; //level in tree (0 is root, leaves are higher)
    T ascendant; //composition of inlabel of all ancestors
    T size;

    void print() {
        std::cout << this->id << ":\t" << this->inlabel << ",\t"  << this->preorder << ",\t" << this->num_rightmost_zeroes << ",\t" << this->level << ",\t" << this->ascendant  << ",\t" << this->size << std::endl;
    }
    
};


//get index of rightmost 1
template<typename T> 
int r1(T x) {

    return static_cast<int>(log2(x-(x&(x-1))));

}

//get index of leftmost 1
//note: cannot plug bitset in here, must plug in number type
template<typename T>
int l1(T x) {
    if ( x < 0) {
        std::cout << "abort getting leftmost 1 of negative #" << std::endl;
        exit(501);
    }
    return floor(log2(x));
}


template<typename T>
void print_augmented_vertices(parlay::sequence<LCAnode<T>>& av) {
    std::cout << "printing augment" << std::endl;
    std::cout << "id,\t inl,\t preo,\t#ri0,\tlvl,\tascd,\t size" << std::endl;
    for (int i =0 ; i < av.size(); i++) {
        av[i].print(); 
    }


}

template<typename T>
void preprocess(parlay::sequence<T>& parent_tree, parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices, parlay::sequence<T>& head) {


    set_preorder(child_tree,root,augmented_vertices);
    set_level(child_tree,root,augmented_vertices);
    set_size(child_tree,root,augmented_vertices);
    set_inlabel(child_tree,root,augmented_vertices);
    set_ascendant(parent_tree,child_tree,root,augmented_vertices);
    set_head(head,augmented_vertices,parent_tree);


} 

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

//traverse the tree in order and mark the preorder field in LCAnode
//tree in children format
//child_tree : tree, in child format
//root : root of tree
//augmented vertices : list of the additional info stored with vertices (preorder, level, etc.)
template<typename T>
void set_preorder(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {

    std::deque<T> stack; //use stack instead of recursion

    stack.push_back(root);
    T count = 1; //TOD2* start at 1 or 0?

    while (stack.size() > 0) {
        T s = stack.back();
        stack.pop_back();

        augmented_vertices[s].preorder = count;
        count += 1;
        //TOD2 does order of accessing children matter here?
        for (T i = child_tree[s].size()-1; i >= 0; i--) {
            stack.push_back(child_tree[s][i]);
        }
    }

}

//set the level (depth in tree) of each vertex
template<typename T>
void set_level(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {
    //pair is (index id, depth)
    std::deque<std::pair<T,T>> stack; //use stack instead of recursion

    stack.push_back(std::make_pair(root,0));

    while (stack.size() > 0) {
        auto s = stack.back();
        stack.pop_back();

        augmented_vertices[s.first].level = s.second;

        for (T i = 0; i < child_tree[s.first].size(); i++) {
            stack.push_back(std::make_pair(child_tree[s.first][i], s.second+1));
        }
    }

}


//set the subtree size of each vertex
template<typename T>
void set_size(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {

    //keeps track if we've seen a node before
    parlay::sequence<bool> first_time(child_tree.size(),true);

    //set initial size to 0
    parlay::parallel_for(0,child_tree.size(), [&] (T i) {
        augmented_vertices[i].size=1;
    });

    std::deque<T> stack; //use stack instead of recursion

    stack.push_back(root);

    while (stack.size() > 0) {
        auto s = stack.back();
        stack.pop_back();
        if (!first_time[s]) {
            for (T i = 0; i < child_tree[s].size(); i++) {
                augmented_vertices[s].size += augmented_vertices[child_tree[s][i]].size;
            }
            
        }
        else {
            first_time[s]=false; //mark that we've seen this node before
            stack.push_back(s); //push back the parent again
            //but look at the children first
            for (T i = 0; i < child_tree[s].size(); i++) {
                stack.push_back(child_tree[s][i]);
            }

        }
        
    }

}


//TOD2 can make much more efficient
//do bottom up computation
template<typename T>
void set_inlabel(parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {
     //keeps track if we've seen a node before
    parlay::sequence<bool> first_time(child_tree.size(),true);

    //calculate the # of rightmost zeroes for each preorder value
    parlay::parallel_for(0,child_tree.size(), [&] (size_t s) {
        augmented_vertices[s].num_rightmost_zeroes = r1(augmented_vertices[s].preorder);

    });

    std::deque<T> stack; //use stack instead of recursion

    stack.push_back(root);

    while (stack.size() > 0) {
        auto s = stack.back();
        stack.pop_back();
        if (!first_time[s]) {            
            //set this vertex's inlabel to be the "max" of the children
            augmented_vertices[s].inlabel = get_inlabel(child_tree,augmented_vertices,s);
         
            std::cout << "set inlabel of " << s << std::endl;
            
        }
        else {
            first_time[s]=false; //mark that we've seen this node before
            stack.push_back(s); //push back the parent again
            //but look at the children first
            for (T i = 0; i < child_tree[s].size(); i++) {
                stack.push_back(child_tree[s][i]);
            }

        }
    }

}

//get the correct inlabel for s, *assuming that all children have had their inlabel correctly assigned*
template<typename T>
T get_inlabel(parlay::sequence<parlay::sequence<T>>& child_tree, parlay::sequence<LCAnode<T>>& augmented_vertices, T s) {
    int inlabel = augmented_vertices[s].preorder;
    int max_zeroes_seen = r1(augmented_vertices[s].preorder);
    for (T i = 0; i < child_tree[s].size(); i++) {
        int cand = r1(augmented_vertices[child_tree[s][i]].inlabel);
        if (cand > max_zeroes_seen) {
            max_zeroes_seen = cand;
            inlabel = augmented_vertices[child_tree[s][i]].inlabel;
        }
    }
    return inlabel;

}


//set the ascendant value for each vertex
template<typename T>
void set_ascendant(parlay::sequence<T>& parent_tree, parlay::sequence<parlay::sequence<T>>& child_tree, T root, parlay::sequence<LCAnode<T>>& augmented_vertices) {

    std::deque<T> stack;
    stack.push_back(root);

    T l = floor(log2(parent_tree.size())); //l (from paper) -- level depth of full binary tree with n nodes

    while (stack.size() > 0) {
        auto s = stack.back();
        stack.pop_back();
        //if we're looking at the root
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
        //add children to stack!
        for (T i = 0; i < child_tree[s].size(); i++) {
            stack.push_back(child_tree[s][i]);
        }


    }


}

//set the head. requires head to be previously allocated (to size n)
template<typename T>
void set_head(parlay::sequence<T>& head, parlay::sequence<LCAnode<T>>& augmented_vertices, parlay::sequence<T>& parent_tree) {
    for (int i = 0; i < parent_tree.size(); i++) {
        //if the inlabel of vertex i is NOT the same as its parent, this means that vertex i is the head of a path (path of same inlabels partition the tree) -- mark this in head
        //root is also a HEAD by def
        if (augmented_vertices[i].inlabel != augmented_vertices[parent_tree[i]].inlabel || i==parent_tree[i]) {
            head[augmented_vertices[i].inlabel] = i;

        }
    }
}

//Step 1 of implementation (Section 4 Schieber & Vishkin)
//imax is a substep in finding b, the binary tree ancestor of the inlabels of x and y
//case1 -- INLABEL(u) ancestor of INLABEL(v)
//case2 -- INLABEL(v) ancestor of INLABEL(u)
//case3 -- INLABEL(u) and INLABEL(v) not ancestors of each other
template<typename T>
std::pair<int,T> get_imax(parlay::sequence<LCAnode<T>>& av, T u, T v) {
    int i1 = r1(av[u].inlabel);
    int i2 = r1(av[v].inlabel);
    int i3 = static_cast<int>(log2(av[u].inlabel ^ av[v].inlabel));
    if (i1 >= i2 && i1 >= i3) {
        return std::make_pair(i1,u);
    }
    return std::make_pair(std::max(i2,i3),v);
   
}

//Step 2 of implementation (Section 4 Schieber & Vishkin)
//separating out for readability
template<typename T>
std::pair<int,int> get_inlabelz(parlay::sequence<LCAnode<T>>& av, T u, T v, int imax) {
    int common = av[u].ascendant & av[v].ascendant;
    //std::cout << "common is " << common << std::endl;
    int ci = (1 << imax) * (common >> imax);
    int jz = r1(ci);
    //std::cout << "j index is  " << jz << std::endl;
    int mask = -1 ^ ((1 << jz)-1);
    return std::make_pair((av[u].inlabel & mask) | (1 << jz),jz);
}

//Step 3 of implementation (section 4 Schieber & Vishkin)
template<typename T>
int get_hat(parlay::sequence<T>& head, parlay::sequence<T>& parent_tree,parlay::sequence<LCAnode<T>>& av, T u, int inlabel_z, int jz) {
    int xhat=-1;
    if (av[u].inlabel == inlabel_z) {
        //std::cout << "inlabels same easy hat" << std::endl;
        xhat=u;
    }
   // additional case not mentioned in paper!? is this needed? TOD2*
    else if (av[parent_tree[u]].inlabel == inlabel_z) {
        //std::cout << "parent is xhat" << std::endl;
        xhat=parent_tree[u];
    }
    else {
        int w = -1;
        int trimmed_ancestor = (1 << jz) - 1;
        //std::cout << "trimmed ancestor : " << trimmed_ancestor << std::endl;
        int k = l1(av[u].ascendant & trimmed_ancestor); //TOD2* is the xor with inlabelz needed?  //((1 << l1(inlabel_z))-1)
        //int k = l1(av[u].ascendant);
        //std::cout << "get hat k is " << k << std::endl;
        int mask =-1 ^ ((1 << k)-1);
        int inlabel_w = (av[u].inlabel & mask) | (1 << k);
       // std::cout << "in label w is " << inlabel_w << std::endl;
        w = head[inlabel_w];
        xhat = parent_tree[w];

    }
    return xhat;

}


//given preprocessed data structure, find an LCA
//av = augmented_vertices
template<typename T>
T query(parlay::sequence<T>& head, parlay::sequence<T>& parent_tree,parlay::sequence<LCAnode<T>>& av, T u, T v) {
    //then u and v are on the same path, so whichever has lower level is the LCA
    if (av[u].inlabel == av[v].inlabel) {
        if (av[u].level <= av[v].level) {
            return u;
        }
        return v;
    }
    //if u and v are not on the same path
    auto i_pair = get_imax(av,u,v);
    int imax = i_pair.first;

    //std::cout << "imax is " << imax << std::endl;

    int mask = -1 ^ ((1 << imax) - 1); //zeroes in positions smaller than imax, 1s everywhere else
    //i_pair.second is the index of the vertex from which we take the base inlabel (the l-i leftmost bits)
    //std::cout << "mask: " << mask << std::endl;
    //std::cout << "chosen inlabel: " << av[i_pair.second].inlabel << std::endl;
    int b = ((av[i_pair.second].inlabel & mask) & (-1 - (1 << imax) + 1)) | (1 << imax);
   // std::cout << "b is " << b << std::endl;
    //TODO* test b is correct; test z is correct; etc.

    auto output = get_inlabelz(av,u,v,imax);
    int inlabel_z = output.first;
    int jz = output.second;


   // std::cout << "z's in label is " << inlabel_z << std::endl;

  
    int xhat=get_hat(head,parent_tree,av,u,inlabel_z,jz);
    int yhat=get_hat(head,parent_tree,av,v,inlabel_z,jz);

    //std::cout << "uhat is " << xhat << ", " << "vhat is " << yhat << std::endl;

    if (av[xhat].level <= av[yhat].level) {
        return xhat;
    }
    return yhat;

}

// //TOD2 parallelize this
// void init_augment(parlay::sequence<LCAnode<T>>& augmented_vertices) {
//     for (int i = 0; i < augmented_vertices.size(); i++) {
//         augmented_vertices[i].id = i;
//         augmented_vertices.inlabel = -1;
//         augmented_vertices.inorder = -1;
//         augmented_vertices.preorder = -1;
//         augmented_vertices.level = -1;
//         augmented_vertices.ascendant = -1;
//         augmented_vertices.size = 0;

//     }
// }

#endif