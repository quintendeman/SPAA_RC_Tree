#ifndef RC_H
#define RC_H

/**
 * Mainly deals with static RC tree creation
*/

#include <atomic>
#include <random>
#include <set>
#include <algorithm>
#include <iostream>
#include <mutex>
#include "../examples/helper/graph_utils.h"
#include <parlay/alloc.h>
#include "cluster.h"
#include "../examples/counting_sort.h"
#include <parlay/random.h>

static const char PRINT_QUERY = 0;


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
template<typename T>
parlay::sequence<T> generate_random_tree(T num_elements, int seed=-1) {
    assert(num_elements > 0);
    parlay::sequence<T> parents = parlay::tabulate(num_elements,[&] (T v) {return (T) 0;});

    std::mt19937 gen = get_rand_gen(seed);//what seed is being used here? TOD2 //check into (if num_elts = 1mil twice, will it give same graph?)

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

    //print sequence
    // for (int i = 0; i < num_elements; i++) {
    //     std::cout << parents[i] << " " << std::endl;
    // }
    // std::cout << std::endl;
    return parents;

}

/*
    Generate a simple, single rooted graph with each node having two children
    Then, randomly, change the parent of each node with a certain probability such that it picks something on the left of it

    Returns an array of parents such that the parents of index V would be parents[V]
*/
template <typename T>
parlay::sequence<T> generate_tree_graph(T num_elements)
{
    assert(num_elements > 0);

    parlay::sequence<T> dummy_initial_parents = parlay::tabulate(num_elements, [&] (T v) {return (T) 0;});

    parlay::random_generator gen(time(0)); //use time(0) as random seed
    std::uniform_real_distribution<double> dis(0, 1);

    parlay::parallel_for(0, num_elements, [&] (T v) {
        auto r = gen[v];
        double random_val = dis(r);

        static const double anywhere_left_weight = 0.1;
        static const double immediate_left_weight = 30; 
        static const double root_weight = 0.0; /* warning, does not represet probability of a null cluster as degree capping may create more forests */

        static const double anywhere_prob = (anywhere_left_weight/(anywhere_left_weight+immediate_left_weight+root_weight));
        static const double root_prob = anywhere_prob + (root_weight/(anywhere_left_weight+immediate_left_weight+root_weight));

        if (random_val <= anywhere_prob && v > 0)
        {
            std::uniform_int_distribution<T> disint(0,v-1);

            dummy_initial_parents[v] = disint(r);
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
void degree_cap_parents(parlay::sequence<T> &parents, const T max_degree)
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


/**
 * Extracts a particular bit (counted from the right) from an element
*/
template <typename T>
inline bool extract_bit(T number, int offset_from_right)
{
    return (number >> offset_from_right) & 1;
}

/*
    Returns index of first different bit from the left
    Sets the value of bit in the boolean bit
    Sets the value of b
    sizeof(T) must be less than 16 bytes
*/
template <typename T>
inline char first_different_bit(const T a, const T b, bool* bit)
{
    T difference = a ^ b;
    char num_bits = sizeof(T) * 8;
    
    for(char i = num_bits-1; i >= 0; i--)
    {
        bool inspected_bit = extract_bit(difference, i);
        if(inspected_bit)
        {
            if(bit) *bit = extract_bit(b, i);
            return i;
        }
    }

    return -1;
}

/*
    returns a char with I_w and C_w(I_w) packed
    I_w is the index of the first different bit
    C_w is the value in neighbour w of this bit
    Also sets the different bit index in the char ptr different_bit_index

    Also, technically I waste the left-most bit in each char
*/

template <typename T>
static unsigned char get_single_colour_contribution(const T vcolour, const T wcolour, char* different_bit_index = NULL)
{
    bool wbit = false;
    char different_bit = first_different_bit(vcolour, wcolour, &wbit);
    char final_returned_character = (different_bit << 1) | wbit;
    if(different_bit_index) *different_bit_index = different_bit;
    return final_returned_character;
}

// /*
//     Given a set of cluster ptrs, colours them.
//     The clusters must have an initial valid colouring stored in their temp_colour folder
//     The clusters should be for "vertices" s.t. one hop from each vertex is a cluster representing an edge
//     and two hops away is a vertex representing a neighbouring node
// */
template<typename T, typename D>
void colour_nodes(parlay::sequence<node<T,D>*> tree_nodes)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;    

    parlay::parallel_for(0, tree_nodes.size(), [&] (T I) {
        auto& cluster_ptr = tree_nodes[I]->cluster_ptr;
        unsigned long local_maximum = cluster_ptr->get_default_colour();
        unsigned long local_minimum = local_maximum;
        unsigned long my_colour = local_maximum;

        auto& node_ptr = tree_nodes[I];

        for(auto& edge_ptr : node_ptr->adjacents)
        {
            if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
                continue;
            cluster<T,D>* other_ptr = get_other_side(node_ptr, edge_ptr)->cluster_ptr;
            unsigned long compared_colour = other_ptr->get_default_colour();
            if(compared_colour > local_maximum)
                local_maximum = compared_colour;
            if(compared_colour < local_minimum)
                local_minimum = compared_colour;   
            
        }
        if(local_maximum == my_colour) // This node is a local maximum, give it a unique colour
        {
            cluster_ptr->colour = local_maximum_colour;
        }
        else if(local_minimum == my_colour)
        {
            cluster_ptr->colour = local_minimum_colour;
        }
        else
        {
            cluster_ptr->colour = 2 + (get_single_colour_contribution(my_colour, local_maximum) / 2); // adding 2 and removing indicator bit
        }
    });
    return;
}

// /*
//     sets a boolean flag in the clusters indicating that they're part of MIS
//     also may change the boolean flag of some other clusters, only consider the clusters in this
//     These clusters must have a maximum degree of 2
// */
template<typename T, typename D>
void set_MIS(parlay::sequence<node<T, D>*>& tree_nodes, bool use_tree_nodes = false)
{
    auto cluster_ptrs = parlay::tabulate(tree_nodes.size(), [&] (T i) {
        return tree_nodes[i]->cluster_ptr;
    });
    colour_nodes(tree_nodes);

    parlay::parallel_for(0, cluster_ptrs.size(), [&] (T v) {
        cluster_ptrs[v]->state &= (~IS_MIS_SET);
        if(use_tree_nodes)
            cluster_ptrs[v]->set_neighbour_mis(false, tree_nodes[v]);
        else
            cluster_ptrs[v]->set_neighbour_mis(false);
    });

    auto colours = parlay::tabulate(cluster_ptrs.size(), [&] (T v) {
        return cluster_ptrs[v]->colour;
    });

    auto vertices = parlay::tabulate(cluster_ptrs.size(), [] (T v) {
        return v;
    });

    auto result = vertices;

    parlay::sequence<unsigned long> offsets = counting_sort(vertices.begin(), vertices.end(), result.begin(), colours.begin(), 8 * sizeof(unsigned long));


    for(uint i = 0; i < offsets.size(); i++)
    {
        T start_index;
        if (i == 0)
            start_index = 0;
        else
            start_index = offsets[i-1];
        T end_index = offsets[i];

        parlay::parallel_for(start_index, end_index, [&] (T i) {
            T v = result[i];
            if(use_tree_nodes)
            {
                if(cluster_ptrs[v]->get_neighbour_mis(tree_nodes[v]) == true)
                {
                    cluster_ptrs[v]->state &= (~IS_MIS_SET);
                    return;
                } 
            }
            else if(cluster_ptrs[v]->get_neighbour_mis() == true)
            {
                cluster_ptrs[v]->state &= (~IS_MIS_SET);
                return;
            }
            cluster_ptrs[v]->state |= IS_MIS_SET;
        });
    }
}

template<typename T, typename D>
void finalize(node<T,D>* contracted_node)
{
    while(contracted_node->cluster_ptr->adjacency.get_tail() != contracted_node)
    {
        auto tail = contracted_node->cluster_ptr->adjacency.get_tail();
        if(tail->state & dont_finalize)
            break;
        contracted_node->cluster_ptr->adjacency.delete_tail();
    }

    return;
}


/**
 * MUST have at most 2 neighbours
 */
template<typename T, typename D>
void contract(node<T,D>* node_ptr, bool affect = false)
{
    
    if(node_ptr->get_num_neighbours_live() == 0)
    {
        node_ptr->state &= (~(unary_cluster | binary_cluster | nullary_cluster));
        node_ptr->state |= contracts_this_round | nullary_cluster;
        node_ptr->state &= (~live);
        node_ptr->state &= (~affected);
        node_ptr->cluster_ptr->parent = nullptr;
        if(affect)
        {
            for(auto& ptr : node_ptr->adjacents)
                ptr = nullptr;
            finalize(node_ptr);
        }
    }
    else if(node_ptr->get_num_neighbours_live() == 1)
    {
        node<T,D>* neighbour_node;
        node<T,D>* edge_node;
        for(auto& ptr : node_ptr->adjacents)
            if(ptr != nullptr && ptr->state & (base_edge | binary_cluster))
            {
                edge_node = ptr;
                break;
            }
        neighbour_node = get_other_side(node_ptr, edge_node);

        // do nothing for myself
        node_ptr->state &= (~live);
        node_ptr->state &= (~(unary_cluster | binary_cluster | nullary_cluster));
        node_ptr->state |= (unary_cluster | contracts_this_round);

        node_ptr->cluster_ptr->parent = neighbour_node->cluster_ptr;
        for(short k = 0; k < max_neighbours; k++)
        {
            auto& ptr = neighbour_node->adjacents[k];
            if(ptr == edge_node)
            {
                neighbour_node->cluster_ptr->children[k] = node_ptr->cluster_ptr; 
                break;
            }
        }

        edge_node->cluster_ptr->parent = node_ptr->cluster_ptr;
        

        for(short k = 0; k < max_neighbours; k++)
        {
            auto& ptr = node_ptr->adjacents[k];
            if(ptr == edge_node)
            {
                node_ptr->cluster_ptr->children[k] = edge_node->cluster_ptr; 
                break;
            }
        }

        for(auto& ptr : neighbour_node->adjacents)
            if(ptr == edge_node)
                ptr = node_ptr;
        if(affect == true)
        {
            finalize(node_ptr);
            finalize(edge_node);
            // for(auto& ptr : node_ptr->adjacents)
            //     ptr = nullptr;
            neighbour_node->state |= (affected | adjacency_changed | live);
            node_ptr->state &= (~affected);
            node_ptr->state &= (~adjacency_changed);
            neighbour_node->state &= (~(binary_cluster | base_edge));
        }
    }
    else if (node_ptr->get_num_neighbours_live() == 2)
    {
        node<T,D>* left_edge;
        node<T,D>* left_node;
        node<T,D>* right_edge;
        node<T,D>* right_node;
        for(auto& ptr : node_ptr->adjacents)
            if(ptr != nullptr && ptr->state & (base_edge | binary_cluster))
            {
                left_edge = ptr;
                left_node = get_other_side(node_ptr, left_edge);
                break;
            }
        for(auto& ptr : node_ptr->adjacents)
            if(ptr != nullptr && ptr->state & (base_edge | binary_cluster) && ptr != left_edge)
            {
                right_edge = ptr;
                right_node = get_other_side(node_ptr, right_edge);
                break;
            }
        node_ptr->state &= (~(unary_cluster | binary_cluster | nullary_cluster));
        node_ptr->state |= binary_cluster | contracts_this_round;
        node_ptr->state &= (~live);
        node_ptr->cluster_ptr->parent = nullptr;
        

        left_edge->cluster_ptr->parent = node_ptr->cluster_ptr;
        right_edge->cluster_ptr->parent = node_ptr->cluster_ptr;
        
        for (short k = 0; k < node_ptr->adjacents.size(); k++)
        {
            auto& ptr = node_ptr->adjacents[k];
            if (ptr == left_edge)
            {
                ptr = left_node;
                node_ptr->cluster_ptr->children[k] = left_edge->cluster_ptr;
            }
            if (ptr == right_edge)
            {
                node_ptr->cluster_ptr->children[k] = right_edge->cluster_ptr;
                ptr = right_node;
            }
        }

        for (short k = 0; k < left_node->adjacents.size(); k++)
        {
            auto& ptr = left_node->adjacents[k];
            if (ptr == left_edge)
            {
                ptr = node_ptr;
            }
        }

        for (short k = 0; k < right_node->adjacents.size(); k++)
        {
            auto& ptr = right_node->adjacents[k];
            if (ptr == right_edge)
            {
                ptr = node_ptr;
            }
        }


        if(affect == true)
        {
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != left_node && ptr != right_node)
                    ptr = nullptr;
            }
            right_node->state |= (affected | adjacency_changed | live);
            left_node->state |= (affected | adjacency_changed | live);
            right_node->state &= (~(binary_cluster | base_edge));
            left_node->state &= (~(binary_cluster | base_edge));
            node_ptr->state &= (~affected);
            node_ptr->state &= (~adjacency_changed);
            finalize(node_ptr);
            finalize(left_edge);
            finalize(right_edge);
        }
    }
}

/**
 * Given an symmetric graph, creates a set of clusters
 * In total, it creates n + m clusters in the array base_clusters
 * The first n are base_vertex clusters
 * it returns an array of n clusters corresponding to the original n vertices.
 * The edge clusters are present as edges.
*/
template <typename T, typename D>
void create_base_clusters(parlay::sequence<parlay::sequence<T>> &G, parlay::sequence<cluster<T, D> > &base_clusters, const T max_size, D defretval = 0.00f)
{

    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    base_clusters = parlay::tabulate(G.size(), [] (T i) {
        cluster<T,D> base_cluster;
        return base_cluster;
    });

    parlay::parallel_for(0, base_clusters.size(), [&] (T i){
        base_clusters[i].index = i;
        base_clusters[i].add_empty_level(live, 0);
    });


    parlay::parallel_for(0, base_clusters.size(), [&] (const T i){
        auto& clstr = base_clusters[i];
        clstr.state = base_vertex | live;
        clstr.data = defretval; //must reset base vertex clusters to default val as well as base edge clusters
        
        const auto& v= i;

        for(const auto& w : G[v])
        {
            if(w < v)
                continue;
            auto edge_cluster = cluster_allocator::alloc();
            edge_cluster->index = -1;
            edge_cluster->state = base_edge | live;
            edge_cluster->data = defretval;
            edge_cluster->add_empty_level(live | base_edge, 0);
            edge_cluster->add_ptr_to_highest_level(clstr.adjacency.get_tail());
            edge_cluster->add_ptr_to_highest_level(base_clusters[w].adjacency.get_tail());
            clstr.add_ptr_to_highest_level(edge_cluster->adjacency.get_tail());
        }
    });

    parlay::parallel_for(0, G.size(), [&] (T v) {
        auto cluster_ptr = &base_clusters[v];
        cluster_ptr->index = v;
        for(const auto& w : G[v])
        {
            if(v < w) // reversed
                continue;
            // find the edge that corresponds/joins to the other side
            const auto& node_ptr_arr = base_clusters[w].adjacency.get_tail()->adjacents;

            node<T,D>* edge_ptr = nullptr;
            bool found = false;
            for(auto& ptr : node_ptr_arr)
            {
                if(ptr != nullptr)
                {
                    const auto& final_adjacency_ptr_arr = ptr->adjacents;
                    for(const auto& v_ptr : final_adjacency_ptr_arr)
                    {
                        if(v_ptr != nullptr && v_ptr->cluster_ptr == cluster_ptr)
                        {    
                            edge_ptr = ptr;
                            found = true;
                            break;
                        }
                    }
                }
                if(found)
                    break;
            }
            base_clusters[v].add_ptr_to_highest_level(edge_ptr);
        }
    });

    return;
}


/**
 * Make sure the base clusters are consistent i.e. that every neighbour of mine has me as their neighbour too
 */
template <typename T, typename D>
void check_consistency(parlay::sequence<node<T, D>*>& tree_nodes)
{
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        if((node_ptr->state & live) == 0)
        {
            return;
        }
        if(node_ptr->state & (unary_cluster | binary_cluster | base_edge))
        {
            return;
        }

        for(auto& ptr : node_ptr->adjacents)
        {
            if(ptr == nullptr || !(ptr->state & (base_edge | binary_cluster)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, ptr);
            if(get_other_side(other_node_ptr, ptr) != node_ptr)
            {
                node_ptr->cluster_ptr->print();
                ptr->cluster_ptr->print();
                other_node_ptr->cluster_ptr->print();
                std::cout << "Base clusters inconsistent" << std::endl;
                exit(1);
            }
        }
    });

    return;
}

template <typename T, typename D>
void check_mis(parlay::sequence<node<T, D>*>& tree_nodes)
{
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        for(auto& ptr : node_ptr->adjacents)
        {
            if(ptr == nullptr || !(ptr->state & (base_edge | binary_cluster)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, ptr);
            if(other_node_ptr->cluster_ptr->state & IS_MIS_SET && node_ptr->cluster_ptr->state & IS_MIS_SET )
            {
                std::cout << red << "Neighbours are not MIS!" << reset << std::endl;
                exit(1);
            }
        }
    });

    return;
}

template <typename T, typename D>
void check_parents_children(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        auto cluster_ptr = &clusters[i];
        if(cluster_ptr->parent != nullptr)
        {
            bool in_parent = false;
            for(auto& par_child : cluster_ptr->parent->children)
            {
                if(par_child == cluster_ptr)
                    in_parent = true;
            }
            if(in_parent == false)
            {
                std::cout << red << "RC tree parents inconsistent!" << reset << std::endl;
                cluster_ptr->print();
                cluster_ptr->parent->print();
                exit(1);
            }
        }
        for(auto& child : cluster_ptr->children)
        {
            if(child == nullptr)
                continue;
            if(child->parent != cluster_ptr)
            {
                std::cout << red << "RC tree children inconsistent!" << reset << std::endl;
                cluster_ptr->print();
                child->print();
                exit(1);
            }
        }
    });
}

template<typename T, typename D>
void check_children_values(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        auto cluster_ptr = &clusters[i];
        D value_from_children = 0.0;
        for(auto& child : cluster_ptr->children)
        {
            if(child == nullptr)
                continue;
            value_from_children+=child->data;
        }
        if(value_from_children != cluster_ptr->data)
        {
            std::cout << red << "RC tree values inconsistent!" << reset << std::endl;
            cluster_ptr->print();
            exit(1);
        }
    });
}

/**
 * Make an exact copy of the current adjacency list, copying over everything (including the state), the only new thing will be the level.
 * When making a new round, pointers to edges will also get renewed, is that ideal? probably not but we don't have a choice do we 
 */
template<typename T, typename D>
void recreate_last_levels(parlay::sequence<node<T,D>*>& tree_nodes)
{
    
    
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        auto& cluster_ptr = node_ptr->cluster_ptr;
        auto& v = cluster_ptr->index;
        cluster_ptr->add_empty_level(node_ptr->state, node_ptr->contraction_level+1);

        for(short e = 0; e < node_ptr->adjacents.size(); e++)
        {
            auto edge_ptr = node_ptr->adjacents[e];
            if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, edge_ptr);
            auto& w = other_node_ptr->cluster_ptr->index;
            if(w < v)
                continue;
            edge_ptr->cluster_ptr->add_empty_level(edge_ptr->state, edge_ptr->contraction_level + 1); 
            node_ptr->next->adjacents[e] = edge_ptr->next;
            // edge_ptr->next->add_ptr(node_ptr->next);
            for(short k = 0; k < edge_ptr->adjacents.size(); k++)
            {
                if(edge_ptr->adjacents[k] == node_ptr)
                    edge_ptr->next->adjacents[k] = node_ptr->next;
            }
        }
        // for(auto& edge_ptr : node_ptr->adjacents)
        // {
        //     if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
        //         continue;
        //     auto other_node_ptr = get_other_side(node_ptr, edge_ptr);
        //     auto& w = other_node_ptr->cluster_ptr->index;
        //     if(w < v)
        //         continue;
            
        //     // std::cout << "Binary added " << edge_ptr->index() << std::endl;

        //     edge_ptr->cluster_ptr->add_empty_level(edge_ptr->state, edge_ptr->contraction_level + 1); 
        //     cluster_ptr->add_ptr_to_highest_level(edge_ptr->next);
        //     edge_ptr->cluster_ptr->add_ptr_to_highest_level(node_ptr->next);
        // }
    });
    parlay::parallel_for(0, tree_nodes.size(), [&] (T i) {
        auto& node_ptr = tree_nodes[i];
        auto& cluster_ptr = node_ptr->cluster_ptr;
        auto& v = cluster_ptr->index;

        for(short e = 0; e < node_ptr->adjacents.size(); e++)
        {
            auto edge_ptr = node_ptr->adjacents[e];
            if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
                continue;
            auto other_node_ptr = get_other_side(node_ptr, edge_ptr);
            auto& w = other_node_ptr->cluster_ptr->index;
            if(v < w)
                continue;
            node_ptr->next->adjacents[e] = edge_ptr->next;
            // edge_ptr->next->add_ptr(node_ptr->next);
            for(short k = 0; k < edge_ptr->adjacents.size(); k++)
            {
                if(edge_ptr->adjacents[k] == node_ptr)
                    edge_ptr->next->adjacents[k] = node_ptr->next;
            }
        }

        // for(auto& edge_ptr : node_ptr->adjacents)
        // {
            
        //     if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
        //         continue;
    
        //     auto other_node_ptr = get_other_side(node_ptr, edge_ptr);
        //     auto& w = other_node_ptr->cluster_ptr->index;
        //     if(v < w)
        //         continue;
        //     cluster_ptr->add_ptr_to_highest_level(edge_ptr->next);
        //     edge_ptr->cluster_ptr->add_ptr_to_highest_level(node_ptr->next);
        // }
    });

    parlay::parallel_for(0, tree_nodes.size(), [&] (T i){
        auto& node_ptr = tree_nodes[i];
        auto& cluster_ptr = node_ptr->cluster_ptr;
        for(short e = 0; e < node_ptr->adjacents.size(); e++)
        {
            auto edge_ptr = node_ptr->adjacents[e];
            if(edge_ptr == nullptr)
                continue;
            if(edge_ptr == nullptr || !(edge_ptr->state & (binary_cluster | base_edge)))
                continue;
            if(edge_ptr->state & unary_cluster && !(edge_ptr->state & (binary_cluster | base_edge)))
            {
                edge_ptr->cluster_ptr->add_empty_level(unary_cluster, edge_ptr->contraction_level+1);
                edge_ptr->next->state = unary_cluster;
                node_ptr->next->adjacents[e] = edge_ptr->next;
            }
        }
        tree_nodes[i] = tree_nodes[i]->next;
    });
}


/**
 * The main workhorse, populates the set base_clusters with internal_clusters
 * As base_vertices aren't useful, it replaces them in place
 * 
 * The accumulation is simple for now -- each child points to the parent cluster
 * This can be used for connectivity tracking by starting at base_clusters[v] and base_clusters[w] and going to the parents
 * until they overlap
 * 
 * A more complex accumulator may be but currently the edges have no weights
 * This takes O(n) work and O(log n log n) span instead of O(log n log log n) span because the filter operatio is exact
 * not an approximate compaction
 * 
 * defretval is the default data stored in a cluster
*/

template <typename T, typename D>
void create_RC_tree(parlay::sequence<cluster<T,D> > &base_clusters, T n, bool randomized = false, D defretval = 0.00f, bool print=true)
{


    auto tree_nodes = parlay::tabulate(base_clusters.size(), [&] (T i){
        return base_clusters[i].adjacency.get_tail();
    });

    // if(base_clusters.size() <= 100)
    //     printTree(base_clusters);

    recreate_last_levels(tree_nodes);

    auto count = 0;
    do
    {
        check_consistency(tree_nodes); // TODO remove
        // break;
        auto eligible = parlay::filter(tree_nodes, [] (auto node_ptr){
            return node_ptr->get_num_neighbours_live() <= 2;
        });

        
        set_MIS(eligible);
        auto candidates = parlay::filter(eligible, [] (auto node_ptr){
            return node_ptr->cluster_ptr->state & IS_MIS_SET;
        });

        check_mis(candidates); // TODO remove

        
        parlay::parallel_for(0, candidates.size(), [&] (T i){
            candidates[i]->cluster_ptr->first_contracted_node = candidates[i]->prev;
            contract(candidates[i]);
        });

        
        tree_nodes = parlay::filter(tree_nodes, [] (auto node_ptr) {
            return node_ptr->state & live;
        });

        
        recreate_last_levels(tree_nodes);

        // if(base_clusters.size() <= 100)
        //     printTree(base_clusters);
    
        if (print)
        std::cout << "[" << count  << "]: " << bold << tree_nodes.size() << reset << std::endl;
        count++;
    }while(tree_nodes.size() > 0 && count < 100);

}



template <typename T, typename D, typename assocfunc>
D PathQuery( cluster<T, D>* v,  cluster<T, D>* w, const D& defretval, assocfunc func)
{
    if(v == w)
        return defretval;

    bool use_v = false;

    cluster<T,D>* prev_boundary_v_l = nullptr;
    cluster<T,D>* prev_boundary_v_r = nullptr;
    D prev_val_till_v_l = defretval;
    D prev_val_till_v_r = defretval;
    
    cluster<T,D>* prev_boundary_w_l = nullptr;
    cluster<T,D>* prev_boundary_w_r = nullptr;
    D prev_val_till_w_l = defretval;
    D prev_val_till_w_r = defretval;


    while(true)
    {
        // Ascend using the lower-height cluster. If either has to ascend above root, they are not connected.
        if(v == nullptr && w == nullptr && prev_boundary_v_l == nullptr && prev_boundary_v_r == nullptr)
            return defretval;   
        else if(v == w)
            break;
        else if (v == nullptr && w != nullptr)
            use_v = false;
        else if (w == nullptr && v != nullptr)
            use_v = true;
        else if(v->get_height() < w->get_height())
            use_v = true;
        else
            use_v = false;
        
        if(use_v)
        {
            cluster<T,D>* &v_cluster = v;

            D lval_print,rval_print = defretval; // TODO remove

            // We are at the start probably can add a bunch of nullptr checks too TODO
            if (prev_val_till_v_l == defretval && prev_val_till_v_r == defretval) 
            {
                v_cluster->find_boundary_vertices(prev_boundary_v_l, prev_val_till_v_l, prev_boundary_v_r, prev_val_till_v_r, defretval);
                lval_print = prev_val_till_v_l;
                rval_print = prev_val_till_v_r;
            } 
            else 
            {
                D lval;
                D rval;
                cluster<T,D>* l;
                cluster<T,D>* r;

                v_cluster->find_boundary_vertices(l, lval, r, rval, defretval);

                // Covers both unary to unary and unary to binary case
                if (prev_boundary_v_l == v && prev_boundary_v_r == v) 
                {
                    lval_print = prev_val_till_w_l;
                    rval_print = prev_val_till_w_r;
                    
                    if(l != r) // ascending into unary
                    {
                        prev_val_till_v_l = func(prev_val_till_v_l, lval);
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                        
                    }
                    else
                    {
                        auto v_node_ptr = v_cluster->first_contracted_node;
                        for(auto& ptr : v_node_ptr->adjacents)
                        {
                            if(ptr != nullptr && ptr->state & (binary_cluster | base_edge))
                            {
                                prev_val_till_v_l = func(prev_val_till_v_l, ptr->cluster_ptr->data);
                                prev_val_till_v_r = func(prev_val_till_v_r, ptr->cluster_ptr->data);

                            }
                        }
                    }
                    prev_boundary_v_l = l;
                    prev_boundary_v_r = r;
                }
                // Covers binary ascending to unary/nullary
                else if (l == r) 
                {
                    if(prev_boundary_v_r == v)
                    {
                        prev_val_till_v_r = prev_val_till_v_l;
                    }
                    else
                        prev_val_till_v_l = prev_val_till_v_r;

                    prev_boundary_v_l = prev_boundary_v_r = l;
                }
                // Binary to binary
                else 
                {
                    if (prev_boundary_v_l == l) 
                    {
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    } 
                    else if (prev_boundary_v_r == l) 
                    {
                        std::swap(prev_boundary_v_r, prev_boundary_v_l);
                        std::swap(prev_val_till_v_l, prev_val_till_v_r);
                        
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    
                    } 
                    else if (prev_boundary_v_l == r) 
                    {
                        std::swap(r, l);
                        std::swap(rval, lval);    
                        
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    } 
                    else // prev_boundary_v_r == r 
                    {     
                        prev_boundary_v_l = l;
                        prev_val_till_v_l = func(prev_val_till_v_l, lval);
                    }
                }
                lval_print = lval;
                rval_print = rval;
            }
            std::cout << bright_magenta << v->index << "[" << (short) v->get_height() << "]" << " ";
            std::cout << prev_boundary_v_l->index << "(" << prev_val_till_v_l << ") " << prev_boundary_v_r->index << "(" << prev_val_till_v_r << ") ";
            std::cout << "lval,rval=" << lval_print << "," << rval_print << " ";
            auto& for_printing = v;
            if(for_printing->adjacency.get_tail()->state & unary_cluster)
                std::cout << " unary ";
            else if (for_printing->adjacency.get_tail()->state & binary_cluster)
                std::cout << " binary ";
            else
                std::cout << " nullary ";
            std::cout << reset << std::endl;
            v = v->get_parent();
        }
        else
        {
            cluster<T,D>* &w_cluster = w;
            D lval_print,rval_print = defretval;

            // We are at the start
            if (prev_val_till_w_l == defretval && prev_val_till_w_r == defretval) 
            {
                w_cluster->find_boundary_vertices(prev_boundary_w_l, prev_val_till_w_l, prev_boundary_w_r, prev_val_till_w_r, defretval);
                lval_print = prev_val_till_w_l;
                rval_print = prev_val_till_w_r;
            } 
            else 
            {
                D lval, rval;
                cluster<T,D>* l;
                cluster<T,D>* r;
                w_cluster->find_boundary_vertices(l, lval, r, rval, defretval);

                // Covers both unary to unary and unary to binary case
                if (prev_boundary_w_l == w && prev_boundary_w_r == w) 
                {

                    if(l != r)
                    {
                        prev_val_till_w_l = func(prev_val_till_w_l, lval);
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    }
                    else
                    {
                        auto w_node_ptr = w_cluster->first_contracted_node;
                        for(auto& ptr : w_node_ptr->adjacents)
                        {
                            if(ptr != nullptr && ptr->state & (binary_cluster | base_edge))
                            {
                                prev_val_till_w_l = func(prev_val_till_w_l, ptr->cluster_ptr->data);
                                prev_val_till_w_r = func(prev_val_till_w_r, ptr->cluster_ptr->data);

                            }
                        }
                    }
                    prev_boundary_w_l = l;
                    prev_boundary_w_r = r;
                }
                // Covers binary ascending to unary/nullary
                else if (l == r) 
                {
                    if(prev_boundary_w_r == w)
                    {
                        prev_val_till_w_r = prev_val_till_w_l;
                    }
                    else
                        prev_val_till_w_l = prev_val_till_w_r;

                    prev_boundary_w_l = prev_boundary_w_r = l;
                }
                // Binary
                else 
                {
                    if (prev_boundary_w_l == l) 
                    {
                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    } 
                    else if (prev_boundary_w_r == l) 
                    {
                        
                        std::swap(prev_boundary_w_l, prev_boundary_w_r);
                        std::swap(prev_val_till_w_l, prev_val_till_w_r);

                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);

                    } 
                    else if (prev_boundary_w_l == r) 
                    {
                        std::swap(l, r);
                        std::swap(lval, rval);

                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    } 
                    else // prev_boundary_w_r == r 
                    { 
                        prev_boundary_w_l = l;
                        prev_val_till_w_l = func(prev_val_till_w_l, lval);
                    }
                }
                lval_print = lval;
                rval_print = rval;
            }
            std::cout << bright_green << w->index << "[" << (short) w->get_height() << "]" << " ";
            std::cout << prev_boundary_w_l->index << "(" << prev_val_till_w_l << ") " << prev_boundary_w_r->index << "(" << prev_val_till_w_r << ") ";
            std::cout << "lval,rval=" << lval_print << "," << rval_print << " ";
            auto& for_printing = w;
            if(for_printing->adjacency.get_tail()->state & unary_cluster)
                std::cout << " unary ";
            else if (for_printing->adjacency.get_tail()->state & binary_cluster)
                std::cout << " binary ";
            else
                std::cout << " nullary ";
            std::cout << reset << std::endl;
            w = w->get_parent();
        }
        
    }
    std::cout << std::endl;

    if(prev_boundary_v_l == nullptr && prev_boundary_w_l == nullptr)
    {
        return defretval;
    }
    if(prev_boundary_v_l == nullptr)
    {
        if(w == prev_boundary_w_r)
            return prev_val_till_w_r;
        else
            return prev_val_till_w_l;   
    }
    if(prev_boundary_w_l == nullptr)
    {
        if(v == prev_boundary_v_l)
            return prev_val_till_v_l;
        else
            return prev_val_till_v_r;
    }
    
    D v_contrib, w_contrib;
    
    if(v == prev_boundary_v_l)
        v_contrib = prev_val_till_v_l;
    else
        v_contrib = prev_val_till_v_r;

    if(w == prev_boundary_w_l)
        w_contrib = prev_val_till_w_l;
    else
        w_contrib = prev_val_till_w_r;

    return func(v_contrib, w_contrib);

}



template<typename T, typename D>
node<T,D>* get_edge(T v, T w, parlay::sequence<cluster<T,D>>& clusters, unsigned char level = 0)
{
    if(v == w)
        return nullptr;
    auto& v_cluster = clusters[v];
    auto v_node_ptr = v_cluster.adjacency.get_head();
    for(auto& edge_ptr : v_node_ptr->adjacents)
    {
        if(edge_ptr == nullptr)
            continue;
        auto w_node_ptr = get_other_side(v_node_ptr, edge_ptr);
        if(w_node_ptr->cluster_ptr->index == w)
            return edge_ptr;
    }

    auto& w_cluster = clusters[w];
    auto w_node_ptr = w_cluster.adjacency.get_head();
    for(auto& edge_ptr : w_node_ptr->adjacents)
    {
        if(edge_ptr == nullptr)
            continue;
        v_node_ptr = get_other_side(w_node_ptr, edge_ptr);
        if(v_node_ptr->cluster_ptr->index == v)
            return edge_ptr;
    }
    return nullptr;
}

// Use chatgpt to convert from recursive to iterative,
// who cares
template<typename T, typename D>
D manual_subtree_sum(cluster<T, D>* root, cluster<T, D>* dir_giver, parlay::sequence<cluster<T,D>>& clusters) {
    D ret_val = 0.0;
    if (root == nullptr || dir_giver == nullptr || root == dir_giver)
        return 0.0;

    // Stack to store pairs of clusters and their corresponding parent
    std::stack<std::pair<cluster<T, D>*, cluster<T, D>*>> stack;
    stack.push({root, dir_giver});

    while (!stack.empty()) {
        auto [current_root, parent] = stack.top();
        stack.pop();

        T r = current_root->index;

        for (auto& edge_ptr : current_root->adjacency.get_head()->adjacents) {
            if (edge_ptr == nullptr)
                continue;

            auto other_cluster = get_other_side(current_root->adjacency.get_head(), edge_ptr)->cluster_ptr;
            if (other_cluster == parent)
                continue;

            T other_index = other_cluster->index;
            ret_val += get_edge(r, other_index, clusters)->cluster_ptr->data;

            // Push the next cluster to the stack with the current cluster as its parent
            stack.push({other_cluster, current_root});
        }
    }

    return ret_val;
}



template<typename T, typename D, typename assocfunc>
D subtree_query(cluster<T, D>* root, cluster<T, D>* dir_giver, D defretval, assocfunc func)
{
    D ret_val = defretval;
    if(root == nullptr || dir_giver == nullptr || root == dir_giver)
        return defretval;
    if(PRINT_QUERY)
        std::cout << bold << bright_yellow << "Root: " << root->index << " dir: " << dir_giver->index << reset << std::endl;

    D lval = defretval;
    D rval = defretval;
    cluster<T,D>* dir_l = nullptr;
    cluster<T,D>* dir_r = nullptr;

    dir_giver->find_boundary_vertices(dir_l, lval, dir_r, rval, defretval);

    cluster<T,D>* root_l = nullptr;
    cluster<T,D>* root_r = nullptr;

    root->find_boundary_vertices(root_l, lval, root_r, lval, defretval);

    if(PRINT_QUERY)
    {
        std::cout << bold << yellow << "Root's boundary vertices ";
        if(root_l == nullptr)
            std::cout << "nl ";
        else
            std::cout << root_l->index << " ";
        if(root_r == nullptr)
            std::cout << "nl ";
        else
            std::cout << root_r->index << " ";
        std::cout << reset << std::endl;

        std::cout << bold << yellow << "dir_giver's boundary vertices ";
        if(dir_l == nullptr)
            std::cout << "nl ";
        else
            std::cout << dir_l->index << " ";
        if(dir_r == nullptr)
            std::cout << "nl ";
        else
            std::cout << dir_r->index << " ";
        std::cout << reset << std::endl;
    }

    if(root == dir_l || root == dir_r)
    {
        if(PRINT_QUERY)
            std::cout << bold << bright_cyan << "root already boundary vertex of dir giver" << reset << std::endl;

        auto child_with_dir = dir_giver;
        while(child_with_dir->parent != root)
            child_with_dir = child_with_dir->parent;

        auto& child_that_contains_u = child_with_dir; // rename

        cluster<T,D>* ign_l = nullptr;
        cluster<T,D>* ign_r = nullptr;
        child_that_contains_u->find_boundary_vertices(ign_l, lval, ign_r, rval, defretval);

        if(PRINT_QUERY)
            std::cout << bold << bright_cyan << "Child to ignore is " << child_with_dir->index << reset << std::endl;

        bool first = true;

        for(const auto& child : root->children)
        {
            if(child == nullptr)
                continue;

            cluster<T,D>* child_l = nullptr;
            cluster<T,D>* child_r = nullptr;
            if(child->adjacency.get_head()->state & base_edge)
                child->find_endpoints(child_l, child_r);
            else
                child->find_boundary_vertices(child_l, lval, child_r, rval, defretval);
            if(PRINT_QUERY)
            {
                std::cout << bold << bright_cyan << "[" << root->index << "] child " << child->index << " boundary vertices ";
                if(child_l == nullptr)
                    std::cout << "nl ";
                else
                    std::cout << child_l->index << " ";
                if(child_r == nullptr)
                    std::cout << "nl ";
                else
                    std::cout << child_r->index << " ";
                std::cout << reset << std::endl;
            }

            if(child == child_with_dir)
                continue;
            
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "Adding value of " << child->index << "/" << child->data << reset << std::endl;
            
            if(first)
            {
                first = false;
                ret_val = child->data;
            }
            else
            {
                ret_val = func(ret_val, child->data);
            }
        }

        if(root_l != nullptr && root_l != ign_l && root_l != ign_r)
        {
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "recursive call to  " << root_l->index << " wrt " << root->index << reset << std::endl;
            auto child_ret_val = subtree_query(root_l, root, defretval, func);
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "returning from recursive call to  " << root_l->index << " wrt " << root->index << ", adding " << child_ret_val << reset << std::endl;
            ret_val = func(ret_val, child_ret_val);
        }
        if(root_r != root_l && root_r != nullptr && root_r != ign_l && root_r != ign_r)
        {
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "recursive call to  " << root_r->index << " wrt " << root->index << reset << std::endl;
            auto child_ret_val = subtree_query(root_r, root, defretval, func);
            if(PRINT_QUERY)
                std::cout << bold << bright_cyan << "returning from recursive call to  " << root_r->index << " wrt " << root->index << ", adding " << child_ret_val << reset << std::endl;
            ret_val = func(ret_val, child_ret_val);    
        }
    }
    else
    {
        if(PRINT_QUERY)
            std::cout << bold << bright_green << "Should only go here once" << reset << std::endl;

        bool first = true;
        for(auto& child : root->children)
        {
            if(child == nullptr)
                continue;
            
            cluster<T,D>* child_l = nullptr;
            cluster<T,D>* child_r = nullptr;
            if(child->adjacency.get_head()->state & base_edge)
                child->find_endpoints(child_l, child_r);
            else
                child->find_boundary_vertices(child_l, lval, child_r, rval, defretval);
            if(PRINT_QUERY)
            {
                std::cout << bold << bright_green << "[" << root->index << "] child " << child->index << " boundary vertices ";
                if(child_l == nullptr)
                    std::cout << "nl ";
                else
                    std::cout << child_l->index << " ";
                if(child_r == nullptr)
                    std::cout << "nl ";
                else
                    std::cout << child_r->index << " ";
                std::cout << reset << std::endl;
            }

            if((child->adjacency.get_head()->state & base_edge) || (child->first_contracted_node->next->state & binary_cluster))
            {
                if(child_l == dir_giver || child_r == dir_giver)
                    continue;
                if(first)
                {
                    first = false;
                    ret_val = child->data;
                }
                else
                {
                    ret_val = func(ret_val, child->data);
                }
                if(child_l != root)
                {               
                    auto child_ret_val = subtree_query(child_l, root, defretval, func);
                    ret_val = func(ret_val, child_ret_val);
                }
                if(child_l != child_r && child_r != root)
                {
                    auto child_ret_val = subtree_query(child_r, root, defretval, func);
                    ret_val = func(ret_val, child_ret_val);
                }
            }
            else
            {
                if(first)
                {
                    first = false;
                    ret_val = child->data;
                }
                else
                {
                    ret_val = func(ret_val, child->data);
                }
            }
            
        }
    }


    return ret_val;
}


template<typename T, typename D, typename assocfunc>
void batchModifyEdgeWeights(const parlay::sequence<std::tuple<T, T, D>>& edges, assocfunc func, parlay::sequence<cluster<T, D>>& clusters, const D& defretval = 0.00f)
{   

    // TODO merge  with up pass
    auto edge_ptrs = parlay::tabulate(edges.size(), [&] (T i){
        T v = std::get<0>(edges[i]);
        T w = std::get<1>(edges[i]);
        D data = std::get<2>(edges[i]);

        auto edge_ptr = get_edge(v, w, clusters);
        if(edge_ptr == nullptr) // TODO REMOVE
        {
            std::cout << red << "Got a nullptr" << reset << std::endl;
            exit(1);
        }
        edge_ptr->cluster_ptr->data = data;

        auto ret_ptr = edge_ptr->cluster_ptr;
        auto cluster_ptr = ret_ptr;

        cluster_ptr->counter = 0;

        while(cluster_ptr != nullptr)
        {
            if(cluster_ptr->counter.fetch_add(1))
                break;
            cluster_ptr = cluster_ptr->get_parent();
        }


        return ret_ptr;
    });



    parlay::parallel_for(0, edge_ptrs.size(), [&] (T i) {
        if(edge_ptrs[i] == nullptr)
            return;

        auto cluster_ptr = edge_ptrs[i]; // base edge, technically not a cluster dw
        
        auto start = edge_ptrs[i];

        // decrement counter
        while(cluster_ptr != nullptr)
        {
            if(cluster_ptr->counter.fetch_add(-1) > 1)
                break;
            if(cluster_ptr != start)
            {
                D final_val = defretval; 
                node<T,D>* node_ptr = cluster_ptr->first_contracted_node;
                bool first_child = true;

                for(auto& child : cluster_ptr->children)
                {
                    if(child == nullptr)
                        continue;
                    
                    if(first_child)
                    {
                        final_val = child->data;
                        first_child = false;
                    }
                    else
                        final_val = func(final_val, child->data);
                
                }             
                cluster_ptr->data = final_val;
            }
            cluster_ptr = cluster_ptr->get_parent();
        }
    });

    return;
}

template<typename T, typename D>
void check_counter(parlay::sequence<cluster<T,D>>& clusters)
{
    parlay::parallel_for(0, clusters.size(), [&] (T i) {
        if(clusters[i].counter)
        {
            std::cout << red << "Counter wasn't reset" << reset << std::endl;
            exit(1);
        }
    });
}

/**
 * Only frees the "floating" edge clusters
*/
template<typename T, typename D>
void deleteRCtree(parlay::sequence<cluster<T, D>> &base_clusters)
{
    using cluster_allocator = parlay::type_allocator<cluster<T,D>>;

    parlay::parallel_for(0, base_clusters.size(), [&] (T i){
        auto& clstr = base_clusters[i];
        auto& node_ptr = clstr.first_contracted_node;
        if(node_ptr->state & binary_cluster)
        {
            node_ptr = node_ptr->prev;
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != nullptr && ptr->state & base_edge)
                {
                    cluster_allocator::destroy(ptr->cluster_ptr);
                }
            }
        }
        else
        {
            for(auto& ptr : node_ptr->adjacents)
            {
                if(ptr != nullptr && ptr->state & base_edge)
                {
                    cluster_allocator::destroy(ptr->cluster_ptr);
                }
            }
        }
    });
}

template<typename T, typename D>
void printTree(parlay::sequence<cluster<T, D>> &base_clusters, unsigned char level = -1)
{
    for(uint i = 0; i < base_clusters.size(); i++)
    {
        base_clusters[i].print(level);
    }   
}

#endif