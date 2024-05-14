#include <atomic>
#include <random>
#include <set>
#include <iostream>
#include <mutex>
#include "cluster.h"


#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

const short empty_type = 0;
const short base_vertex = 1;
const short base_edge = 2;
const short unary_cluster = 4;
const short binary_cluster = 8;
const short nullary_cluster = 16;
const short live = 256;
const short carrying_weight = 512;
const short internal = 8192;



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

    std::mt19937 gen(std::random_device{}());

    parlay::parallel_for(0, num_elements, [&] (T v) {
        double lambda = 1.0 / (((double) v) * 0.1);
        std::exponential_distribution<double> dist(lambda);
        double value = dist(gen);
        T T_value = (T) value;
        if (T_value > v)
            T_value = v;
        dummy_initial_parents[v] = T_value;
    });  
    
    return dummy_initial_parents;
}

/*
    Converts the parents array into an assymetric graph
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
        if(counts[parents[v]].load() >= (max_degree))
        {
            parents[v] = v;
        }
        T parent_count = counts[parents[v]].fetch_add(1);
        if(parent_count >= (max_degree - 1))
        {
            parents[v] = v;
        }
    });

}




/**
 * Deletes assymetric pairs in a nested sequence representing a graph
*/
template <typename vertex>
void delete_assymetric_pairs(parlay::sequence<parlay::sequence<vertex>>& G)
{
     auto vertices = parlay::tabulate<vertex>(G.size(), [&] (vertex i) {return i;});

     parlay::sequence<parlay::sequence<bool> >  keep_edges_graph = parlay::map(vertices, [&] (vertex v) {
        auto edge = G[v];
        auto keep_edge = parlay::tabulate<bool>(edge.size(), [&] (vertex i) {return (bool) false;});
        
        vertex starting_node = v;
        parlay::parallel_for(0, edge.size(), [&] (vertex e){
            vertex ending_node = edge[e];
            parlay::sequence<vertex> ending_edge_list = G[ending_node];
            
            // is starting node in ending node?
            parlay::parallel_for(0, ending_edge_list.size(), [&] (vertex w) {
                if (ending_edge_list[w] == starting_node)
                    keep_edge[e] = true;
            });
        });
        return keep_edge;
     });

     parlay::parallel_for(0, G.size(), [&] (vertex v) {
        auto edge = G[v];
        parlay::parallel_for(0, edge.size(), [&] (vertex w){
            edge[w] = keep_edges_graph[v][w] ? edge[w] : -1;
        }
        );

        edge = parlay::filter(edge, [&] (vertex w){
            return w != -1;
        });

        G[v] = edge;
     });

     
}


/**
 * Extracts a particular bit (counted from the right) from an element
*/
template <typename T>
static inline bool extract_bit(T number, int offset_from_right)
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
//     Given a set of clusters, colours them.
//     The clusters must have an initial valid colouring stored in their temp_colour folder
//     The clusters should be for "vertices" s.t. one hop from each vertex is a cluster representing an edge
//     and two hops away is a vertex representing a neighbouring node
// */
template<typename T>
void colour_clusters(parlay::sequence<cluster<T>*> clusters)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;    

    parlay::parallel_for(0, clusters.size(), [&] (T v) {
        unsigned long local_maximum = clusters[v]->get_colour();
        unsigned long local_minimum = clusters[v]->get_colour();

        for(uint i = 0; i < clusters[v]->neighbours.size(); i++)
        {

            auto edge_ptr = clusters[v]->neighbours[i].load();
            
            if (edge_ptr == nullptr)
                continue;
            
            auto other_node_ptr = edge_ptr->get_other_side(clusters[v]);

            unsigned long compared_colour = other_node_ptr->get_colour();
            if(compared_colour > local_maximum)
                local_maximum = compared_colour;
            if(compared_colour < local_minimum)
                local_minimum = compared_colour;   
        }

        if(local_maximum == clusters[v]->get_colour()) // This node is a local maximum, give it a unique colour
        {
            clusters[v]->colour = local_maximum_colour;
        }
        else if(local_minimum == clusters[v]->get_colour())
        {
            clusters[v]->colour = local_minimum_colour;
        }
        else
        {
            clusters[v]->colour = 2 + (get_single_colour_contribution(clusters[v]->get_colour(), local_maximum) / 2); // bit shifting right to eliminate that pesky indicator bit
        }

    });

    return;

}

// /*
//     sets a boolean flag in the clusters indicating that they're part of MIS
//     also may change the boolean flag of some other clusters, only consider the clusters in this
//     These clusters must have a maximum degree of 2
// */
template<typename T>
void set_MIS(parlay::sequence<cluster<T>*> clusters)
{

    colour_clusters(clusters);

    parlay::parallel_for(0, clusters.size(), [&] (T v) {
        clusters[v]->is_MIS = false;
        clusters[v]->set_neighbour_mis(false);
    });

    auto colours = parlay::tabulate(clusters.size(), [&] (T v) {
       return clusters[v]->colour;
    });

    auto vertices = parlay::tabulate(clusters.size(), [&] (T v) {
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
            if(clusters[v]->get_neighbour_MIS() == true)
            {
                clusters[v]->is_MIS = false;
                return;
            }
            clusters[v]->is_MIS = true;
        });
    }

}


/**
 * Given an ASSYMETRIC graph, creates a set of clusters
 * In total, it creates n + m clusters in the array base_clusters
 * The first n are base_vertex clusters
 * And the last m are base_edge clusters
 * These clusters are linked to point to each other appropriately
 * it returns an array of n clusters corresponding to the original n vertices.
 * The edge clusters are present as edges.
*/
template <typename T>
void create_base_clusters(parlay::sequence<parlay::sequence<T>> &G, parlay::sequence<cluster<T> > &base_clusters, const T max_size)
{

    T n = G.size();

    base_clusters = parlay::tabulate(n, [&] (T v) {
        cluster<T> base_cluster;
        base_cluster.index = v;
        base_cluster.state = base_vertex | live;
        base_cluster.neighbours = parlay::sequence<std::atomic<cluster<T>*>>(max_size);
        base_cluster.children = parlay::sequence<std::atomic<cluster<T>*>>(2*max_size);
        parlay::parallel_for(0, max_size, [&] (T d){
            base_cluster.neighbours[d] = nullptr;
            base_cluster.children[d] = nullptr;
            base_cluster.counter[0] = 0;
        });
        return base_cluster;
    });


    // // populate base edge clusters
    parlay::parallel_for(0, n, [&] (T v) {    
        for(uint i = 0; i < G[v].size(); i++)
        {
            cluster<T>* edge_cluster = new cluster<T>;
            edge_cluster->index = -1;
            edge_cluster->state = base_edge | live;
            edge_cluster->neighbours = parlay::sequence<std::atomic<cluster<T>*>>(2);
            edge_cluster->neighbours[0] = &base_clusters[v];
            edge_cluster->neighbours[1] = &base_clusters[G[v][i]];
            edge_cluster->counter[0] = 0;
            base_clusters[v].neighbours[i] = edge_cluster;
            
        }
    });


    
    parlay::parallel_for(0, n, [&] (T v) {
        for(uint i = 0; i < G[v].size(); i++) // It is fine because constant degree graph so the lock isn't too bad
        {
            auto edge_cluster = base_clusters[v].neighbours[i].load();
            auto other_side = edge_cluster->neighbours[1].load();
            other_side->add_neighbour(edge_cluster);
            
        }
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
*/

template <typename T>
void create_RC_tree(parlay::sequence<cluster<T> > &base_clusters, T n)
{
    

    parlay::sequence<cluster<T>*> all_cluster_ptrs = parlay::tabulate(n, [&] (T v) {
        return &base_clusters[v];
    });

    parlay::sequence<cluster<T>*> forest, candidates;

    // Initially the forest of live nodes is all live nodes
    forest = parlay::filter(all_cluster_ptrs, [&] (cluster<T>* C) {
        return ((C->state & live));
    });

    do
    {
    // Shrink the forst as live nodes decrease
    forest = parlay::filter(forest, [&] (cluster<T>* C) {
        return ((C->state & live));
    });

    std::cout << "forest.size(): " << forest.size() << std::endl;

    // Eligible nodes are those with 0, 1 or 2 neighbours
    auto eligible = parlay::filter(forest, [&] (cluster<T>* C) {
        return (C->get_neighbour_count() <= 2);
    });

    std::cout << "eligible size is " << eligible.size() << std::endl;
    
    // Set the flag is_MIS amongst them
    set_MIS(eligible);

    // Filter out an MIS of eligible nodes
    candidates = parlay::filter(eligible, [&] (cluster<T>* C) {
        return (C->is_MIS);
    });

    std::cout << "candidates.size(): " << candidates.size() << std::endl;

    // do rake and compress
    parlay::parallel_for(0, candidates.size(), [&] (T v) {
        cluster<T>* cluster_ptr = candidates[v];
        // rake
        if(cluster_ptr->get_neighbour_count() == 0)
        {
            cluster_ptr->state&=(~live);
            cluster_ptr->state|=(nullary_cluster);
            cluster_ptr->state|=internal;
        }
        if(cluster_ptr->get_neighbour_count() == 1)
        {
            cluster<T>* edge_ptr = nullptr;
            for(T i = 0; i < cluster_ptr->neighbours.size(); i++)
            {
                if(cluster_ptr->neighbours[i] != nullptr)
                    edge_ptr = cluster_ptr->neighbours[i];
            }
            
            // find the other side of the edge
            cluster<T>* other_side = edge_ptr->get_other_side(cluster_ptr);
            
            // delete neighbours of edgeptr
            // edge_ptr->remove_neighbour(cluster_ptr);
            // edge_ptr->remove_neighbour(other_side);

            // delete edge_ptr as a neighbour of this
            cluster_ptr->remove_neighbour(edge_ptr);

            // delete edge_ptr as a neighbour for the other side
            other_side->remove_neighbour(edge_ptr);

            // now add these two as children

            other_side->add_child(edge_ptr);
            other_side->add_child(cluster_ptr);

            // make other_side be the parent for both
            edge_ptr->parent = cluster_ptr;
            cluster_ptr->parent = other_side;

            // mark both of these as not live
            edge_ptr->state&=(~live);
            cluster_ptr->state&=(~live);

            cluster_ptr->state|=unary_cluster;
            cluster_ptr->state|=internal;
        }
        else 
        if (cluster_ptr->get_neighbour_count() == 2)
        {
            // find left and right vertices/nodes
            cluster<T>* left_edge_ptr;
            cluster<T>* right_edge_ptr;

            cluster_ptr->get_two_neighbouring_edges(&left_edge_ptr, &right_edge_ptr);

            if(left_edge_ptr == NULL || right_edge_ptr == NULL)
            {
                std::cout << "This should never happen" << std::endl;
            }
            
            auto left_node_ptr = left_edge_ptr->get_other_side(cluster_ptr);
            auto right_node_ptr = right_edge_ptr->get_other_side(cluster_ptr);
            
            left_node_ptr->remove_neighbour(left_edge_ptr);
            left_node_ptr->add_neighbour(cluster_ptr);

            right_node_ptr->remove_neighbour(right_edge_ptr);
            right_node_ptr->add_neighbour(cluster_ptr);
            
            left_edge_ptr->parent = cluster_ptr;
            right_edge_ptr->parent = cluster_ptr;

            cluster_ptr->add_child(left_edge_ptr);
            cluster_ptr->add_child(right_edge_ptr);

            left_edge_ptr->state&=(~live);
            right_edge_ptr->state&=(~live);

            cluster_ptr->remove_neighbour(left_edge_ptr);
            cluster_ptr->add_neighbour(left_node_ptr);

            cluster_ptr->remove_neighbour(right_edge_ptr);
            cluster_ptr->add_neighbour(right_node_ptr);


            cluster_ptr->state&=(~live);

            cluster_ptr->state|=binary_cluster;
            cluster_ptr->state|=internal;
            
        }
    });
    }while(candidates.size());

}

/**
 * Make sure edges are undirected
*/

template <typename T, typename func>
void adjust_weights(parlay::sequence<cluster<T> > &clusters, parlay::sequence<std::tuple<T, T, T>> weighted_edges, func lambdafunc)
{

    parlay::sequence<cluster<T>*> relevant_edges = parlay::sequence<cluster<T>*>(weighted_edges.size(), nullptr);

    parlay::parallel_for(0, weighted_edges.size(), [&] (T i) {
        T v = std::get<0>(weighted_edges[i]);
        T w = std::get<1>(weighted_edges[i]);
        T weight = std::get<2>(weighted_edges[i]);

        if(v > w)
            std::swap(v, w);
    
        cluster<T>* v_cluster = &clusters[v];
        cluster<T>* w_cluster = &clusters[w];

        // either v or w will have the v-w edge as a child
        cluster<T>* v_w_edge = nullptr;

        for(T i = 0; i < v_cluster->children.size(); i++)
        {
            if(v_cluster->children[i].load() != nullptr)
            {
                if(v_cluster->children[i].load()->is_neighbour(v_cluster) && v_cluster->children[i].load()->is_neighbour(w_cluster))
                {
                    v_w_edge = v_cluster->children[i].load();
                    break;
                }
            }
        }
        if(v_w_edge == nullptr)
        {
            for(T i = 0; i < w_cluster->children.size(); i++)
            {
                if(w_cluster->children[i].load() != nullptr)
                {
                    if(w_cluster->children[i].load()->is_neighbour(v_cluster) && w_cluster->children[i].load()->is_neighbour(w_cluster))
                    {
                        v_w_edge = w_cluster->children[i].load();
                        break;
                    }
                }
            }
        }
        if(v_w_edge == nullptr)
        {
            std::cout << "This shouldn't happen, did you feed an invalid edge? " << v <<  " " << w <<std::endl;
            v_cluster->print();
            w_cluster->print();
            relevant_edges[i] = nullptr;
            return;
        }
        v_w_edge->data = weight;
        relevant_edges[i] = v_w_edge;
    });

    parlay::parallel_for(0, relevant_edges.size(), [&] (T v){
        auto node = relevant_edges[v];
        T ret_val;
        while(node->parent != nullptr && (ret_val = node->parent->counter[0].fetch_add(1)) == 0)
        {
            node = node->parent;
        }
    });

    std::cout << "counts: " << std::endl; 
    for(uint i = 0; i < clusters.size(); i++)
    {
        std::cout << clusters[i].index << ":" << clusters[i].counter[0].load() << std::endl;
    }
    std::cout << std::endl;

    parlay::parallel_for(0, relevant_edges.size(), [&] (T v){
        auto node = relevant_edges[v];
        T ret_val;
        while(node->parent != nullptr && (ret_val = node->parent->counter[0].fetch_add(-1)) == 1)
        {
            
            // All the children are ready
            

            node = node->parent;
        }
    });

    std::cout << "counts: " << std::endl; 
    for(uint i = 0; i < clusters.size(); i++)
    {
        std::cout << clusters[i].index << ":" << clusters[i].counter[0].load() << std::endl;
    }
    std::cout << std::endl;

    



}

/**
 * Since edges are dynamically allocated, delete them
 * takes the vertices! 
*/
template<typename T>
void delete_RC_Tree_edges(parlay::sequence<cluster<T> > &base_clusters)
{

    parlay::parallel_for(0, base_clusters.size(), [&] (T v) {
        for(T i = 0; i < base_clusters[v].neighbours.size(); i++)
        {
            
            if(base_clusters[v].neighbours[i].load() != nullptr && base_clusters[v].neighbours[i].load()->state & base_edge)
                base_clusters[v].neighbours[i].load()->counter[0].fetch_add(1);
        }
        for(T i = 0; i < base_clusters[v].children.size(); i++)
        {
            
            if(base_clusters[v].children[i].load() != nullptr && base_clusters[v].children[i].load()->state & base_edge)
                base_clusters[v].children[i].load()->counter[0].fetch_add(1);
        }
    });

    parlay::parallel_for(0, base_clusters.size(), [&] (T v) {
        for(T i = 0; i < base_clusters[v].neighbours.size(); i++)
        {
            
            if(base_clusters[v].neighbours[i].load() != nullptr && base_clusters[v].neighbours[i].load()->state & base_edge)
            {

                auto ret_count = base_clusters[v].neighbours[i].load()->counter[0].fetch_add(-1);
                if(ret_count == 1)
                {
                    delete base_clusters[v].neighbours[i].load();
                }
            }
        }
        for(T i = 0; i < base_clusters[v].children.size(); i++)
        {
            
            if(base_clusters[v].children[i].load() != nullptr && base_clusters[v].children[i].load()->state & base_edge)
            {

                auto ret_count = base_clusters[v].children[i].load()->counter[0].fetch_add(-1);
                if(ret_count == 1)
                {
                    delete base_clusters[v].children[i].load();
                }
            }
        }
    });
    
}

