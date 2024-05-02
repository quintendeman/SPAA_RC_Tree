#include <atomic>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/examples/helper/graph_utils.h"
#include "/scratch/parlaylib/examples/counting_sort.h"
#include <random>
#include <set>
#include <iostream>
#include <mutex>

const short empty_type = 0;
const short base_vertex = 1;
const short base_edge = 2;
const short unary_cluster = 4;
const short binary_cluster = 8;
const short nullary_cluster = 16;
const short needs_colouring = 32;
const short unaffected = 64;
const short affected = 128;
const short live = 256;
const short eligible = 512;
const short uninitialized = 2048;
const short initialized = 4096;
const short internal = 8192;


template <typename T>
struct cluster
{
    public:
        T index = -1;
        short state = empty_type; // unaffected, affected or update eligible
        parlay::sequence<cluster<T>*> data;
        cluster<T>* parent;
        T temp_colour;
        T final_colour = -1;
        T is_MIS = false;
        
};

/*
    Generate a simple, single rooted graph with each node having two children
    Then, randomly, change the parent of each node with a certain probability such that it picks something on the left of it
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
    Converts parents into a graph
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
 * 
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




/*
    Basic workflow
    1) Find edges that are assymetric (do a parfor on the edge list of the target array)
    2) Mark them in a global, boolean graph
    then
    3) filter them out
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


template<typename T, typename graph>
void print_graph(graph G, parlay::sequence<T> colours)
{
    for(uint i = 0; i < G.size(); i++)
    {
        std::cout << i << ":" << colours[i] << std::endl << "   ";
        for(uint j = 0; j < G[i].size(); j++)
        {
            std::cout << G[i][j] << ":" << colours[G[i][j]] << " ";
        }
        std::cout << std::endl;
    }


    return;
}


/*
    WARNING, ONLY WORKS ON CHAINS WITH 2 EDGES
*/
template<typename T>
void colour_clusters(parlay::sequence<cluster<T>*> clusters)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;    

    parlay::parallel_for(0, clusters.size(), [&] (T v) {
        T local_maximum = clusters[v]->temp_colour;
        T local_minimum = clusters[v]->temp_colour;

        for(uint i = 0; i < clusters[v]->data.size(); i++)
        {
            auto edge_ptr = clusters[v]->data[i];
            auto other_node_ptr = edge_ptr->data[0];
            if(other_node_ptr == clusters[v])
                other_node_ptr = edge_ptr->data[1];

            T compared_colour = other_node_ptr->temp_colour;
            if(compared_colour > local_maximum)
                local_maximum = compared_colour;
            if(compared_colour < local_minimum)
                local_minimum = compared_colour;   
        }

        if(local_maximum == clusters[v]->temp_colour) // This node is a local maximum, give it a unique colour
        {
            clusters[v]->final_colour = local_maximum_colour;
        }
        else if(local_minimum == clusters[v]->temp_colour)
        {
            clusters[v]->final_colour = local_minimum_colour;
        }
        else
        {
            clusters[v]->final_colour = 2 + (get_single_colour_contribution(clusters[v]->temp_colour, local_maximum) / 2); // bit shifting right to eliminate that pesky indicator bit
        }

    });

    return;

}

/*
    sets a boolean flag in the clusters indicating that they're part of MIS
    also may change the boolean flag of some other clusters, only consider the clusters in this
    These clusters must have a maximum degree of 2
*/
template<typename T>
void set_MIS(parlay::sequence<cluster<T>*> clusters)
{

    colour_clusters(clusters);

    parlay::parallel_for(0, clusters.size(), [&] (T v) {
        clusters[v]->is_MIS = false;
        for(uint i = 0; i < clusters[v]->data.size(); i++)
        {
            clusters[v]->data[i]->data[0]->is_MIS = false;
            clusters[v]->data[i]->data[1]->is_MIS = false;
        }
    });

    auto colours = parlay::tabulate(clusters.size(), [&] (T v) {
       return clusters[v]->final_colour;
    });

    auto vertices = parlay::tabulate(clusters.size(), [&] (T v) {
        return v;
    });

    auto result = vertices;

    parlay::sequence<unsigned long> offsets = counting_sort(vertices.begin(), vertices.end(), result.begin(), colours.begin(), 8 * sizeof(T));

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
            bool keep_this_node = true;
            for(uint w = 0; w < clusters[v]->data.size(); w++)
            {
                if(clusters[v]->data[w]->data[0]->is_MIS == true || clusters[v]->data[w]->data[1]->is_MIS == true)
                {
                    keep_this_node = false;
                    break;
                }
            }
            clusters[v]->is_MIS = keep_this_node;
        });
    }

}


/*
    sets a boolean flag in the clusters indicating that they're part of MIS
    also may change the boolean flag of some other clusters, only consider the clusters in this
    These clusters must have a maximum degree of 2
*/
template<typename T>
void check_MIS(parlay::sequence<cluster<T>*> clusters)
{

    bool is_valid_MIS = true;

    T example_index;
    cluster<T>* W;

    for(T v = 0; v < clusters.size(); v++)
    {
        if(clusters[v]->is_MIS)
        {
            // check if any neighbours are valid MIS
            for(uint i = 0; i < clusters[v]->data.size();i++)
            {
                auto edge_ptr = clusters[v]->data[i];
                auto other_node = edge_ptr->data[0];
                if (other_node->index == clusters[v]->index)
                {
                    other_node = edge_ptr->data[1];
                }
                
                if(other_node->is_MIS && other_node->data.size()<=2)
                {
                    is_valid_MIS = false;
                    example_index = v;
                    W = other_node;
                }
            }
        }
    }

    if(!is_valid_MIS)
    {
        std::cout << "Not a valid MIS: " << std::endl;
        auto V = clusters[example_index];
        std::cout << V->index << " " << V->final_colour << " ";
        for(uint i = 0; i < V->data.size(); i++)
        {
            auto other_ptr = V->data[i]->data[0];
            if(other_ptr == V)
            {
                other_ptr = V->data[i]->data[1];
            }
            std::cout << other_ptr->index << ":" << other_ptr->final_colour << " ";
        }
        std::cout << std::endl;
        std::cout << W->index << " " << W->final_colour << " ";
        for(uint i = 0; i < W->data.size(); i++)
        {
            auto other_ptr = W->data[i]->data[0];
            if(other_ptr == W)
            {
                other_ptr = W->data[i]->data[1];
            }
            std::cout << other_ptr->index << ":" << other_ptr->final_colour << " ";
        }
        std::cout << std::endl;
    }

}


/**
* Input: G, an ASSYMETRIC graph
*/
template <typename T>
void create_base_clusters(parlay::sequence<parlay::sequence<T>> &G, parlay::sequence<cluster<T> > &base_clusters)
{

    T n = G.size();
    auto [sums, m] = parlay::scan(parlay::tabulate(G.size(), [&] (T v) {
        return G[v].size();
    }));

    base_clusters = parlay::tabulate(n+m, [&] (T v) {
        cluster<T> base_cluster;
        base_cluster.index = v;
        base_cluster.temp_colour = v;
        if(v < n)
            base_cluster.state = base_vertex | live;
        else
            base_cluster.state = base_edge | live;
        return base_cluster;
    });


    // populate base edge clusters
    parlay::parallel_for(0, n, [&] (T v) {
        
        for(uint i = 0; i < G[v].size(); i++)
        {
            auto location = n + sums[v] + i;
            base_clusters[location].data.push_back(&base_clusters[v]);
            base_clusters[location].data.push_back(&base_clusters[G[v][i]]);
        }
    });


    // connect outgoing edges
    parlay::parallel_for(0, n, [&] (T v) {
        for(uint i = 0; i < G[v].size(); i++) // It is fine because constant degree graph
        {
            auto edge_location = n + sums[v] + i;
            base_clusters[v].data.push_back(&base_clusters[edge_location]);
        }
    });


    std::mutex* mutexes = new std::mutex[n];

    // connect incoming edges
    parlay::parallel_for(0, n, [&] (T v) {
        for(uint i = 0; i < G[v].size(); i++) // It is fine because constant degree graph
        {
            auto edge_location = n + sums[v] + i;
            auto w = G[v][i]; // destination
            std::unique_lock<std::mutex> lock(mutexes[w]);
            base_clusters[w].data.push_back(&base_clusters[edge_location]);
            lock.unlock();
        }
    });

    delete[] mutexes;

}

template <typename T>
void create_RC_tree(parlay::sequence<cluster<T> > &base_clusters, T n)
{
    

    parlay::sequence<cluster<T>*> all_cluster_ptrs = parlay::tabulate(n, [&] (T v) {
        return &base_clusters[v];
    });

    parlay::sequence<cluster<T>*> candidates;

    uint tester = 0;

    do
    {
    
    // TODO: Convert this to approximate compaction!!!
    candidates = parlay::filter(all_cluster_ptrs, [&] (cluster<T>* C) {
        return (C->data.size() <= 2 && (C->state & live));
    });

    set_MIS(candidates);

    candidates = parlay::filter(all_cluster_ptrs, [&] (cluster<T>* C) {
        return (C->is_MIS);
    });

    // check_MIS(candidates);

    // do rake and compress
    parlay::parallel_for(0, candidates.size(), [&] (T v) {
        cluster<T>* cluster_ptr = candidates[v];
        if(cluster_ptr->is_MIS == false)
            return;
        // rake
        if(cluster_ptr->data.size() == 0)
        {
            cluster_ptr->state&=(~live);
            cluster_ptr->state|=(nullary_cluster);
        }
        if(cluster_ptr->data.size() == 1)
        {
            cluster<T>* edge_ptr = cluster_ptr->data[0];
            
            // find the other side of the edge
            cluster<T>* other_side = edge_ptr->data[0];
            if(other_side == cluster_ptr)
                other_side = edge_ptr->data[1];
            
            if(other_side == NULL)
                std::cout << "This should never happen" << std::endl;
            
            // delete this edge in the other_side
            for(uint i = 0; i < other_side->data.size(); i++)
            {
                if(other_side->data[i] == edge_ptr)
                {
                    other_side->data[i] = NULL;
                    break;
                }
            }

            // make other_side be the parent for both
            // TODO: replace with some other operation that could be performed here

            edge_ptr->parent = cluster_ptr;
            cluster_ptr->parent = other_side;

            // mark both of these as not live

            edge_ptr->state&=(~live);
            cluster_ptr->state&=(~live);

            cluster_ptr->state|=unary_cluster;
        }
        else if (cluster_ptr->data.size() == 2)
        {
            auto left_edge_ptr = cluster_ptr->data[0];
            auto right_edge_ptr = cluster_ptr->data[1];

            auto left_node_ptr = left_edge_ptr->data[0];
            if(left_node_ptr == cluster_ptr)
                left_node_ptr =  left_edge_ptr->data[1];

            auto right_node_ptr = right_edge_ptr->data[0];
            if(right_node_ptr == cluster_ptr)
                right_node_ptr = right_edge_ptr->data[1];
            
            // realign left_edge_ptr in left_node_ptr
            for(uint i = 0; i < left_node_ptr->data.size(); i++)
            {
                if(left_node_ptr->data[i] == left_edge_ptr)
                {
                    left_node_ptr->data[i] = cluster_ptr;
                     break;
                }
            }
            // realign right_edge_ptr in right_node_ptr
            for(uint i = 0; i < right_node_ptr->data.size(); i++)
            {
                if(right_node_ptr->data[i] == right_edge_ptr)
                {
                    right_node_ptr->data[i] = cluster_ptr;
                     break;
                }
            }

            left_edge_ptr->parent = cluster_ptr;
            right_edge_ptr->parent = cluster_ptr;

            left_edge_ptr->state&=(~live);
            right_edge_ptr->state&=(~live);

            cluster_ptr->data[0] = left_node_ptr;
            cluster_ptr->data[1] = right_node_ptr;

            cluster_ptr->state&=(~live);

            cluster_ptr->state|=binary_cluster;

        }
    });

    // remove any extra NULLs
    parlay::parallel_for(0, n, [&] (T v) {
        base_clusters[v].data = parlay::filter(base_clusters[v].data, [&] (cluster<T>* C) {
            return C!=NULL;
        });
    });

    // find MIS amongst these
    tester+=1;
    }while(tester < 3);
    

}