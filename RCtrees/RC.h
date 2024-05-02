#include <atomic>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/examples/helper/graph_utils.h"
#include <random>
#include <set>
#include <iostream>
#include <mutex>

static const short empty_type = 0;
static const short base_vertex = 1;
static const short base_edge = 2;
static const short unary_cluster = 4;
static const short binary_cluster = 8;
static const short nullary_cluster = 16;
static const short needs_colouring = 32;
static const short unaffected = 64;
static const short affected = 128;
static const short live = 256;
static const short eligible = 512;
static const short dead = 1024;
static const short uninitialized = 2048;
static const short initialized = 4096;


template <typename T>
struct cluster
{
    public:
        T index = -1;
        short state = empty_type; // unaffected, affected or update eligible
        parlay::sequence<cluster<T>*> data;
        cluster<T>* parent;
        T temp_colour;
        T final_colour;
        T is_MIS;
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
void colour_clusters(parlay::sequence<cluster<T>>& clusters)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;    

    parlay::parallel_for(0, clusters.size(), [&] (T v) {
        T local_maximum = clusters[v].temp_colour;
        T local_minimum = clusters[v].temp_colour;

        for(uint i = 0; i < clusters[v].data.size(); i++)
        {
            T compared_colour = clusters[v].data[i]->temp_colour;
            if(compared_colour > local_maximum)
                local_maximum = compared_colour;
            if(compared_colour < local_minimum)
                compared_colour = local_minimum;   
        }

        if(local_maximum == clusters[v].temp_colour) // This node is a local maximum, give it a unique colour
        {
            clusters[v].final_colour = local_maximum_colour;
        }
        else if(local_minimum == clusters[v].temp_colour)
        {
            clusters[v].final_colour = local_minimum_colour;
        }
        else
        {
            clusters[v].final_colour = 2 + (get_single_colour_contribution(clusters[v].temp_colour, local_maximum) / 2); // bit shifting right to eliminate that pesky indicator bit
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
void set_MIS(parlay::sequence<cluster<T>>& clusters)
{

    colour_clusters(clusters);

    parlay::parallel_for(0, clusters.size(), [&] (T v) {
        clusters[v].is_MIS = false;
        for(uint i = 0; i < clusters[v].data.size(); i++)
            clusters[v].data[i]->is_MIS = false;
    });

    auto colours = parlay::tabulate(clusters.size(), [&] (T v) {
        clusters[v].final_colour;
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


        parlay::parallel_for(start_index, end_index, [&] (T v) {
            bool keep_this_node = true;
            for(uint w = 0; w < clusters[v].data.size(); w++)
            {
                if(clusters[v].data[w]->is_MIS == true)
                {
                    keep_this_node = false;
                    break;
                }
            }
            clusters[v].is_MIS = keep_this_node;
        });
    }


}

/**
* Input: G, an ASSYMETRIC graph
*/
template <typename T>
parlay::sequence<cluster<T>> create_RC_Tree(parlay::sequence<parlay::sequence<T>> &G, const T max_degree)
{

    T n = G.size();
    auto [sums, m] = parlay::scan(parlay::tabulate(G.size(), [&] (T v) {
        return G[v].size();
    }));

    parlay::sequence<cluster<T>> base_clusters = parlay::tabulate(n+m, [&] (T v) {
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






    return base_clusters;

}


