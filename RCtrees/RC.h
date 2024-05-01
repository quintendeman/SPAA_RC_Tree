#include <atomic>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/examples/helper/graph_utils.h"
#include <set>
#include <iostream>

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


template <typename T>
struct cluster
{
    public:
        T index = -1;
        short state = empty_type; // unaffected, affected or update eligible
        parlay::sequence<cluster*> data;
        
};

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


/*
    Only works on symmetric graphs with no redundancies
    Like the one returned from rmat_symmetric_graph
*/
template <typename vertex>
auto degree_cap_graph(parlay::sequence<parlay::sequence<vertex>>& G, const vertex max_degree)
{   
    vertex n = G.size();
    auto vertices = parlay::tabulate<vertex>(n, [&] (vertex i) {return n;});

    parlay::parallel_for(0, vertices.size(), [&] (vertex i) {
        if(G[i].size() > max_degree)
            G[i] = G[i].subseq(0, max_degree); // Does this result in gaps?
    });

    delete_assymetric_pairs(G);

    return;
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
    Works on a CHAIN GRAPH
*/
template <typename T, typename graph>
parlay::sequence<T> colour_chains_to_logn(graph G, const T max_degree = 2)
{
    static const T local_maximum_colour = (T) 0;
    static const T local_minimum_colour = (T) 1;
    

    parlay::sequence<T> colours = parlay::tabulate(G.size(), [&]  (T v) {
        return v;
    });

    auto initial_colours = colours;

      if(max_degree > 2)
    {
        std::cout << "Error: Need a chain" << std::endl;  
        return colours;  
    }
 

    parlay::parallel_for(0, G.size(), [&] (T v) {
        // auto edgelist = G[v]; // TODO replace with G[v] to remove extra writes
        T local_maximum = initial_colours[v];
        T local_minimum = initial_colours[v];

        for(uint i = 0; i < std::min(max_degree, (T) G[v].size()); i++)
        {
            // T w = G[v][i];
            if(initial_colours[G[v][i]] > local_maximum)
            {
                local_maximum = initial_colours[G[v][i]];
            }
            if (initial_colours[G[v][i]] < local_minimum)
            {
                local_minimum = initial_colours[G[v][i]];
            }
        }

        if(local_maximum == initial_colours[v]) // This node is a local maximum, give it a unique colour
        {
            colours[v] = local_maximum_colour;
        }
        else if(local_minimum == initial_colours[v])
        {
            colours[v] = local_minimum_colour;
        }
        else
        {
            colours[v] = 2 + (get_single_colour_contribution(initial_colours[v], local_maximum) / 2); // bit shifting right to eliminate that pesky indicator bit
        }

    });

    return colours;
} 



template <typename T, typename graph>
bool is_valid_colouring(graph G, parlay::sequence<T> colours)
{

    parlay::sequence<T> wrongs = parlay::tabulate(G.size(), [&]  (T v) {
        return (T) -1;
    });

    bool is_valid = true;


    parlay::parallel_for(0, G.size(), [&] (T v) {
        auto edgelist = G[v];
        parlay::parallel_for(0, edgelist.size(), [&] (T ind) {
            auto w = edgelist[ind];
            if(colours[v] == colours[w])
            {
                is_valid = false;
                wrongs[v] = w;
            }
        });
    });

    for(uint v = 0; v < G.size(); v++)
    {
        if (wrongs[v] != -1)
        {
            std::cout << "Vertex " << v <<" has a wrong colour: " << colours[v] << std::endl;
            std::cout << v << "\'s neighbour colours are: ";
            for(uint i = 0; i < G[v].size(); i++)
                std::cout << "(" << G[v][i] << ")" << ":" << colours[G[v][i]] << " ";
            std::cout << std::endl;  
            std::cout << "Conflicting node's ID: " << wrongs[v] << std::endl;
            std::cout << "Conflicting node's neighbour colours: ";
            for(uint i = 0; i < G[wrongs[v]].size(); i++)
                std::cout << "(" << G[wrongs[v]][i] << ")"  << ":" << colours[G[wrongs[v]][i]] << " ";
            std::cout << std::endl << std::endl;
        }
    }
    return is_valid;

}


template <typename T, typename graph>
parlay::sequence<bool> get_MIS(graph G, parlay::sequence<T> sorted_vertices, parlay::sequence<unsigned long> offsets)
{
    auto n = G.size();

    parlay::sequence<bool> is_in_MIS = parlay::tabulate(n, [&] (T v) {
        return false;
    });

    // iterate over the colours
    for(T i = 0; i < offsets.size(); i++)
    {
        T start_index;
        if (i == 0)
            start_index = 0;
        else
            start_index = offsets[i-1];
        T end_index = offsets[i];

        parlay::parallel_for(start_index, end_index, [&] (T v) {
            bool keep_this_node = true;

            for(T w = 0; w < G[v].size(); w++)
            {
                if(is_in_MIS[G[v][w]])
                    keep_this_node = false;
            }

            is_in_MIS[v] = keep_this_node;
        });
    }
    return is_in_MIS;
}

/*
    Generate a simple, single rooted graph with each node having two children
    Then, randomly, change the parent of each node with a certain probability such that it picks something on the left of it
*/
template <typename T>
parlay::sequence<T> generate_tree_graph(T num_elements, bool randomized = true, bool sequential = false)
{
    assert(num_elements > 0);

    parlay::sequence<T> dummy_initial_parents = parlay::tabulate(num_elements, [&] (T v) {return (T) 0;});

    if(sequential)
    {
       parlay::parallel_for(0, num_elements, [&] (T i) {
            
            if(i)
                dummy_initial_parents[i] = i-1;
            else 
                dummy_initial_parents[0] = 0; // root's parent is itself

        }); 
        return dummy_initial_parents;
    }

    parlay::parallel_for(0, num_elements, [&] (T i) {
        dummy_initial_parents[i] = i/2;
    }); 
    if(!randomized)
        return dummy_initial_parents;

    parlay::sequence<T> vertices = parlay::tabulate(num_elements, [&] (T v) {return v;});
    parlay::sequence<T> random_index = parlay::random_shuffle(vertices);

    parlay::parallel_for(0, num_elements, [&] (T v) {
        T picked_parent = random_index[v];
        if(picked_parent < v)
            dummy_initial_parents[v] = picked_parent;
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
        parlay::sequence<T> temp = parlay::tabulate(1, [&] (T i) {return i;});
        temp[0] = parents[v];
        return temp;
    });

    G = graph_utils<T>::symmetrize(G);

    return G;
}

template <typename T>
parlay::sequence<cluster<T>> create_RC_Tree(parlay::sequence<parlay::sequence<T>> G, const T max_degree)
{

    // create base clusters with no edges
    parlay::sequence<cluster<T>> base_clusters = parlay::tabulate(G.size(), [&] (T v) {
        cluster<T> base_cluster;
        base_cluster.index = v;
        base_cluster.state = base_vertex | live;
        return base_cluster;
    });

    // add edges
    // If a base_cluster directly pointers to another base_cluster, this represents an edge
    parlay::parallel_for(0, base_clusters.size(), [&] (T v) {
        base_clusters[v].data = parlay::tabulate(G[v].size(), [&] (T w) {
            return &base_clusters[G[v][w]];
        });
    });

    parlay::sequence<cluster<T>> internal_clusters = parlay::tabulate(G.size(), [&] (T v) {
        cluster<T> base_cluster;
        base_cluster.index = v;
        base_cluster.state = empty_type;
        return base_cluster;
    });

    return base_clusters;

}


