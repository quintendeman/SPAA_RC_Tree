#include <atomic>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"




template <typename T>
void printBits(T num) {
    int numBits = sizeof(T) * 8; // Total number of bits in an integer

    // Start from the most significant bit and iterate through each bit
    for (int i = numBits - 1; i >= 0; i--) {
        // Check if the bit at position i is set
        int mask = 1 << i;
        if (num & mask) {
            std::cout << "1";
        } else {
            std::cout << "0";
        }

        // Print a space after every 4 bits for better readability
        if (i % 4 == 0) {
            std::cout << " ";
        }
    }
}

void print_string(std::string input_string)
{
    for(uint i = 0; i < input_string.length(); i++)
    {
        std::cout << (int) input_string[i] << " ";
    }
    std::cout << std::endl;
}

/*
    Basic workflow
    1) Find edges that are assymetric (do a parfor on the edge list of the target array)
    2) Mark them in a global, boolean graph
    then
    3) filter them out
*/
template <typename graph>
void delete_assymetric_pairs(graph& G)
{
     auto vertices = parlay::tabulate<int>(G.size(), [&] (int i) {return i;});

     parlay::sequence<parlay::sequence<bool> >  keep_edges_graph = parlay::map(vertices, [&] (int v) {
        auto edge = G[v];
        auto keep_edge = parlay::tabulate<bool>(edge.size(), [&] (int i) {return (bool) false;});
        
        int starting_node = v;
        parlay::parallel_for(0, edge.size(), [&] (int e){
            int ending_node = edge[e];
            parlay::sequence<int> ending_edge_list = G[ending_node];
            
            // is starting node in ending node?
            parlay::parallel_for(0, ending_edge_list.size(), [&] (int w) {
                if (ending_edge_list[w] == starting_node)
                    keep_edge[e] = true;
            });
        });
        return keep_edge;
     });

     parlay::parallel_for(0, G.size(), [&] (int v) {
        auto edge = G[v];
        parlay::parallel_for(0, edge.size(), [&] (int w){
            edge[w] = keep_edges_graph[v][w] ? edge[w] : -1;
        }
        );

        edge = parlay::filter(edge, [&] (int w){
            return w != -1;
        });

        G[v] = edge;
     });

     
}


/*
    Only works on symmetric graphs with no redundancies
    Like the one returned from rmat_symmetric_graph
*/
template <typename graph>
auto return_degree_capped_graph(graph& G, const int max_degree)
{   
    int n = G.size();
    auto vertices = parlay::tabulate<int>(n, [&] (int i) {return n;});

    parlay::parallel_for(0, vertices.size(), [&] (int i) {
        if(G[i].size() > max_degree)
            G[i] = G[i].subseq(0, max_degree); // Does this result in gaps?
    });

    delete_assymetric_pairs(G);

    return;
}


/*
    Returns first different bit from the left

    sizeof(T) must be less than 16 bytes
*/
template <typename T>
inline char first_different_bit(const T a, const T b, bool* bit)
{
    T difference = a ^ b;
    char num_bits = sizeof(T) * 8;
    
    for(char i = num_bits - 1; i >= 0; i--)
    {
        T temp = difference >> i;
        if (temp)
        {   
            *bit = (b >> i) & 1;
            return num_bits - 1 - i;
        }
    }

    return -1;
}

/*
    returns a char with I_w and C_w(I_w) packed
    I_w is the index of the first different bit
    C_w is the value in neighbour w of this bit

    Also, technically I waste the left-most bit in each char
*/

static char get_single_colour_contribution(const int vcolour, const int wcolour)
{
    bool wbit = false;
    char different_bit = first_different_bit(vcolour, wcolour, &wbit);
    char final_returned_character = (different_bit << 1) | wbit;
    return final_returned_character;
}

std::string get_vertex_colour_string(int mycolour, parlay::sequence<int> adjacent_colours, const int max_degree)
{
    // calculate the size of one vertex_colour
    // This is equal to max_degree * ( log(max colour variety) + 1) 
    std::string final_string(max_degree, '\0');

    if(adjacent_colours.size() > max_degree )
        std::cout << "AHA!" << std::endl;

    parlay::parallel_for(0, adjacent_colours.size(), [&] (int v) {
        final_string[v] = get_single_colour_contribution(mycolour, adjacent_colours[v]);
    });


    return final_string;
}

template <typename graph>
auto colour_graph(graph& G, const int max_degree)
{
    auto initial_colour = parlay::tabulate(G.size(),[&] (int i) {
        return std::make_pair(get_vertex_colour_string(i, G[i], max_degree), i);        
    } );

    return initial_colour;
}


template<typename graph>
parlay::sequence<int> generate_MIS(graph& G, parlay::sequence<std::pair<std::string, int>> &string_values)
{
    parlay::sequence<int> mis_vertices = parlay::tabulate(G.size(), [&] (int v){return v;});
    parlay::sequence<bool> deleted = parlay::tabulate(G.size(), [&] (int v){return false;});
    
    int colour_index = 0;
    std::string top_colour = string_values[string_values.size()-1].first;
    
    while(string_values.size() > 0)
    {
        while(string_values.size() > 0 && string_values[string_values.size()-1].first == top_colour)
        {
            string_values.pop_back();
        }
        top_colour = string_values[string_values.size()-1].first;

        //TODO don't parfor for all the nodes, just the ones we know match the colour
        parlay::parallel_for(0, G.size(), [&] (int v){
            if(string_values[v].first != top_colour)
                return;
            int node = string_values[v].second;
            // is node deleted?
            if(deleted[node])
                return;
            // delete the neighbours of node
            auto edge_list = G[node];
            parlay::parallel_for(0, edge_list.size(), [&] (int cnt) {
                deleted[edge_list[cnt]] = true;
            });
        });

    }

    mis_vertices = parlay::filter(mis_vertices, [&] (int v) {
        return !deleted[v];
    });

    return mis_vertices;
}
