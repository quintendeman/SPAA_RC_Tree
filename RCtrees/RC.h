#include <atomic>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/examples/helper/graph_utils.h"
#include <set>
#include <iostream>

template <typename T>
struct ColourIndexPair {
    public:
        T colour;
        T index;
        ColourIndexPair(T incolour, T inindex)
        {
            this->colour = incolour;
            this->index = inindex;
        }

        bool operator==(const ColourIndexPair& other) const {
            return this->colour == other->colour;
        }
        bool operator<(const ColourIndexPair& other) const {
            return this->colour < other->colour;
        }
        bool operator>(const ColourIndexPair& other) const {
            return this->colour > other->colour;
        }
        bool operator>=(const ColourIndexPair& other) const {
            return this->colour >= other->colour;
        }
        bool operator<=(const ColourIndexPair& other) const {
            return this->colour <= other->colour;
        }
        ColourIndexPair& operator=(const ColourIndexPair& other) {
        if (this != &other) { // Avoid self-assignment
            this->colour = other.colour;
            this->index = other.index;
        }
        return *this; // Return a reference to the modified object
        }

        operator long int() const {
            return this->colour; // Return the colour
        }


        
        friend std::ostream& operator<<(std::ostream& os, const ColourIndexPair& obj) {
        os << obj.colour << ":" << obj.index;
        return os;
        }

};


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

      if(max_degree > 2)
    {
        std::cout << "Error: Need a chain" << std::endl;  
        return colours;  
    }
 

    parlay::parallel_for(0, G.size(), [&] (T v) {
        auto edgelist = G[v]; // TODO replace with G[v] to remove extra writes
        T local_maximum = v;
        T local_minimum = v;

        for(uint i = 0; i < std::min(max_degree, (T) edgelist.size()); i++)
        {
            T w = edgelist[i];
            if(w > local_maximum)
            {
                local_maximum = w;
            }
            if (w < local_minimum)
            {
                local_minimum = w;
            }
        }

        if(local_maximum == v) // This node is a local maximum, give it a unique colour
        {
            colours[v] = local_maximum_colour;
        }
        else if(local_minimum == v)
        {
            colours[v] = local_minimum_colour;
        }
        else
        {
            colours[v] = 2 + (get_single_colour_contribution(v, local_maximum) / 2); // bit shifting right to eliminate that pesky indicator bit
        }

    });

    // print_graph(G, colours);

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

template <typename T>
parlay::sequence<ColourIndexPair<T>> colours_to_pairs( parlay::sequence<T> colours)
{
    parlay::sequence<ColourIndexPair<T>> ret_pairs = parlay::tabulate(colours.size(), [&] (T v) {
        return ColourIndexPair(colours[v], v);
    });
    
    return ret_pairs;
}