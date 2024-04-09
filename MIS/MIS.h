#include <atomic>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"


/**
 * @brief Make sure the graph has a maximum of max_degree edges
 * 
 * @param[in] The graph
 * @param[in] The maximum degree as an integer
 * 
 * Takes the adjacency list of each node and removes the edges until it has a size of max_degree
 * It should work in-place
 * @return nothing
*/

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

template <typename graph> 
auto remove_higher_degree(graph& G, const int max_degree)
{   
    auto n = G.size();
    auto vertices = parlay::tabulate<int>(n, [&] (int i) {return n;});

    parlay::parallel_for(0, vertices.size(), [&] (int i) {
        if(G[i].size() > max_degree)
            G[i] = G[i].subseq(0, max_degree); // Does this result in gaps?
    });
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

    parlay::parallel_for(0, adjacent_colours.size(), [&] (int v) {
        final_string[v] = get_single_colour_contribution(mycolour, adjacent_colours[v]);
    });


    return final_string;
}

template <typename graph>
auto colour_graph(graph& G, const int max_degree)
{
    auto initial_colour = parlay::tabulate(G.size(),[&] (int i) {
        return get_vertex_colour_string(i, G[i], max_degree);        
    } );

    return initial_colour;
}

