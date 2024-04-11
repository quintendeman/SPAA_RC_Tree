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
static unsigned char get_single_colour_contribution(const T vcolour, const T wcolour, char* different_bit_index)
{
    bool wbit = false;
    char different_bit = first_different_bit(vcolour, wcolour, &wbit);
    char final_returned_character = (different_bit << 1) | wbit;
    if(different_bit_index) *different_bit_index = different_bit;
    return final_returned_character;
}


/**
 * @in[G] A sequence of integers representing the parent of each node
 * 
 * Note that each node has a parent which is the ID of another node.
 * If it is a root node, the parent ID is itself
 * 
 * @returns A sequence of 8-bit chars representing the colour of each node
*/
template <typename T>
parlay::sequence<uint8_t> six_colour_rooted_tree(parlay::sequence<int> parents, parlay::sequence<T> initial_colours)
{

    parlay::sequence<uint8_t> colouring = parlay::tabulate(parents.size(), [&] (T v) {return (uint8_t) 0;});

    parlay::parallel_for(0, parents.size(),[&] (int v) {
        int parent_id = parents[v];
        if (parent_id == v) // root node
            colouring[v] = (uint8_t) 0;
        else
        {
            uint8_t my_colouring = get_single_colour_contribution(initial_colours[v], initial_colours[parent_id], NULL);

            colouring[v] =  my_colouring;
        }
    });
    return colouring;
}


/*
    Generate a simple, single rooted graph with each node having two children
    Then, randomly, change the parent of each node with a certain probability such that it picks something on the left of it
*/
parlay::sequence<int> generate_tree_graph(int num_elements, bool randomized = true, bool sequential = false)
{
    assert(num_elements > 0);

    parlay::sequence<int> dummy_initial_parents = parlay::tabulate(num_elements, [&] (int v) {return 0;});

    if(sequential)
    {
       parlay::parallel_for(0, num_elements, [&] (int i) {
            
            dummy_initial_parents[i] = i-1;

        }); 
        return dummy_initial_parents;
    }

    parlay::parallel_for(0, num_elements, [&] (int i) {
        dummy_initial_parents[i] = i/2;
    }); 
    if(!randomized)
        return dummy_initial_parents;

    parlay::sequence<int> vertices = parlay::tabulate(num_elements, [&] (int v) {return v;});
    parlay::sequence<int> random_index = parlay::random_shuffle(vertices);

    parlay::parallel_for(0, num_elements, [&] (int v) {
        int picked_parent = random_index[v];
        if(picked_parent < v)
            dummy_initial_parents[v] = picked_parent;
    });

    
    return dummy_initial_parents;
}