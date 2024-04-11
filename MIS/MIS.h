#include <atomic>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/examples/helper/graph_utils.h"




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
template <typename graph, typename vertex>
void delete_assymetric_pairs(graph& G)
{
     auto vertices = parlay::tabulate<vertex>(G.size(), [&] (vertex i) {return i;});

     parlay::sequence<parlay::sequence<bool> >  keep_edges_graph = parlay::map(vertices, [&] (vertex v) {
        auto edge = G[v];
        auto keep_edge = parlay::tabulate<bool>(edge.size(), [&] (vertex i) {return (bool) false;});
        
        int starting_node = v;
        parlay::parallel_for(0, edge.size(), [&] (vertex e){
            int ending_node = edge[e];
            parlay::sequence<int> ending_edge_list = G[ending_node];
            
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
template <typename graph, typename vertex>
auto return_degree_capped_graph(graph& G, const int max_degree)
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
auto six_colour_rooted_tree(parlay::sequence<T> parents, parlay::sequence<T> initial_colours)
{

    parlay::sequence<T> colouring = parlay::tabulate(parents.size(), [&] (T v) {return (T) 0;});

    auto same_count = 0;
    double L = (double) (sizeof(long) * 8);
    do
    {
        parlay::parallel_for(0, parents.size(),[&] (int v) {
        long parent_id = parents[v];
        if (parent_id == v) // root node
            colouring[v] = (uint8_t) 0;
        else
        {
            uint8_t my_colouring = get_single_colour_contribution(initial_colours[v], initial_colours[parent_id], NULL);

            colouring[v] =  my_colouring;
        }
        });
        initial_colours = colouring;
        if(L == ceil(log2(L) + 1))
            same_count++;
        L = ceil(log2(L) + 1);
        
    }while(same_count != 3);
    
    return colouring;
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
            
            dummy_initial_parents[i] = i-1;

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

template<typename T>
bool isBadColour(T colour, parlay::sequence<T> unique_colours)
{
    // assert(unique_colours.size() == 6);
    if (colour == unique_colours[0])
        return false;
    if (colour == unique_colours[1])
        return false;
    if (colour == unique_colours[2])
        return false;
    if (colour == unique_colours[3])
        return true;
    if (colour == unique_colours[4])
        return true;
    if (colour == unique_colours[5])
        return true;
    return false;
}

/*
    Convert 6 colour tree into 3 colour tree
*/
template <typename graph, typename T>
parlay::sequence<T> convert_6_to_3_tree(graph G, parlay::sequence<T> parents, parlay::sequence<T> original_colours)
{   



    static const T uncoloured = (T) -1;

    parlay::sequence<T> colours = parlay::tabulate(parents.size(), [&] (T v) {
        return original_colours[v];
    });

    // sample_sort(colours);  // TODO find a better way to get the 6 unique colours
    auto unique_colours = parlay::tabulate(6, [&] (T v) {return v;});



    // each node is recolored as its parent
    parlay::parallel_for(0, parents.size(), [&] (T v){
        if(parents[v] == v)
        {
            if (original_colours[v] == unique_colours[0]) // root recolours as any other colour
            {
                colours[v] = unique_colours[1];
            }
            else
                colours[v] = unique_colours[0];
            return;
        }
        colours[v] = original_colours[parents[v]];
    });

    // each node that is a bad colour and has a neighbour with a good colour
    // gets no colour
    parlay::parallel_for(0, parents.size(), [&] (T v){

        T my_original_colour = colours[v];
        bool current_is_bad = isBadColour(my_original_colour, unique_colours);
        if(!current_is_bad)
        {
            return;
        }

        auto neighbours = G[v];
        bool good_neighbour_exists = false;
        parlay::parallel_for(0, neighbours.size(),[&] (int e) {
            T neighbouring_node = neighbours[e];
            T neighbouring_node_colour = colours[neighbouring_node];
            if(!isBadColour(neighbouring_node_colour, unique_colours))
                good_neighbour_exists = true;
        });
        T parent_colour = colours[parents[v]];
        if(!isBadColour(parent_colour,unique_colours))
            colours[v] = uncoloured;

        if(good_neighbour_exists)
            colours[v] = uncoloured;

    });

    // recolour bad colours to colour mod 3
    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(colours[v] == uncoloured)
            return;
        if(isBadColour(colours[v], unique_colours))// 4 5 6 are the bad colours
            colours[v] = colours[v] % 3;
    });


    // recolour uncoloured nodes
    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(colours[v] != uncoloured)
            return;
        
        // pick the parent colour
        // This should work for root too
        T parent_colour = colours[parents[v]];

        // pick a child colour, all children should have the same colour

        T child_colour;

        if(G[v].size() > 0) // a child exists
        {
            child_colour = colours[G[v][0]]; //first child's colour
        }
        else
            child_colour = (parent_colour + 1) % 3;

        T my_colour = uncoloured;
        
        for(uint i = 0; i < 3; i++)
        {
           if(unique_colours[i] != child_colour && unique_colours[i] != parent_colour)
            {
                my_colour = unique_colours[i]; 
            }
        }

        colours[v] = my_colour;

    });

    




    return colours;
}


/*
    Converts parents into a graph
*/
template <typename graph, typename T>
graph convert_parents_to_graph(graph G, parlay::sequence<T> parents)
{
    parlay::sequence<T> vertices = parlay::tabulate(parents.size(), [&] (T v) {return v;});

    G = parlay::map(vertices, [&] (T v) {
        parlay::sequence<T> temp = parlay::tabulate(1, [&] (T v) {return v;});
        temp[0] = parents[v];
        return temp;
    });

    G = graph_utils<T>::symmetrize(G);
    /*
        // delete parent pointers
        parlay::parallel_for(0, vertices.size(), [&] (T v) {
            T parent = parents[v];
            auto edgelist = G[v];
            for(T i = 0; i < edgelist.size(); i++)
            {
                if(edgelist[i] == parent)
                {
                    std::swap(edgelist[i], edgelist[edgelist.size()-1]);
                    edgelist.pop_back();
                    return;
                }
            }
        });
    */

    return G;
}