#include <atomic>
#include "/scratch/parlaylib/include/parlay/primitives.h"
#include "/scratch/parlaylib/include/parlay/sequence.h"
#include "/scratch/parlaylib/examples/helper/graph_utils.h"
#include <set>

// template <typename T>
// class bad_set // A horrible implementation, only useful for very very small sets (such as colours)
// {
//     private:
//         parlay::sequence<T> elements;

//         int get_element_index(T element)
//         {
//             for(int i = 0; i < this->elements.size(); i++)
//             {
//                 if(this->elements[i] == element)
//                     return i;
//             }
//             return -1;
//         }

//     public:
//         bad_set(void)
//         {
//             return;
//         }
//         bad_set(parlay::sequence<T> input_elements)
//         {
//             this->elements = input_elements;
//             return;
//         }
//         bad_set(T arr[])
//         {
//             auto size_arr = sizeof(arr)/sizeof(arr[0]); 
//             for(auto i = 0;i < size_arr; i++)
//             {
//                 this->elements.push_back(arr[i]);
//             }
//         }

//         bool is_member(T element)
//         {
//             bool retval = false;
//             for(uint i = 0; i < this->num_elementsize(); i++)
//             {
//                 if(element == this->elements[i])
//                 {
//                     return true;
//                 }
//             }
//             return false;
//         }

//         void add_member(T element)
//         {
//             this->elements.push_back(element);
//         }

//         void remove_member(T element)
//         {
//             int index = this->get_element_index(element);
//             int last_index = this->elements.size() - 1;
//             if(index >= 0 && last_index >= 0)
//             {
//                 std::swap(this->elements[index], this->elements[last_index]);
//                 this->elements.pop_back();
//             }
//         }

//         void add_set(parlay::sequence<T> new_elements)
//         {
//             for(auto i = 0; i < new_elements.size(); i++)
//             {
//                 this->add_member(new_elements[i]);
//             }
//         }

//         void subtract_set(parlay::sequence<T> new_elements)
//         {
//             for(auto i = 0; i < new_elements.size(); i++)
//             {
//                 this->remove_member(new_elements[i]);
//             }
//         }

//         parlay::sequence<T> get_elements(void)
//         {
//             return this->elements;
//         }

//         auto get_cound(void)
//         {
//             return this->elements.size();
//         }
// };


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


// if (verbose)
    // {
    //     for(uint i = 0; i < G.size(); i++)
    //     {
    //         std::cout << i << "(" << parents[i] << ")\t\t";
    //     }
    //     std::cout << std::endl;
    // }

    // if (verbose)
    // {
    //     for(uint i = 0; i < G.size(); i++)
    //     {
    //         std::cout << original_colours[i] << "\t\t";
    //     }
    //     std::cout << std::endl;
    // }

template <typename T>
std::string coloured(T colour)
{
    switch (colour) {
        case 0:
            return "\033[1;34m0\033[0m"; // Blue
        case 1:
            return "\033[1;32m1\033[0m"; // Green
        case 2:
            return "\033[1;36m2\033[0m"; // Cyan
        case 3:
            return "\033[1;33m3\033[0m"; // Yellow
        case 4:
            return "\033[1;36m4\033[0m"; // Cyan
        case 5:
            return "\033[1;35m5\033[0m"; // Magenta
        case 6:
            return "\033[1;31m6\033[0m"; // Red
        default:
            return "\033[0m7\033[0m"; // Unknown colour (reset colour)
    }

}

/*
    Convert 6 colour tree into 3 colour tree
*/
template <typename graph, typename T>
parlay::sequence<T> convert_6_to_3_tree(graph G, parlay::sequence<T> parents, parlay::sequence<T> original_colours, bool verbose = false)
{   

    static const int invalid_index = 9087;    

    static const T uncoloured = (T) -1;

    parlay::sequence<T> colours = parlay::tabulate(parents.size(), [&] (T v) {
        return original_colours[v];
    });

    std::cout << "Invalid index: " << invalid_index << " parentclr: " << coloured(colours[parents[invalid_index]]) << " colour: " << coloured(colours[invalid_index]) << " ";
    std::cout << "Children: ";
    for(uint i = 0; i < G[invalid_index].size(); i++)
        std:: cout << coloured(colours[G[invalid_index][i]]) << " ";
    std::cout << std::endl;


    // each node (except root) is recolored as its parent
    parlay::parallel_for(0, parents.size(), [&] (T v){
        if(v != parents[v])
            colours[v] = original_colours[parents[v]];
    });

    std::cout << "Invalid index: " << invalid_index << " parentclr: " << coloured(colours[parents[invalid_index]]) << " colour: " << coloured(colours[invalid_index]) << " ";
    std::cout << "Children: ";
    for(uint i = 0; i < G[invalid_index].size(); i++)
        std:: cout << coloured(colours[G[invalid_index][i]]) << " ";
    std::cout << std::endl;

    parlay::sequence<T> colours_copy = parlay::tabulate(parents.size(), [&] (T v) {
        return colours[v];
    });

    parlay::parallel_for(0, parents.size(), [&] (T v){
        if(v == parents[v]) //its a root
        {
            std::set<int> Set;
            Set.insert(0);
            Set.insert(1);
            Set.insert(2);

            // All children must have the same colour
            auto children = G[v];

            auto children_colour = colours_copy[children[0]];

            Set.erase(children_colour);

            colours[v] = *Set.begin(); // sets are ordered, first is minimum
            
        }    
    }); //so far so good

    std::cout << "Invalid index: " << invalid_index << " parentclr: " << coloured(colours[parents[invalid_index]]) << " colour: " << coloured(colours[invalid_index]) << " ";
    std::cout << "Children: ";
    for(uint i = 0; i < G[invalid_index].size(); i++)
        std:: cout << coloured(colours[G[invalid_index][i]]) << " ";
    std::cout << std::endl;

    parlay::sequence<bool> Vprime = parlay::tabulate(parents.size(), [&] (T v) {return false;});

    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(colours[v] > 2) // V2
        {
            auto edgelist = G[v];
            edgelist.push_back(parents[v]); // just for ease
            parlay::parallel_for(0, edgelist.size(), [&] (T edgenum) {
                auto w = edgelist[edgenum]; 
                if(colours[w] <= 2) //It is in V1
                {
                    Vprime[v] = true;
                }
            });
        }
    });

    // Cv <- uncoloured

    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(Vprime[v])
            colours[v] = 6; // 6 now represents uncoloured
    });

    colours_copy = parlay::tabulate(parents.size(), [&] (T v) {
        return colours[v];
    });

    std::cout << "Invalid index: " << invalid_index << " parentclr: " << coloured(colours[parents[invalid_index]]) << " colour: " << coloured(colours[invalid_index]) << " ";
    std::cout << "Children: ";
    for(uint i = 0; i < G[invalid_index].size(); i++)
        std:: cout << coloured(colours[G[invalid_index][i]]) << " ";
    std::cout << std::endl;

    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if (Vprime[v])
        {
            // if(Vprime[parents[v]] == false) // if father of v is not in V prime
            if(colours_copy[parents[v]] != 6)
            {

                bool invalid_colours[7] = {false, false, false, false, false, false, false};
                
                parlay::parallel_for(0, G[v].size(), [&] (T child_id) {
                    
                    invalid_colours[colours_copy[G[v][child_id]]] = true;
                });
                invalid_colours[colours_copy[parents[v]]] = true;

                std::set<int> Set;
                Set.insert(0);
                Set.insert(1);
                Set.insert(2);

                for(uint i = 0; i < 7; i++)
                {
                    if(invalid_colours[i])
                        Set.erase(i);
                }

                colours[v] = *Set.begin();

                Vprime[v] = false;
            }

        }
    });

    colours_copy = parlay::tabulate(parents.size(), [&] (T v) {
        return colours[v];
    });

    std::cout << "Invalid index: " << invalid_index << " parentclr: " << coloured(colours[parents[invalid_index]]) << " colour: " << coloured(colours[invalid_index]) << " ";
    std::cout << "Children: ";
    for(uint i = 0; i < G[invalid_index].size(); i++)
        std:: cout << coloured(colours[G[invalid_index][i]]) << " ";
    std::cout << std::endl;

    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if (Vprime[v])
        {
            

                bool invalid_colours[7] = {false, false, false, false, false, false, false};
                
                parlay::parallel_for(0, G[v].size(), [&] (T child_id) {
                    invalid_colours[colours_copy[G[v][child_id]]] = true;
                });
                invalid_colours[colours_copy[parents[v]]] = true;

                std::set<int> Set;
                Set.insert(0);
                Set.insert(1);
                Set.insert(2);

                for(uint i = 0; i < 7; i++)
                {
                    if(invalid_colours[i])
                        Set.erase(i);
                }

                colours[v] = *Set.begin();

            

        }
    });

    std::cout << "Invalid index: " << invalid_index << " parentclr: " << coloured(colours[parents[invalid_index]]) << " colour: " << coloured(colours[invalid_index]) << " ";
    std::cout << "Children: ";
    for(uint i = 0; i < G[invalid_index].size(); i++)
        std:: cout << coloured(colours[G[invalid_index][i]]) << " ";
    std::cout << std::endl;


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
        parlay::sequence<T> temp = parlay::tabulate(1, [&] (T i) {return i;});
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

template <typename T>
bool verify_colouring(parlay::sequence<T> parents, parlay::sequence<T> colours, T* invalid_index)
{

    bool is_valid = true;

    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(v != parents[v] && colours[v] == colours[parents[v]])
        {
                is_valid = false;
                *invalid_index = v;
        }
    });

    return is_valid;
}

/*
    returns edges, you gotta subtract them yourself later
*/
template<typename graph, typename T>
parlay::sequence<T> find_forest(graph G, const int max_degree)
{
    parlay::sequence<T> edges = parlay::tabulate(G.size(), [&] (T v) {
        return -1; // an edge of -1 means its not in the forest
    });

    parlay::sequence<T> roots = parlay::tabulate(G.size(), [&] (T v) {
        return -1; // -1 means this vertex is note a root
    });

    parlay::sequence<T> vertices = parlay::tabulate(G.size(), [&] (T v) {
        return v; // -1 means this vertex is note a root
    });

    parlay::parallel_for(0, G.size(), [&] (T v) {
        auto edgelist = G[v];

        bool local_maxima = true;

        T max_neighbour = v;

        for(uint i = 0; i < std::min(edgelist.size(), max_degree); i++) // max_degree so it doesn't mess up asymptotic spans etc
        {
            if (edgelist[i] > max_neighbour)
            {
                max_neighbour = edgelist[i];
                local_maxima = false;
            }
        }

        if(local_maxima)
        {
            roots[v] = v;
        }
        else
            edges[v] = max_neighbour;
        
    });

    parlay::parallel_for(0, G.size(), [&] (T v) {
        if(roots[v])
        {
            auto edgelist = G[v];
            if(edgelist.size())
            {
                edges[v] = edgelist[0];
            }
        }
    });

    return edges;
}

template<typename graph, typename T>
graph subtract_edges(graph G, parlay::sequence<T> edges)
{

    parlay::parallel_for(0, G.size(), [&] (T v) {
        if(edges[v] != -1)
        {
            G[v] = parlay::filter(G[v], (T w) {
                return w != edges[v];
            });
            T w = edges[v];
            G[w] = parlay::filter(G[w], (T vp) {
                return vp != v;
            });
        }
    });

    return G;
}


template <typename graph, typename T>
bool verify_colouring(graph G, parlay::sequence<T> parents, parlay::sequence<T> colours, T* invalid_index)
{

    bool is_valid = true;

    parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(v != parents[v] && colours[v] == colours[parents[v]])
            {
                is_valid = false;
                *invalid_index = v;
            }
        auto edgelist = G[v];
        parlay::parallel_for(0, edgelist.size(), [&] (int i) {
            auto outgoing_node = edgelist[i];
            if(colours[v] == colours[outgoing_node])
            {    
                is_valid = false;
                *invalid_index = v;
            }
        });
    });

    return is_valid;
}


template <typename graph, typename T>
parlay::sequence<T> colour_graph(graph G, const int max_degree = 8)
{
    parlay::sequence<T> colour;

    parlay::sequence<T> vertex_id = parlay::tabulate(0, G.size(), [&] (T v) {
        return v;
    });

    parlay::sequence<T> parents = find_forest(G, max_degree); 

    parlay::sequence<T> valid_parent = parlay::map(parents, [&] (T v) {
        return v != -1;
    });

    parents = parlay::parallel_for(0, parents.size(), [&] (T v) {
        if(parents[v] == -1)
            parents[v] = v; // just make each unconnected root's parent itself, so we don't break the implementations of colouring graphs
    });



    colour = six_colour_rooted_tree(parents, vertex_id);

    // TODO remove edges to parents from graph
    
    // TODO 6 colour graph




    return colour;
}