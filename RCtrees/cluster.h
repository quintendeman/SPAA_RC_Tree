#include <atomic>
#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../examples/counting_sort.h"

const short empty_type = 0;
const short base_vertex = 1;
const short base_edge = 2;
const short unary_cluster = 4;
const short binary_cluster = 8;
const short nullary_cluster = 16;
const short live = 256;
const short carrying_weight = 512;
const short internal = 8192;

const char neighbour_type = 1;
const char parent_type = 2;
const char child_type = 4;
const char edge_type = 8;

const int max_neighbours = 3;

/*
    This represents a cluster in an RC tree.
    All these variables might be excessive but these flags and such are necessary since will be relying on pointer chasing
*/
#include <atomic>
#include <cstring> // for std::memcpy if needed

template <typename T>
struct cluster
{
public:
    cluster<T>* ptrs[max_neighbours * 2]; 
    T index = -1;
    T colour = -1;
    T data = 0;
    std::atomic<T> height;
    std::atomic<T> counter; 
    static const short size = max_neighbours*2;
    short state = 0; // unaffected, affected or update eligible
    bool is_MIS = false;
    char types[max_neighbours * 2] = {};

    // Default constructor
    cluster() : height(0) {
        for(uint i = 0; i < this->size; i++)
        {
            this->ptrs[i] = nullptr;
            this->types[i] = 0;
        }
    }

    // Copy constructor only for initialization when creating a sequence
    cluster(const cluster& other)
        : index(other.index),
          colour(other.colour),
          data(other.data),
          height(other.height.load()), // Load the atomic value
          counter(other.counter.load()), // Load the atomic value
          state(other.state),
          is_MIS(other.is_MIS)
    {    
        for(uint i = 0; i < this->size; i++)
        {
            ptrs[i] = nullptr;
            types[i] = 0;
        }
        return;
    }

    // Method to get colour based off memory address
    unsigned long get_default_colour(void) const
    {
        return reinterpret_cast<T>(this);
    }

    // To be used only for connecting the base_edges to base_nodes
    void add_initial_neighbours(cluster<T>* V,  cluster<T>* W)
    {
        this->ptrs[0] = V;
        this->ptrs[1] = W;
        this->types[0] = neighbour_type;
        this->types[1] = neighbour_type;
    }

    // index of neighbour if successful, -1 if not
    // Note, if that pointer already exists, set its type to neighbour
    // and returns that index
    // index i+1 will contain edge for easier management
    short add_neighbour(cluster<T>* neighbour, cluster<T>* edge)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] == nullptr || this->ptrs[i] == neighbour)
            {
                this->ptrs[i] = neighbour;
                this->ptrs[i+1] = edge;
                this->types[i] |= neighbour_type;
                this->types[i+1] |= edge_type;
                return i;
            }
        }
        return -1;
    }

    /**
     * Removes a node as a neighbour and returns the index
    */
    short remove_neighbour(cluster<T>* neighbour)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == neighbour)
            {
                this->types[i]&=(~neighbour_type);
                return i;
            }
        }
        return -1;
    }

    short change_to_neighbour(cluster<T>* neighbour)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == neighbour)
            {
                this->types[i]&=(~child_type);
                this->types[i]&=(~edge_type);
                this->types[i]|=neighbour_type;
                return i;
            }
        }
        return -1;
    }

    // changes a ptr from neighbour_type or edge_type to child_type
    short change_to_child(cluster<T>* child)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == child)
            {
                this->types[i]&=(~neighbour_type);
                this->types[i]&=(~edge_type);
                this->types[i]|=child_type;
                return i;
            }
        }
        return -1;
    }

    /**
     * designates a neighbour as a parent
    */
    short set_parent(cluster<T>* node)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == node)
            {
                this->types[i]=parent_type;
                this->types[i]&=(~neighbour_type);
                return i;
            }
        }
        return -1;
    }

    cluster<T>* get_parent(void)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->types[i]&parent_type)
            {
                return this->ptrs[i];
            }
        }
        return nullptr;
    }

    // set the neighbours MIS to val
    // useful when colouring
    void set_neighbour_mis(bool val)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i])
                this->ptrs[i]->is_MIS = val;
        }
    }

    void set_neighbour_colour(T input_colour)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i])
                this->ptrs[i]->colour = input_colour;
        }   
    }

    bool is_max_neighbour_colour(void)
    {
        T max_colour = this->colour;
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] && this->ptrs[i]->colour > max_colour)
                max_colour = this->ptrs[i]->colour;
        }
        if(max_colour == this->colour)
            return true;
        return false;
    }

    /*
        return true if any neighbours are in MIS
        else return false
    */
    bool get_neighbour_MIS(void)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->ptrs[i] && this->ptrs[i]->is_MIS)
                return true;
        }
        return false;
    }

    short get_neighbour_count(void)
    {
        short num_neighbours = 0;
        for(short i = 0; i < this->size; i++)
        {
            if(this->types[i] & neighbour_type)
                num_neighbours++;
        }

        return num_neighbours;
    }

    void get_two_neighbours_edges(cluster<T>*& neighbour1, cluster<T>*& edge1, cluster<T>*& neighbour2, cluster<T>*& edge2)
    {
        bool first_neighbour_found = false;
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->types[i]&neighbour_type)
            {
                if(first_neighbour_found == false)
                {
                    first_neighbour_found = true;
                    neighbour1 = this->ptrs[i];
                    edge1 = this->ptrs[i+1];
                }
                else
                {
                    neighbour2 = this->ptrs[i];
                    edge2 = this->ptrs[i+1];
                    return;   
                }
            }
        }
    }

    /*
        Mainly used when doing compress

    */
    short overwrite_neighbour(cluster<T>* old_neighbour, cluster<T>* new_neighbour, cluster<T>* new_edge)
    {

        for(short i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] == old_neighbour)
            {
                this->ptrs[i] = new_neighbour;
                this->types[i]|=neighbour_type;
                
                this->ptrs[i+1] = new_edge;
                this->types[i+1]|= edge_type;
                this->types[i+1]&= (~neighbour_type);
                
                return i;
            }
        }

        return -1;
    }

    void print_as_edge(void)
    {
        if(this == nullptr)
            return;
        for(uint i = 0; i < this->size; i++)
        {
            if(this->ptrs[i] != nullptr)
                std::cout << " " << this->ptrs[i]->index;
        }
        std::cout << ((this->state & base_edge) ? ""  : ("[" + std::to_string(this->index) + "]"));
        return;
    }

    void print(void)
    {
        std::cout << "Index: " << this->index;
        std::cout << "  ";
        if(this->state&live)
            std::cout << "live ";
        else
            std::cout << "dead ";
        for(uint i = 0; i < this->size; i+=2)
        {
            auto cluster = this->ptrs[i];
            if(cluster == nullptr)
                continue;
            std::cout << cluster->index;
            std::cout << " ";
            auto ptr_type = this->types[i];
            if(ptr_type & neighbour_type)
                std::cout << "N";
            if(ptr_type & edge_type)
                std::cout << "E";
            if(ptr_type & parent_type)
                std::cout << "P";
            if(ptr_type & child_type)
                std::cout << "C";
            std::cout << "(";
            this->ptrs[i+1]->print_as_edge();
            std::cout << " ";
            ptr_type = this->types[i+1];
            if(ptr_type & neighbour_type)
                std::cout << "N";
            if(ptr_type & edge_type)
                std::cout << "E";
            if(ptr_type & parent_type)
                std::cout << "P";
            if(ptr_type & child_type)
                std::cout << "C";
            std::cout << ")   ";
        }
        // std::cout << " colour(" << this->colour << ") ";
        std::cout << std::endl;
    }
};
