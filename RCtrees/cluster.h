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
const short is_MIS = 1024;
const short is_candidate = 2048;
const short internal = 8192;


const char neighbour_type = 1;
const char parent_type = 2;
const char child_type = 4;
const char edge_type = 8;

const int max_neighbours = 6;

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
    T indices[max_neighbours * 2]; 
    T index = -1;
    T colour = -1;
    T data = 0;
    std::atomic<T> height;
    static const short size = max_neighbours*2;
    short state = 0; // unaffected, affected or update eligible
    char types[max_neighbours * 2] = {};

    // Default constructor
    cluster() : height(0) {}

    // Copy constructor only for initialization when creating a sequence
    cluster(const cluster& other)
        : index(other.index),
          colour(other.colour),
          data(other.data),
          height(other.height.load()), // Load the atomic value
          state(other.state)
    {    
        for(uint i = 0; i < this->size; i++)
        {
            indices[i] = -1;
            types[i] = 0;
        }
        return;
    }

    // Method to get colour based off memory address
    unsigned long get_default_colour(void) const
    {
        return this->index;
    }

    // To be used only for connecting the base_edges to base_nodes
    void add_initial_neighbours(T V,  T W)
    {
        this->indices[0] = V;
        this->indices[1] = W;
        this->types[0] =neighbour_type;
        this->types[1] =neighbour_type;
    }

    // index of neighbour if successful, -1 if not
    // Note, if that pointer already exists, set its type to neighbour
    // and returns that index
    // index i+1 will contain edge for easier management
    short add_neighbour(T neighbour, T edge)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->indices[i] == -1 || this->indices[i] == neighbour)
            {
                this->indices[i] = neighbour;
                this->indices[i+1] = edge;
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
    short remove_neighbour(T neighbour)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->indices[i] == neighbour)
            {
                this->types[i]&=(~neighbour_type);
                return i;
            }
        }
        return -1;
    }

    short change_to_neighbour(T neighbour)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->indices[i] == neighbour)
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
    short change_to_child(T child)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->indices[i] == child)
            {
                this->types[i]=child_type;
                return i;
            }
        }
        return -1;
    }

    /**
     * designates a neighbour as a parent
    */
    short set_parent(T node)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->indices[i] == node)
            {
                this->types[i]=parent_type;
                this->types[i]&=(~neighbour_type);
                return i;
            }
        }
        return -1;
    }

    T get_parent(void)
    {
        for(short i = 0; i < this->size; i++)
        {
            if(this->types[i]&parent_type)
            {
                return this->indices[i];
            }
        }
        return -1;
    }

    void set_MIS(bool val)
    {
        static const short flags = 0xffff & is_MIS;
        if(val == true)
            this->state |= flags;
        else
            this->state &= (~flags);
    }

    // set the neighbours MIS to val
    // useful when colouring
    void set_neighbour_mis(bool val, parlay::sequence<cluster<T>>& all_clusters)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->indices[i] > -1)
                all_clusters[this->indices[i]].set_MIS(val);
        }
    }

    void set_neighbour_colour(T input_colour,  parlay::sequence<cluster<T>>& all_clusters)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->indices[i] > -1)
                all_clusters[this->indices[i]].colour = input_colour;
        }   
    }

    bool is_max_neighbour_colour( parlay::sequence<cluster<T>>& all_clusters)
    {
        T max_colour = this->colour;
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->indices[i] > -1 && all_clusters[this->indices[i]].colour >= max_colour && all_clusters[this->indices[i]].state & is_candidate)
                max_colour = all_clusters[this->indices[i]].colour;
        }
        if(max_colour == this->colour)
            return true;
        return false;
    }

    /*
        return true if any neighbours are in MIS
        else return false
    */
    bool get_neighbour_MIS( parlay::sequence<cluster<T>>& all_clusters)
    {
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->indices[i] > -1 && all_clusters[this->indices[i]].state & is_MIS && all_clusters[this->indices[i]].state & is_candidate)
                return true;
        }
        return false;
    }

    short get_neighbour_count(void)
    {
        short num_neighbours = 0;
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->types[i] & neighbour_type)
                num_neighbours++;
        }

        return num_neighbours;
    }

    void get_two_neighbours_edges(T& neighbour1, T& edge1, T& neighbour2, T& edge2)
    {
        bool first_neighbour_found = false;
        for(short i = 0; i < this->size; i+=2)
        {
            if(this->types[i]&neighbour_type)
            {
                if(first_neighbour_found == false)
                {
                    first_neighbour_found = true;
                    neighbour1 = this->indices[i];
                    edge1 = this->indices[i+1];
                }
                else
                {
                    neighbour2 = this->indices[i];
                    edge2 = this->indices[i+1];
                    return;   
                }
            }
        }
    }

    /*
        Mainly used when doing compress

    */
    short overwrite_neighbour(T old_neighbour, T new_neighbour, T new_edge)
    {

        for(short i = 0; i < this->size; i+=2)
        {
            if(this->indices[i] == old_neighbour)
            {
                this->indices[i] = new_neighbour;
                this->types[i]|=neighbour_type;
                
                this->indices[i+1] = new_edge;
                this->types[i+1]|= edge_type;
                this->types[i+1]&= (~neighbour_type);
                
                return i;
            }
        }

        return -1;
    }

    void print( parlay::sequence<cluster<T>>& all_clusters)
    {
        std::cout << "Index: " << this->index;
        std::cout << "  ";
        if(this->state&live)
            std::cout << "live ";
        else
            std::cout << "dead ";
        
        if(this->state & base_vertex)
        {
            for(uint i = 0; i < this->size; i+=2)
                {
                    if(this->indices[i] == -1)
                        continue;
                    auto cluster = &all_clusters[this->indices[i]];
                    std::cout << cluster->index;
                    std::cout << "(";
                    std::cout << (&all_clusters[this->indices[i+1]])->index;
                    std::cout << ") ";
                    auto ptr_type = this->types[i];
                    if(ptr_type & neighbour_type)
                        std::cout << "N";
                    if(ptr_type & edge_type)
                        std::cout << "C";
                    if(ptr_type & parent_type)
                        std::cout << "P";
                    if(ptr_type & child_type)
                        std::cout << "C";
                    std::cout << "   ";
                }
        }
        else
        {
            for(uint i = 0; i < this->size; i+=1)
            {
                if(this->indices[i] == -1)
                    continue;
                auto cluster = &all_clusters[this->indices[i]];
                std::cout << cluster->index;
                std::cout << "(";
                auto ptr_type = this->types[i];
                if(ptr_type & neighbour_type)
                    std::cout << "N";
                if(ptr_type & edge_type)
                    std::cout << "C";
                if(ptr_type & parent_type)
                    std::cout << "P";
                if(ptr_type & child_type)
                    std::cout << "C";
                std::cout << "   ";
            }
        }
        // std::cout << " colour(" << this->colour << ") ";
        if(this->state & binary_cluster)
            std::cout << "compress ";
        if(this->state & unary_cluster)
            std::cout << "rake ";
        if(this->state & nullary_cluster)
            std::cout << "root ";

        std::cout << std::endl;
    }
};
