#ifndef TERNARIZER_H
#define TERNARIZER_H

#include <iostream>
#include <sstream>
#include <cmath>
#include <parlay/primitives.h>
#include <parlay/hash_table.h>
#include <parlay/alloc.h>
#include <cstdint>






const double load_factor = 2.0f; // load factor for the hash table
const long extra_tern_node_factor = 2;

const bool debug_print = false;

const char edge_marked_for_one_sided = 1;
const char edge_marked_for_deletion = 1 << 1;
const char edge_marked_as_tail = 1 << 2;

void checkDuplicates(const parlay::sequence<std::pair<long,long>>& del_edges) {
    parlay::sequence<std::atomic<char>> flags = parlay::tabulate<std::atomic<char>>(4 * del_edges.size(), [] (long i) {return 0;});
    parlay::sequence<std::pair<long,long>> vals = parlay::sequence<std::pair<long,long>>(flags.size());

    parlay::parallel_for(0, del_edges.size(), [&] (long i) {
        long v = del_edges[i].first;
        long w = del_edges[i].second;
        if(v < w)
            std::swap(v, w);
        size_t expected_location = ((v * 1000003l) ^ (w * 10000000019l)) % flags.size();
        
        while(true) {
            std::atomic<char>& flag = flags[expected_location];
            if(!flag) {
                char expected = 0;
                char desired = 1;
                bool success = flag.compare_exchange_strong(expected, desired);
                if(success){
                    vals[expected_location] = {v , w};
                    break;
                }
            } else if(vals[expected_location].first == v && vals[expected_location].second == w) {
                char expected = 1;
                char desired = 2;
                bool success = flag.compare_exchange_strong(expected, desired);
                break;
            } else
                expected_location = (expected_location + 1) % flags.size();
        
        }

    });

    parlay::parallel_for(0, flags.size(), [&] (long i) {
        if(flags[i] >= 2) {
            std::cout << vals[i].first << " " << vals[i].second << " is present at least twice" << std::endl;
        }
        assert(flags[i] < 2);
    });
}

void checkDuplicates(const parlay::sequence<long>& input) {
    parlay::sequence<std::atomic<char>> flags = parlay::tabulate<std::atomic<char>>(4 * input.size(), [] (long i) {return 0;});
    parlay::sequence<long> vals = parlay::sequence<long>(flags.size());

    parlay::parallel_for(0, input.size(), [&] (long i) {
        
        long v = input[i];

        size_t expected_location = ((v * 1000003l) ^ (v * 10000000019l)) % flags.size();
        
        while(true) {
            std::atomic<char>& flag = flags[expected_location];
            if(!flag) {
                char expected = 0;
                char desired = 1;
                bool success = flag.compare_exchange_strong(expected, desired);
                if(success){
                    vals[expected_location] = v;
                    break;
                }
            } else if(vals[expected_location] == v) {
                char expected = 1;
                char desired = 2;
                bool success = flag.compare_exchange_strong(expected, desired);
                break;
            } else
                expected_location = (expected_location + 1) % flags.size();
        
        }

    });

    parlay::parallel_for(0, flags.size(), [&] (long i) {
        if(flags[i] >= 2) {
            std::cout << vals[i] << " is present at least twice" << std::endl;
        }
        assert(flags[i] < 2);
    });
}


namespace { // we already have a tern_node in adjacency linked list but.. I couldn't come up with a better name
    

    template <typename T, typename D>
    struct tern_node{
        // parlay::sequence<std::tuple<T,T,D>> adds;
        // parlay::sequence<std::pair<T,T>> dels;
        T dummy_index = -1; // the actual node that I am
        T owner = -1; // The owner i.e. whose linked list I am in
        T tail_outgoing_index = -1;
        T tail_node_index = -1;
        T new_tail_node_index = -1;
        T adjacency_count = 0;
        T spill = 0;
        D tail_data;
        std::array<T, 3> outgoing_edgeindices = {-1};
        std::array<D, 3> outgoing_weights;
        char random_colour = 0.0;
        bool real = true; // TODO unnecessary?
        bool marked_for_deletion = false;
        bool shortlisted = false; // pun not intended
        bool last_to_contract = false;
        bool tail_changed = false;
        std::array<bool, 3> new_conn = {false, false, false};
        
        void clear_data(D identity)
        {
            this->owner = -1;
            // assert(!real); // TODO remove
            this->marked_for_deletion = false;
            this->shortlisted = false;
            this->random_colour = 0;
            this->outgoing_edgeindices[0] = this->outgoing_edgeindices[1] = this->outgoing_edgeindices[2] = -1;
            this->outgoing_weights[0] = this->outgoing_weights[1] = this->outgoing_weights[2] = identity;

            // this->head_index = -1;
            this->tail_outgoing_index = -1;
            this->tail_node_index = -1;
            this->new_tail_node_index = -1;
            this->adjacency_count = 0;
            this->spill = 0;
            this->tail_changed = false;
            this->last_to_contract = false;
            new_conn[0] = false;
            new_conn[1] = false;
            new_conn[2] = false;
        }

        short find_min_valid_index(void) {
            for(short a = 0; a < 3; ++a) {
                if(this->outgoing_edgeindices[a] == -1)
                    return a;
            }
            return -1;
        }

        friend std::ostream& operator<<(std::ostream& os, const tern_node& obj) {
            os << "[" << obj.dummy_index << "/" << obj.owner << "]: " << obj.dummy_index << " ";
            os << obj.outgoing_edgeindices[0] << "/" << obj.outgoing_weights[0] << " ";
            os << obj.outgoing_edgeindices[1] << "/" << obj.outgoing_weights[1] << " ";
            os << obj.outgoing_edgeindices[2] << "/" << obj.outgoing_weights[2] << " ";
            os << (int) obj.random_colour << " ";
            os << "adj: " << obj.adjacency_count << " ";
            os << obj.tail_node_index << " " << obj.new_tail_node_index << " " << obj.tail_outgoing_index << " ";
            if(obj.real)
                os << "REAL ";
            else
                os << "DMMY ";
            if(obj.marked_for_deletion)
                os << "M";
            if(obj.shortlisted)
                os << "S";
            os << "]";
            return os;
        }
    };

    // basic edge in hash table that helps us find dummy vertices
    // helps map an edge to two nodes in simplified tree
    template <typename T>
    struct basic_edge {
        T v = -2;
        T w = -1;
        T dummy_v = -2;
        T dummy_w = -1;
        char ind = -3;
        std::atomic<char> state;


        friend std::ostream& operator<<(std::ostream& os, const basic_edge<T>& edge) {
        return os << "basic_edge { v: " << edge.v
                << ", w: " << edge.w
                << ", dummy_v: " << edge.dummy_v
                << ", dummy_w: " << edge.dummy_w
                << ", state: " << (int) edge.state 
                << " }";
        }
    };

    template<typename T, typename D>
    struct hash_nodic {
        using kType = basic_edge<T>; // Key type: the tuple
        using eType = void*; // Value type: pointer to the tuple

        // Define an empty value as nullptr
        eType empty() { return nullptr; }

        // Extract the key (tuple) from the value (pointer)
        kType getKey(eType v) {
            return *static_cast<kType*>(v); // Dereference pointer to get the tuple
        }

            // Hash function for the tuple
        size_t hash(kType V) {
            // Hash each element of the tuple and combine
            assert(V.v != -1 && V.w != -1); 
            T v = V.v;
            T w = V.w;
            // if(v < w)
            //     std::swap(v, w);
            static const unsigned long random_prime1 = 1000003l;
            static const unsigned long random_prime2 = 10000000019l;

            size_t h1 = v * random_prime1;
            size_t h2 = w * random_prime2;
            // std::cout << "[HASH] returning " << (h1 ^ h2) << " for " << v << " -- " << w << std::endl;

            return h1 ^ h2; // Combine hashes using XOR 
            // return v + w; // Combine hashes using XOR 
        }

        // Comparison function for tuples
        int cmp(kType v1, kType v2) {
            // Extract values
            T val1 = v1.v, wal1 = v1.w;
            T val2 = v2.v, wal2 = v2.w;

            // Ensure (val, wal) is always in sorted order
            // if (val1 < wal1) std::swap(val1, wal1);
            // if (val2 < wal2) std::swap(val2, wal2);

            // std::cout << "cmp" << std::endl;

            // Lexicographic comparison: (val1, wal1) vs (val2, wal2)
            // if (val1 != val2) return (val1 < val2) ? -1 : 1;
            // if (wal1 != wal2) return (wal1 < wal2) ? -1 : 1;
            if(val1 != val2 || wal1 != wal2) return -1;
            return 0; // Equal
        }

        // Replace query: always return false (no replacement logic here)
        bool replaceQ(eType, eType) { return false; }

        // Update: return the current value
        eType update(eType v, eType) { return v; }

        // Atomic compare-and-swap (CAS) operation
        bool cas(eType* p, eType o, eType n) {
            return std::atomic_compare_exchange_strong_explicit(
            reinterpret_cast<std::atomic<eType>*>(p), &o, n, std::memory_order_relaxed, std::memory_order_relaxed);
        }
    };
}

size_t simple_hash(long v, long w) {
    static const size_t random_prime1 = 1000003ul;
    static const size_t random_prime2 = 10000000019ul;

    size_t h1 = v * random_prime1;
    size_t h2 = w * random_prime2;

    size_t x = h1 ^ h2;

    // Mix bits — inspired by MurmurHash finalizer
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;

    return x;
}


long roundUpPow2(long v) {
    if (v <= 0) return 1;
    // compute ceil(log2(v)) and then 2^that
    double lg  = std::log2(static_cast<double>(v));
    double clg = std::ceil(lg);
    long ret_val = static_cast<long>(std::pow(2.0, clg));
    if(ret_val < 8)
        ret_val = 8;
    return ret_val;
}

template <typename T, typename D>
class ternarizer{
    using wedge = std::tuple<T,T,D>;
    private:
        T max_index;
        T num_vertices; // including dummy graphs
        D identity;
        parlay::sequence<T> free_list; // free dummy indices;
        parlay::sequence<tern_node<T,D>> simplified_tree;
        // parlay::hashtable<hash_nodic<T,D>> HT;
        parlay::sequence<std::atomic<T>> tail_indices;
        // parlay::sequence<T> adjacents_count;
        parlay::sequence<parlay::sequence<std::pair<T, D>>> new_connections;
        parlay::sequence<basic_edge<T>> basic_edges;
        // parlay::sequence<T> basic_edge_free_list;
        // parlay::sequence<parlay::sequence<std::pair<T, D>>> new_connections_other_side;
        long seed = 0;
        // parlay::sequence<std::atomic<bool>> is_empty;
        parlay::sequence<std::atomic<T>> atomic_owner;

        parlay::sequence<std::atomic<bool>> key_valid;
        parlay::sequence<basic_edge<T>> key_edges;

        
        T cursor = 0; // location of free_list entries;

        T hash_table_size = 1024;
        // T basic_edge_cursor = 0;
    public:

        void print_state(void)
        {

            for(auto& n : this->simplified_tree)
            { 
                if(n.owner != n.dummy_index || n.real) {
                    if(n.outgoing_edgeindices[0] != -1 || n.outgoing_edgeindices[1] != -1 || n.outgoing_edgeindices[2] != -1)
                        std::cout << n << std::endl;
                }
            }

            for(long i = 0; i < this->hash_table_size; ++i) {
                if(this->key_valid[i])
                    std::cout << i << " " << this->key_edges[i] << std::endl;
            }
            return;
        }

        ternarizer(T max_index, D identity) 
        {
            this->identity = identity;
            this->max_index = max_index;
            this->num_vertices = extra_tern_node_factor * max_index;
            this->free_list = parlay::tabulate(this->num_vertices - max_index, [max_index] (T i) {
                return i + max_index; // N to 2N-1
            });
            this->simplified_tree = parlay::tabulate(this->num_vertices, [max_index, identity] (T i) {
                auto ret_tern_node = tern_node<T,D>();
                ret_tern_node.clear_data(identity);
                ret_tern_node.dummy_index = i;
                ret_tern_node.real = i >= max_index ? false : true;
                ret_tern_node.owner = i >= max_index ? -1 : i;
                // ret_tern_node.head_index = i >= max_index ? -1 : 0;
                // ret_tern_node.tail_index = i >= max_index ? -1 : 2;
                ret_tern_node.tail_node_index = i;
                return ret_tern_node;
            });

            // this->adjacents_count = parlay::sequence<T>(this->max_index, 0);

            this->new_connections = parlay::sequence<parlay::sequence<std::pair<T,D>>>(this->max_index);
            // this->is_empty = parlay::sequence<std::atomic<bool>>(this->max_index);
            this->atomic_owner = parlay::tabulate<std::atomic<T>>(this->num_vertices, [&] (T i) {return i;});

            T size_of_hashtable = roundUpPow2(this->num_vertices * 2);
            this->hash_table_size = size_of_hashtable;

            this->key_valid = parlay::tabulate(size_of_hashtable, [&] (T i) {return std::atomic<bool>(false);});

            this->key_edges = parlay::sequence<basic_edge<T>>(size_of_hashtable);
            
            this->cursor = 0; // start at N  (this just is for the free list)

            return;
        }

        ~ternarizer()
        {
            
        }

        auto get_owner(T v)
        {
            return this->simplified_tree[v].owner;
        }

        auto group_by_first(parlay::sequence<wedge>& wedge_list, bool ascending = true)
        {
            auto keyed_sequence = parlay::delayed_tabulate(wedge_list.size(), [&] (T i) {
                T& v = std::get<0>(wedge_list[i]);
                T& w = std::get<1>(wedge_list[i]);
                if(ascending) {
                    if(v < w)
                        std::swap(v, w);
                } else {
                    if(w < v)
                        std::swap(v, w);
                }
                    
                return std::tuple<T,wedge>(std::get<0>(wedge_list[i]), wedge_list[i]);
            });

            return std::move(group_by_key(keyed_sequence));
        }

        auto group_by_second(parlay::sequence<wedge>& wedge_list, bool ascending = true)
        {
            auto keyed_sequence = parlay::delayed_tabulate(wedge_list.size(), [&] (T i) {
                T& v = std::get<0>(wedge_list[i]);
                T& w = std::get<1>(wedge_list[i]);
                if(ascending) {
                    if(v < w)
                        std::swap(v, w);
                } else {
                    if(w < v)
                        std::swap(v, w);
                }
                return std::tuple<T,wedge>(std::get<1>(wedge_list[i]), wedge_list[i]);
            });
            return std::move(group_by_key(keyed_sequence));
        }

        auto group_by_first(parlay::sequence<std::pair<T,T>>& edge_list, bool ascending = true)
        {
            auto keyed_sequence = parlay::delayed_tabulate(edge_list.size(), [&] (T i) {
                T& v = edge_list[i].first;
                T& w = edge_list[i].second;
                if(ascending) {
                    if(v < w)
                        std::swap(v, w);
                } else {
                    if(w < v)
                        std::swap(v, w);
                }
                return std::tuple<T,std::pair<T,T>>(std::get<0>(edge_list[i]), edge_list[i]);
            });

            return std::move(group_by_key(keyed_sequence));
        }

        auto group_by_second(parlay::sequence<std::pair<T,T>>& edge_list, bool ascending = true)
        {
            auto keyed_sequence = parlay::delayed_tabulate(edge_list.size(), [&] (T i) {
                T& v = edge_list[i].first;
                T& w = edge_list[i].second;
                if(ascending) {
                    if(v < w)
                        std::swap(v, w);
                } else {
                    if(w < v)
                        std::swap(v, w);
                }
                return std::tuple<T,std::pair<T,T>>(std::get<1>(edge_list[i]), edge_list[i]);
            });

            return std::move(group_by_key(keyed_sequence));
        }

        parlay::sequence<T> collect_endpoints(const parlay::sequence<std::pair<T,T>>& edge_list, const bool mark = false) {
            parlay::parallel_for(0, edge_list.size(), [&] (T i) {
                const T& v = edge_list[i].first;
                const T& w = edge_list[i].second;
                this->atomic_owner[v] = i * 2;
                this->atomic_owner[w] = i * 2 + 1;
            });

            return parlay::filter(parlay::flatten(parlay::delayed_tabulate(edge_list.size(), [&] (T i) {
                const T& v_in = edge_list[i].first;
                const T& w_in = edge_list[i].second;
                T v = v_in;
                T w = w_in;
                if(this->atomic_owner[v] != i * 2)
                    v = -1;
                if(this->atomic_owner[w] != (i * 2 + 1))
                    w = -1;
                return std::array<T, 2>{v, w};
            })), [&] (T i) {
                if(mark && i != -1)
                    this->simplified_tree[i].marked_for_deletion = true;
                return i != -1;});
        }

        parlay::sequence<T> collect_endpoints(const parlay::sequence<wedge>& edge_list, const bool mark = false) {
            parlay::parallel_for(0, edge_list.size(), [&] (T i) {
                const T& v = std::get<0>(edge_list[i]);
                const T& w = std::get<1>(edge_list[i]);
                this->atomic_owner[v] = i * 2;
                this->atomic_owner[w] = i * 2 + 1;
            });

            return parlay::filter(parlay::flatten(parlay::delayed_tabulate(edge_list.size(), [&] (T i) {
                const T& v_in = std::get<0>(edge_list[i]);
                const T& w_in = std::get<1>(edge_list[i]);
                T v = v_in;
                T w = w_in;
                if(this->atomic_owner[v] != i * 2)
                    v = -1;
                if(this->atomic_owner[w] != (i * 2 + 1))
                    w = -1;
                return std::array<T, 2>{v, w};
            })), [&] (T i) {
                if(mark && i != -1)
                    this->simplified_tree[i].marked_for_deletion = true;
                return i != -1;});
        }

        basic_edge<T>& find_basic_edge_in_ht(T v, T w) // smaller helper function
        {
            // std::cout << "finding " + std::to_string(v) + " " + std::to_string(w) + "\n";
            size_t expected_location = simple_hash(v, w) & (hash_table_size - 1);
            size_t start = expected_location; 
            while(true) {
                if(!this->key_valid[expected_location]){
                    expected_location = (expected_location + 1) & (hash_table_size - 1);
                }
                else if(this->key_edges[expected_location].v == v && this->key_edges[expected_location].w == w)
                    return this->key_edges[expected_location];
                else
                    expected_location = (expected_location + 1) & (hash_table_size - 1);
                // if(expected_location == start)
                //     std::cout  << "Couldn't find " + std::to_string(v) + " " + std::to_string(w) + "\n";
            }
            
        }

        void insert_basic_edge_in_ht(T v, T w, T dummy_v, T dummy_w, short ind) 
        {
            // std::cout << "inserting " + std::to_string(v) + " " + std::to_string(w) + "\n";
            size_t expected_location = simple_hash(v, w) & (hash_table_size - 1);
            while(true) {
                std::atomic<bool>& valid_loc = this->key_valid[expected_location];
                if(valid_loc == false) {
                    bool expected = false;
                    bool desired = true;
                    bool success = valid_loc.compare_exchange_strong(expected, desired);
                    if(success) {
                        this->key_edges[expected_location].v = v;
                        this->key_edges[expected_location].w = w;
                        this->key_edges[expected_location].dummy_v = dummy_v;
                        this->key_edges[expected_location].dummy_w = dummy_w;
                        this->key_edges[expected_location].ind = ind;
                        this->key_edges[expected_location].state = 0;
                        
                        return;
                    }
                }
                expected_location = (expected_location + 1) & (hash_table_size - 1);
            }

        }

        void delete_basic_edge_in_ht(T v, T w) 
        {
            // std::cout << "deleting " + std::to_string(v) + " " + std::to_string(w) + "\n";
            size_t expected_location = simple_hash(v, w) & (hash_table_size - 1);
            while(true) {
                std::atomic<bool>& valid_loc = this->key_valid[expected_location];
                if(valid_loc == true) {
                    if(this->key_edges[expected_location].v == v &&
                        this->key_edges[expected_location].w == w) {
                        valid_loc = false; // should never be attempting to delete at the same time
                        return;
                    }
                }
                expected_location = (expected_location + 1) & (hash_table_size - 1);
            }
        }



        std::pair<parlay::sequence<std::pair<T,T>>, parlay::sequence<wedge>> add_edges(parlay::sequence<wedge>& add_edges)
        {

            parlay::sequence<wedge> ret_wedges;
            parlay::sequence<std::pair<T,T>> ret_del_pairs;

            auto by_first = group_by_first(add_edges, false);            
            parlay::parallel_for(0, by_first.size(), [&] (T I){
                auto& by_first_group_pair = by_first[I];
                T& first_index = by_first_group_pair.first;
                // this->is_empty[first_index] = false;
                // this->adjacents_count[first_index]+=by_first_group_pair.second.size();
                parlay::sequence<std::pair<T,D>> new_connections_up = parlay::sequence<std::pair<T,D>>(by_first_group_pair.second.size());
                parlay::parallel_for(0, by_first_group_pair.second.size(), [&] (T J) {
                    T& second_index = std::get<1>(by_first_group_pair.second[J]);
                    D& weight = std::get<2>(by_first_group_pair.second[J]);
                    new_connections_up[J] = {second_index, weight};
                });
                this->new_connections[first_index].append(new_connections_up);
            });
            by_first = group_by_first(add_edges, true);
            parlay::parallel_for(0, by_first.size(), [&] (T I){
                auto& by_first_group_pair = by_first[I];
                T& first_index = by_first_group_pair.first;
                // this->is_empty[first_index] = false;
                // this->adjacents_count[first_index]+=by_first_group_pair.second.size();
                parlay::sequence<std::pair<T,D>> new_connections_up = parlay::sequence<std::pair<T,D>>(by_first_group_pair.second.size());
                parlay::parallel_for(0, by_first_group_pair.second.size(), [&] (T J) {
                    T& second_index = std::get<1>(by_first_group_pair.second[J]);
                    D& weight = std::get<2>(by_first_group_pair.second[J]);
                    new_connections_up[J] = {second_index, weight};
                });
                this->new_connections[first_index].append(new_connections_up);
            });

            const auto endpoints = this->collect_endpoints(add_edges);

            const auto new_entries_counts = parlay::tabulate(endpoints.size(), [&] (T i) {
                const T& endpoint = endpoints[i];

                // T old_count = this->adjacents_count[endpoint];
                // T new_count = this->adjacents_count[endpoint] + this->new_connections[endpoint].size();
                T old_count = this->simplified_tree[endpoint].adjacency_count;
                T new_count = this->simplified_tree[endpoint].adjacency_count + this->new_connections[endpoint].size();
                T spill = new_count - old_count;
                if(new_count <= 3)
                    spill = 0;
                else if(old_count >= 3)
                    spill = new_count - old_count;
                else {
                    // assert(new_count > 3 && old_count < 3); // todo
                    spill = new_count - 3;
                    // spill = new_count - 3 - (3 - old_count);
                    // spill = new_count - (3 - old_count);
                }
                this->simplified_tree[endpoint].spill = spill;
    
                return spill;
            });

            const auto sums_pair = parlay::scan(new_entries_counts);
            const T total_count = sums_pair.second;
            const parlay::sequence<T>& offsets = sums_pair.first;

            assert(this->cursor + total_count < this->free_list.size());

            // collect identity edges
            auto identity_edges = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {
                const T& endpoint = endpoints[i];
                parlay::sequence<wedge> ret_iden_edges;
                const T& my_adjacents_count = this->simplified_tree[endpoint].adjacency_count;
                if(simplified_tree[endpoint].spill == 0){
                    // since spill = 0, can be at most 3 edges
                    for(int a = 0; a < this->new_connections[endpoint].size(); ++a) {
                        const T& other_side = this->new_connections[endpoint][a].first;
                        short ind = this->simplified_tree[endpoint].find_min_valid_index();
                        this->simplified_tree[endpoint].outgoing_edgeindices[ind] = other_side;
                        this->simplified_tree[endpoint].outgoing_weights[ind] = this->new_connections[endpoint][a].second;
                        // if(ind == -1) { // todo
                            // std::cout << simplified_tree[endpoint].spill << " " << this->new_connections[endpoint].size() << std::endl;
                        // }
                        // assert(ind >= 0); // todo
                        this->insert_basic_edge_in_ht(endpoint, other_side, endpoint, -1, ind);
                    }
                    this->simplified_tree[endpoint].tail_node_index = endpoint;
                    this->simplified_tree[endpoint].new_tail_node_index = endpoint;
                    // TODO, should we finalize the tail already?
                    return ret_iden_edges; // There is nothing to do
                }
                // // there was a spill
                
                T actual_end = this->new_connections[endpoint].size();
                T local_offset = this->cursor + offsets[i];
                T new_tail_free_index = local_offset + this->simplified_tree[endpoint].spill - 1;
                this->simplified_tree[endpoint].new_tail_node_index = this->free_list[new_tail_free_index];
                auto& ntni = this->simplified_tree[endpoint].new_tail_node_index;
                ret_iden_edges = parlay::tabulate(this->simplified_tree[endpoint].spill, [&] (T j) {
                    const T& other_side = this->new_connections[endpoint][j].first;
                    const D& other_weight = this->new_connections[endpoint][j].second;
                    T free_node_index = this->free_list[local_offset + j];
                    T prev_node_index = j > 0 ? this->free_list[local_offset + j - 1] : this->simplified_tree[endpoint].tail_node_index;
                    this->insert_basic_edge_in_ht(endpoint, other_side, free_node_index, -1, 1);
                    this->simplified_tree[free_node_index].outgoing_edgeindices[0] = prev_node_index;
                    this->simplified_tree[free_node_index].outgoing_weights[1] = other_weight;
                    this->simplified_tree[free_node_index].outgoing_weights[0] = this->identity;
                    // assert(free_node_index != prev_node_index); // todo
                    this->simplified_tree[free_node_index].owner = endpoint;
                    if(j) {
                        this->simplified_tree[prev_node_index].outgoing_edgeindices[2] = free_node_index;
                        this->simplified_tree[prev_node_index].outgoing_weights[2] = this->identity;
                    } 
                    // std::cout << "Went here" << std::endl;
                    return std::tuple<T,T,D>(free_node_index, prev_node_index, this->identity);
                });
                bool check_ = true;
                for(int a = this->simplified_tree[endpoint].spill; a < this->new_connections[endpoint].size(); ++a) {
                    // assert(check_); // todo
                    const T& other_side = this->new_connections[endpoint][a].first;
                    const D& other_weight = this->new_connections[endpoint][a].second;
                    short ind = this->simplified_tree[endpoint].find_min_valid_index();
                    // std::cout << "For " << endpoint << " min valid index " << ind << " with other_side/weight " << other_side << "/" << other_weight << std::endl;
                    // if(ind < 0) { // todo
                    //     std::cout << this->simplified_tree[endpoint].spill << " " << this->new_connections[endpoint].size() << " " <<  a << std::endl;
                    //     std::cout << this->simplified_tree[endpoint];
                    // }
                    // assert(ind >= 0); // todo
                    if(ind == 0 || ind == 1) {
                        // assert(endpoint != other_side); // todo
                        this->simplified_tree[endpoint].outgoing_edgeindices[ind] = other_side;
                        this->simplified_tree[endpoint].outgoing_weights[ind] = other_weight;
                        this->insert_basic_edge_in_ht(endpoint, other_side, endpoint, -1, ind);
                    } else if (ind == 2) {
                        // This would have been the last, now will be replaced
                        // assert(my_adjacents_count < 3); // todo
                        // assert(endpoint != this->free_list[local_offset]); // todo
                        // assert(ntni != other_side); // todo
                        this->simplified_tree[endpoint].outgoing_edgeindices[ind] = this->free_list[local_offset];
                        this->simplified_tree[endpoint].outgoing_weights[ind] = this->identity;
                        this->simplified_tree[ntni].outgoing_edgeindices[ind] = other_side;
                        this->simplified_tree[ntni].outgoing_weights[ind] = other_weight;
                        this->insert_basic_edge_in_ht(endpoint, other_side, ntni, -1, 2);
                        // this->simplified_tree[endpoint].tail_changed = true; // signify that a deletion is necessary
                    }  
                }

                return ret_iden_edges;
            }));



            parlay::parallel_for(0, endpoints.size(), [&] (T i) {
                const T& endpoint = endpoints[i];
                const T& my_adjacents_count = this->simplified_tree[endpoint].adjacency_count;
                auto& ntni = this->simplified_tree[endpoint].new_tail_node_index;
                if(my_adjacents_count >= 3) {
                    // meaning there was a tail already
                    // just update the entry in the HT for now, we are going to be deleting the OG tail anyways
                    // find the relevant edge
                    T v = endpoint;
                    T w = this->simplified_tree[this->simplified_tree[this->simplified_tree[endpoint].tail_node_index].outgoing_edgeindices[2]].owner;
                    // assert(v >= 0 && w >= 0); // todo
                    // if(v == w) {
                    //     std::cout << this->simplified_tree[endpoint] << std::endl;
                    //     std::cout << this->simplified_tree[this->simplified_tree[endpoint].tail_node_index] << std::endl;
                    // }
                    // assert(v != w);

                    if(w < v)
                        std::swap(v, w);
                    // std::cout << "Searching for " << std::to_string(v) + " " + std::to_string(w) + "\n";
                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w); // TODO can you insert and search for something that already exists in parallel
                    // std::cout << "Changed " << relevant_edge;
                    // assert(relevant_edge.dummy_v != -1 && relevant_edge.dummy_w != -1); // todo
                    if(this->simplified_tree[relevant_edge.dummy_v].owner == endpoint)
                        relevant_edge.dummy_v = this->simplified_tree[endpoint].new_tail_node_index;
                    else if(this->simplified_tree[relevant_edge.dummy_w].owner == endpoint)
                        relevant_edge.dummy_w = this->simplified_tree[endpoint].new_tail_node_index;
                    else {
                        std::cout << relevant_edge << std::endl;
                        assert(false && "Should never reach here"); // todo
                    }
                    // tail changed AND was populated, need to upgrade it
                    this->simplified_tree[ntni].outgoing_edgeindices[2] = this->simplified_tree[this->simplified_tree[endpoint].tail_node_index].outgoing_edgeindices[2];
                    this->simplified_tree[ntni].outgoing_weights[2] = this->simplified_tree[this->simplified_tree[endpoint].tail_node_index].outgoing_weights[2];
                    this->simplified_tree[endpoint].tail_changed = true; // signify that a deletion is necessary
                    // std::cout << " to " << relevant_edge << std::endl;
                }
            });

            // I have made the identity edges, and outgoing connections
            // lets connect the simple edges together 
            auto simple_edges = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {
                const T& endpoint = endpoints[i];
                parlay::sequence<wedge> ret_wedge_contrib = parlay::flatten(parlay::tabulate(this->new_connections[endpoint].size(), [&] (T j) {
                    parlay::sequence<wedge> ret_cont;
                    const T& other_side = this->new_connections[endpoint][j].first;
                    const D& other_weight = this->new_connections[endpoint][j].second;
                    if(other_side < endpoint)
                        return ret_cont;
                    // connect both sides blindly
                    auto& outgoing = this->find_basic_edge_in_ht(endpoint, other_side);
                    auto& incoming = this->find_basic_edge_in_ht(other_side, endpoint);

                    // connect
                    auto& this_node_index = outgoing.dummy_v;
                    auto& that_node_index = incoming.dummy_v;
                    auto& this_node = this->simplified_tree[this_node_index];
                    auto& that_node = this->simplified_tree[that_node_index];

                    const auto& this_nodes_outgoing_index = outgoing.ind;
                    // assert(this_nodes_outgoing_index != -1); // todo
                    const auto& that_nodes_outgoing_index = incoming.ind;
                    // assert(that_nodes_outgoing_index != -1); // todo

                    this_node.outgoing_edgeindices[this_nodes_outgoing_index] = that_node_index;
                    that_node.outgoing_edgeindices[that_nodes_outgoing_index] = this_node_index;

                    outgoing.dummy_w = incoming.dummy_v;
                    // std::cout << "returning " << outgoing.dummy_v << " " << outgoing.dummy_w << " " << this_node.outgoing_weights[this_nodes_outgoing_index] << std::endl;
                    ret_cont.push_back({outgoing.dummy_v, outgoing.dummy_w, this_node.outgoing_weights[this_nodes_outgoing_index]});
                    return ret_cont;

                }));

                return ret_wedge_contrib;
            }));
            // delete the edges in the HT
            parlay::parallel_for(0, endpoints.size(), [&] (T i) {
                const T& endpoint = endpoints[i];
                parlay::parallel_for(0, this->new_connections[endpoint].size(), [&] (T j) {
                    const T& other_side = this->new_connections[endpoint][j].first;
                    if(other_side < endpoint)
                        return;
                    this->delete_basic_edge_in_ht(other_side, endpoint);
                });
            });

            // collect delete edges
            auto delete_edges = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {
                parlay::sequence<std::pair<T,T>> ret_del_edges;

                const T& endpoint = endpoints[i];


                if(!this->simplified_tree[endpoint].tail_changed) {
                    return ret_del_edges;
                }

                const T& my_old_tail = this->simplified_tree[endpoint].tail_node_index;
                const T& other_side = this->simplified_tree[my_old_tail].outgoing_edgeindices[2];
                const T& other_endpoint = this->simplified_tree[other_side].owner;
                this->simplified_tree[endpoint].tail_outgoing_index = this->simplified_tree[my_old_tail].outgoing_edgeindices[2];
                this->simplified_tree[endpoint].tail_data = this->simplified_tree[my_old_tail].outgoing_weights[2];

                if(this->simplified_tree[other_endpoint].tail_changed) {
                    const auto& my_node = this->simplified_tree[endpoint];
                    const auto& their_node = this->simplified_tree[other_endpoint];
                    // was I connected to their tail
                    if(this->simplified_tree[their_node.tail_node_index].outgoing_edgeindices[2] == my_old_tail) {
                        // we need a tie break
                        if(other_endpoint < endpoint)
                            return ret_del_edges;
                        // std::cout << endpoint << " " << other_endpoint << std::endl;
                        // std::cout << my_node << " -------- " << their_node << std::endl;
                        // std::cout << "generated " << my_old_tail << " -- " << other_side <<  " del after tiebreak" << std::endl;
                        ret_del_edges.push_back({my_old_tail, other_side});
                    } else {
                        // std::cout << endpoint << " " << other_endpoint << std::endl;
                        // std::cout << my_node << " -------- " << their_node << std::endl;
                        // std::cout << "generated " << my_old_tail << " -- " << other_side <<  " del.." << std::endl;
                        ret_del_edges.push_back({my_old_tail, other_side});
                    }
                } else {
                    const auto& my_node = this->simplified_tree[endpoint];
                    const auto& their_node = this->simplified_tree[other_endpoint];
                    // std::cout << endpoint << " " << other_endpoint << std::endl;
                    // std::cout << my_node << " -------- " << their_node << std::endl;
                    // std::cout << "generated " << my_old_tail << " -- " << other_side <<  " del" << std::endl;
                    ret_del_edges.push_back({my_old_tail, other_side});
                }
                return ret_del_edges;
            }));

            // collect new tail
            auto new_tails = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {
                parlay::sequence<wedge> ret_seq;

                const T& endpoint = endpoints[i];

                if(!this->simplified_tree[endpoint].tail_changed) {
                    return ret_seq; // nothing to do
                }

                const T& my_old_tail = this->simplified_tree[endpoint].tail_node_index;
                const T& other_side = this->simplified_tree[endpoint].tail_outgoing_index;
                const T& other_endpoint = this->simplified_tree[other_side].owner;
                const D& rel_weight = this->simplified_tree[endpoint].tail_data;

                assert(other_side != -1); // TODO
                assert(other_endpoint != endpoint);
                
                T v = endpoint;
                T w = other_endpoint;
                
                

                if(w < v)
                    std::swap(v, w);


                auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                
                if(this->simplified_tree[other_endpoint].tail_changed) {
                    // std::cout << "Both had their tails changed" << std::endl;
                    const auto& my_node = this->simplified_tree[endpoint];
                    const auto& their_node = this->simplified_tree[other_endpoint];
                    // was I connected to their tail
                    if(their_node.tail_outgoing_index == my_old_tail) {
                        // std::cout << "And I was their tail" << std::endl;
                        // we need a tie break
                        if(other_endpoint < endpoint)
                            return ret_seq;
                        T rel_dummy_v = relevant_edge.dummy_v;
                        T rel_dummy_w = relevant_edge.dummy_w;
                        
                        if(this->simplified_tree[rel_dummy_v].owner != endpoint)
                            std::swap(rel_dummy_v, rel_dummy_w);

                        if(my_node.new_tail_node_index == rel_dummy_w || their_node.new_tail_node_index == rel_dummy_v) {//todo
                            std::cout << relevant_edge << std::endl;
                            std::cout << my_node << std::endl;
                            std::cout << their_node << std::endl;
                            std::cout << this->simplified_tree[my_node.new_tail_node_index] << std::endl;
                            std::cout << this->simplified_tree[their_node.new_tail_node_index] << std::endl;
                            
                        }
                        // assert(my_node.new_tail_node_index != rel_dummy_w); // todo
                        // assert(their_node.new_tail_node_index != rel_dummy_v);
                        
                        this->simplified_tree[my_node.new_tail_node_index].outgoing_edgeindices[2] = rel_dummy_w;
                        this->simplified_tree[my_node.new_tail_node_index].outgoing_weights[2] = my_node.tail_data;
                        this->simplified_tree[their_node.new_tail_node_index].outgoing_edgeindices[2] = rel_dummy_v;
                        this->simplified_tree[their_node.new_tail_node_index].outgoing_weights[2] = my_node.tail_data;
                        ret_seq.push_back({rel_dummy_v, rel_dummy_w, my_node.tail_data});
                    } else { // their tail did change but it doesn't affect me
                        // std::cout << "But I wasn't connected to their tail" << std::endl;
                        const auto& my_node = this->simplified_tree[endpoint];
                        auto& other_side_node = this->simplified_tree[other_side];
                        short matching_ind = -1;
                        for(short a = 0; a < 3; ++a) {
                            if(other_side_node.outgoing_edgeindices[a] == my_old_tail) {
                                matching_ind = a;
                                break;
                            }
                        }
                        // assert(matching_ind >= 0); // todo

                        other_side_node.outgoing_edgeindices[matching_ind] = my_node.new_tail_node_index;
                        this->simplified_tree[my_node.new_tail_node_index].outgoing_edgeindices[2] = other_side;
                        this->simplified_tree[my_node.new_tail_node_index].outgoing_weights[2] = my_node.tail_data;

                        // std::cout << "Returning " << my_node.new_tail_node_index << " " << other_side << " " << my_node.tail_data << std::endl;
    
                        ret_seq.push_back({my_node.new_tail_node_index, other_side, my_node.tail_data});    
                    }
                } else {
                    // std::cout << "I am connected to something that didn't even grow a tail" << std::endl;
                    const auto& my_node = this->simplified_tree[endpoint];
                    auto& other_side_node = this->simplified_tree[other_side];

                    short matching_ind = -1;
                    for(short a = 0; a < 3; ++a) {
                        if(other_side_node.outgoing_edgeindices[a] == my_old_tail) {
                            matching_ind = a;
                            break;
                        }
                    }
                    // assert(matching_ind != -1); // todo
                    other_side_node.outgoing_edgeindices[matching_ind] = my_node.new_tail_node_index;
                    this->simplified_tree[my_node.new_tail_node_index].outgoing_edgeindices[2] = other_side;
                    this->simplified_tree[my_node.new_tail_node_index].outgoing_weights[2] = my_node.tail_data;

                    // std::cout << "Returning " << my_node.new_tail_node_index << " " << other_side << " " << my_node.tail_data << std::endl;

                    ret_seq.push_back({my_node.new_tail_node_index, other_side, my_node.tail_data});

                }
                return ret_seq;
            }));



            // /****************************************Free*******************************************************************************/

            // // free
            parlay::parallel_for(0, endpoints.size(), [&] (T i) {
                const T& endpoint = endpoints[i];

                if(this->simplified_tree[endpoint].tail_changed) {
                    T local_offset = this->cursor + offsets[i];
                    T head_index = this->free_list[local_offset];
                    auto& old_tail_node = this->simplified_tree[this->simplified_tree[endpoint].tail_node_index];
                    old_tail_node.outgoing_edgeindices[2] = head_index;
                    old_tail_node.outgoing_weights[2] = this->identity;
                }

                this->simplified_tree[endpoint].tail_changed = false;
                this->simplified_tree[endpoint].tail_node_index = this->simplified_tree[endpoint].new_tail_node_index;
                this->simplified_tree[endpoint].adjacency_count += this->new_connections[endpoint].size();
                this->new_connections[endpoint].clear();

            });

            this->cursor += total_count;

            // std::cout << "Cursor changed to " << this->cursor << std::endl;
           

            ret_wedges.append(identity_edges);
            ret_wedges.append(simple_edges);
            ret_wedges.append(new_tails);

            return {delete_edges, ret_wedges};
        }

        void remove_duplicates(parlay::sequence<std::pair<T,T>>& dels, parlay::sequence<wedge>& adds) {
            parlay::sequence<std::atomic<char>> flags = parlay::tabulate<std::atomic<char>>(roundUpPow2(dels.size() * 1.5f), [] (T i) {return 0;});
            parlay::sequence<std::pair<T,T>> vals = parlay::sequence<std::pair<T,T>>(flags.size());

            parlay::parallel_for(0, dels.size(), [&] (T i) {
                if(dels[i].second < dels[i].first)
                    std::swap(dels[i].first, dels[i].second);
                long expected_location = simple_hash(dels[i].first, dels[i].second) & (flags.size() - 1);
                while(true) {
                    if(flags[expected_location] == 0) {
                        char expected = 0;
                        char desired = 1;
                        bool success = flags[expected_location].compare_exchange_strong(expected, desired);
                        if(success) {
                            vals[expected_location].first = dels[i].first;
                            vals[expected_location].second = dels[i].second;
                            return;
                        }
                    }
                    expected_location = (expected_location + 1) & (flags.size() - 1);
                }
            });

            adds = parlay::filter(adds, [&] (auto& e) {
                T v = std::get<0>(e);
                T w = std::get<1>(e);
                if(w < v)
                    std::swap(v, w);
                long expected_location = simple_hash(v, w) & (flags.size() - 1);
                while(true) {
                    if(flags[expected_location]) {
                        if(vals[expected_location].first == v && vals[expected_location].second == w) {
                            char expected = 1;
                            char desired = 2;
                            flags[expected_location].compare_exchange_strong(expected, desired);
                            return false;
                        }
                        expected_location = (expected_location + 1) & (flags.size() - 1);                        
                    } else
                        return true; // empty location means no hit
                }
            });

            auto indices = parlay::delayed_tabulate(flags.size(), [&] (T i){
                return i;
            });
            auto good_indices = parlay::filter(indices, [&] (T i) {
                return flags[indices[i]] == 1;
            });
            dels = parlay::tabulate(good_indices.size(), [&] (T i) {
                return vals[good_indices[i]];
            });
        }

        std::pair<parlay::sequence<std::pair<T,T>>, parlay::sequence<wedge>> delete_edges(parlay::sequence<std::pair<T,T>> del_edges)
        {

            parlay::sequence<wedge> ret_wedges;
            parlay::sequence<std::pair<T,T>> ret_pairs;

            auto by_first = group_by_first(del_edges, false);            
            parlay::parallel_for(0, by_first.size(), [&] (T I){
                auto& by_first_group_pair = by_first[I];
                T& first_index = by_first_group_pair.first;
                parlay::sequence<std::pair<T,D>> new_connections_up = parlay::sequence<std::pair<T,D>>(by_first_group_pair.second.size());
                parlay::parallel_for(0, by_first_group_pair.second.size(), [&] (T J) {
                    T& second_index = std::get<1>(by_first_group_pair.second[J]);
                    D& weight = this->identity;
                    new_connections_up[J] = {second_index, weight};
                });
                this->new_connections[first_index].append(new_connections_up);
            });
            by_first = group_by_first(del_edges, true);
            parlay::parallel_for(0, by_first.size(), [&] (T I){
                auto& by_first_group_pair = by_first[I];
                T& first_index = by_first_group_pair.first;
                parlay::sequence<std::pair<T,D>> new_connections_up = parlay::sequence<std::pair<T,D>>(by_first_group_pair.second.size());
                parlay::parallel_for(0, by_first_group_pair.second.size(), [&] (T J) {
                    T& second_index = std::get<1>(by_first_group_pair.second[J]);
                    D& weight = this->identity;
                    new_connections_up[J] = {second_index, weight};
                });
                this->new_connections[first_index].append(new_connections_up);
            });

            const auto endpoints = this->collect_endpoints(del_edges);
            
            // mark the relevant indices
            auto ternarized_edges = parlay::tabulate(del_edges.size(), [&] (T i) { // Definitely must be deleted -- also, fortunately, covers the cross edges
                T v = del_edges[i].first;
                T w = del_edges[i].second;
                if(w < v)
                    std::swap(v, w);
                auto& relevant_edge = this->find_basic_edge_in_ht(v, w);

                // in dummy_v, find the edge going to dummy_w and set it to -1
                // in dummy_w, find the edge going to dummy_v and set it to -1
                for(short a = 0; a < 3; a++) {
                    // assert(a != 3); // todo
                    T& outgoing_index = this->simplified_tree[relevant_edge.dummy_v].outgoing_edgeindices[a];
                    if(outgoing_index == relevant_edge.dummy_w) {
                        outgoing_index = -1;
                        break;
                    }
                }

                for(short a = 0; a < 3; a++) {
                    // assert(a != 3); // todo
                    T& outgoing_index = this->simplified_tree[relevant_edge.dummy_w].outgoing_edgeindices[a];
                    if(outgoing_index == relevant_edge.dummy_v) {
                        outgoing_index = -1;
                        break;
                    }
                }
                this->delete_basic_edge_in_ht(v, w);

                return std::pair<T,T>(relevant_edge.dummy_v, relevant_edge.dummy_w);
            });

            auto all_dummy_endpoints = this->collect_endpoints(ternarized_edges, true);

            auto del_iden_edges = parlay::flatten(parlay::tabulate(all_dummy_endpoints.size(), [&] (T i) {
                parlay::sequence<std::pair<T,T>> ret_pairs;
                auto& dummy_node = this->simplified_tree[all_dummy_endpoints[i]];
                dummy_node.shortlisted = true;
                if(dummy_node.real || this->simplified_tree[dummy_node.owner].tail_node_index == dummy_node.dummy_index)
                    return ret_pairs;
                auto& left_node = this->simplified_tree[dummy_node.outgoing_edgeindices[0]];
                // std::cout << blue << "That " << dummy_node.dummy_index << " " << dummy_node.outgoing_edgeindices[0] << reset << std::endl;
                ret_pairs.push_back({dummy_node.dummy_index, dummy_node.outgoing_edgeindices[0]});

                const T& right = dummy_node.outgoing_edgeindices[2];
                if(right == this->simplified_tree[dummy_node.owner].tail_node_index || this->simplified_tree[right].marked_for_deletion == false) {
                    // std::cout << red << "that weird addition " << dummy_node.dummy_index << " " << right << reset << std::endl;
                    ret_pairs.push_back({dummy_node.dummy_index, right});
                }

                // dummy_node.outgoing_edgeindices[0] = -1;
                // left_node.outgoing_edgeindices[2] = -1;
                // if(right != -1 && this->simplified_tree[right].owner == dummy_node.owner && this->simplified_tree[right].marked_for_deletion) {// isn't a cross edge
                //     std::cout << red << "that weird addition " << dummy_node.dummy_index << " " << right << reset << std::endl;
                //     ret_pairs.push_back({dummy_node.dummy_index, right});
                // }
                return ret_pairs;
            })); // These are all the identity edges that we will delete

            auto contraction_nodes = all_dummy_endpoints;

            long temp = 0;

            while(contraction_nodes.size()) {
                temp++;
                // std::cout << temp << " " << contraction_nodes.size() << std::endl;

                if(contraction_nodes.size() < 10 && temp > 50) {
                    for(auto& cn : contraction_nodes) {
                        std::cout << this->simplified_tree[cn] << std::endl;
                    }
                } 

                parlay::parallel_for(0, contraction_nodes.size(), [&] (T i) {
                    auto& curr_node = this->simplified_tree[contraction_nodes[i]];
                    curr_node.random_colour = simple_hash(temp, i) & 1;
                    return;
                });
                parlay::parallel_for(0, contraction_nodes.size(), [&] (T i) {
                    auto& curr_node = this->simplified_tree[contraction_nodes[i]];
                    if(curr_node.real) { // head
                        curr_node.shortlisted = false;
                        return;
                    }
                    if(this->simplified_tree[curr_node.owner].tail_node_index == curr_node.dummy_index) { // tail
                        curr_node.shortlisted = false;
                        return;
                    }

                    // am I the only one with random_colour == 0 in my left and right neighbours
                    T& left_index = curr_node.outgoing_edgeindices[0];
                    T& right_index = curr_node.outgoing_edgeindices[2];
                    auto left_colour = 0;
                    if(left_index == -1)
                        left_colour = 1;
                    else if(this->simplified_tree[left_index].marked_for_deletion == false || !this->simplified_tree[left_index].shortlisted)
                        left_colour = 1;
                    else
                        left_colour = this->simplified_tree[left_index].random_colour;
                    auto right_colour = 0;
                    if(right_index == -1)
                        right_colour = 1;
                    else if(this->simplified_tree[right_index].marked_for_deletion == false || !this->simplified_tree[right_index].shortlisted)
                        right_colour = 1;
                    else
                        right_colour = this->simplified_tree[right_index].random_colour;
                    
                    
                    
                    if(left_colour == 1 && right_colour == 1 && curr_node.random_colour == 0) {
                        // merge
                        // std::cout << "Shortlisted " << curr_node << std::endl;
                        // std::cout << "Changed  " << this->simplified_tree[left_index] << " to ";
                        this->simplified_tree[left_index].outgoing_edgeindices[2] = right_index;
                        // std::cout << this->simplified_tree[left_index] << std::endl;
                        
                        // std::cout << "Changed " << this->simplified_tree[right_index] << " to ";
                        this->simplified_tree[right_index].outgoing_edgeindices[0] = left_index;
                        // std::cout << this->simplified_tree[right_index] << std::endl;
                        
                        if((!this->simplified_tree[left_index].marked_for_deletion || this->simplified_tree[left_index].real) && (!this->simplified_tree[right_index].marked_for_deletion || this->simplified_tree[curr_node.owner].tail_node_index == right_index)) {
                            // std::cout << curr_node << " was last to contract " << std::endl;
                            curr_node.last_to_contract = true;
                            curr_node.shortlisted = false; // Do a special contraction later
                            return;
                        } else {
                            // std::cout << curr_node << " was not the last to contract " << std::endl;
                            // std::cout << (!this->simplified_tree[left_index].marked_for_deletion) << " " << (this->simplified_tree[left_index].real) << "   " << (!this->simplified_tree[right_index].marked_for_deletion) << " " << (curr_node.tail_node_index == right_index) << std::endl;
                            // std::cout << curr_node.tail_node_index << " " << right_index << std::endl;
                        }
                        curr_node.clear_data(this->identity);
                        curr_node.shortlisted = false;
                        return;
                    } else {
                        // if(contraction_nodes.size() <= 50) {
                        //     std::cout << curr_node << " lost on tie breaks and couldn't contract" << std::endl;
                        //     std::cout << left_colour << " " << right_colour << std::endl;
                        // }
                    }

                });
                contraction_nodes = parlay::filter(contraction_nodes, [&] (T i) {
                    return this->simplified_tree[i].shortlisted;
                });

                assert(temp < 256); // prevent infinite loops
            }

            auto linked_list_edges = parlay::flatten(parlay::tabulate(all_dummy_endpoints.size(), [&] (T i) {
                parlay::sequence<wedge> ret_wedges;
                auto& dummy_node = this->simplified_tree[all_dummy_endpoints[i]];
                if(!dummy_node.last_to_contract)
                    return ret_wedges;
                ret_wedges.push_back({dummy_node.outgoing_edgeindices[0], dummy_node.outgoing_edgeindices[2], this->identity});
                this->simplified_tree[dummy_node.outgoing_edgeindices[0]].outgoing_edgeindices[2] =  dummy_node.outgoing_edgeindices[2];
                this->simplified_tree[dummy_node.outgoing_edgeindices[2]].outgoing_edgeindices[0] =  dummy_node.outgoing_edgeindices[0];
                dummy_node.clear_data(this->identity);
                return ret_wedges;
            })); 

            // this->print_state(); 

            // Consolidate deficencies 
            parlay::parallel_for(0, endpoints.size(), [&] (T i) {
                const T& endpoint = endpoints[i];

                // first check tail
                auto& head_node = this->simplified_tree[endpoint];
                auto& my_tail_node = this->simplified_tree[head_node.tail_node_index];
                if(my_tail_node.dummy_index == endpoint) {
                    return; // Nothing to do, I didn't contract as my tail is already me
                } else if(head_node.adjacency_count <= 3) { // Nothing fancy needs to be done -- we had not spilled
                    return;
                }
                auto& left_of_tail_node = this->simplified_tree[my_tail_node.outgoing_edgeindices[0]];

                bool tail_top_available = true;
                bool tail_right_available = true;
                bool left_of_tail_available = true;
                bool left_of_left_of_tail_available = true;

                // std::cout << "Changed " << head_node << " to ";

                if(head_node.outgoing_edgeindices[0] == -1) {
                    // There was a deletion
                    if(my_tail_node.outgoing_edgeindices[1] != -1) {
                        // update edge at this side
                        // try to steal from tail
                        T v = endpoint;
                        T w = this->simplified_tree[my_tail_node.outgoing_edgeindices[1]].owner;
                        if(w < v)
                            std::swap(v, w);
                        // head_node.adds.push_back({endpoint, this->simplified_tree[my_tail_node.outgoing_edgeindices[1]].owner, my_tail_node.outgoing_weights[1]});
                        auto& relevant_edge = this->find_basic_edge_in_ht(v, w); 
                        char old_val = relevant_edge.state.fetch_add(1);
                        // assert(old_val < 2); // todo
                        head_node.outgoing_edgeindices[0] = my_tail_node.outgoing_edgeindices[1];
                        head_node.outgoing_weights[0] = my_tail_node.outgoing_weights[1];
                        my_tail_node.outgoing_edgeindices[1] = -1;
                        tail_top_available = false;
                        head_node.new_conn[0] = true;
                    } else if(my_tail_node.outgoing_edgeindices[2] != -1) { // try to steal from tail
                        T v = endpoint;
                        T w = this->simplified_tree[my_tail_node.outgoing_edgeindices[2]].owner;
                        if(w < v)
                            std::swap(v, w);
                        // head_node.adds.push_back({endpoint, this->simplified_tree[my_tail_node.outgoing_edgeindices[2]].owner, my_tail_node.outgoing_weights[2]});
                        auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                        char old_val = relevant_edge.state.fetch_add(1);
                        // assert(old_val < 2); // todo
                        head_node.outgoing_edgeindices[0] = my_tail_node.outgoing_edgeindices[2];
                        head_node.outgoing_weights[0] = my_tail_node.outgoing_weights[2];
                        my_tail_node.outgoing_edgeindices[2] = -1;
                        head_node.new_conn[0] = true;
                        tail_right_available = false;
                    } else if(left_of_tail_node.dummy_index != endpoint && left_of_tail_node.outgoing_edgeindices[1] != -1) {
                        T v = endpoint;
                        T w = this->simplified_tree[left_of_tail_node.outgoing_edgeindices[1]].owner;
                        if(w < v)
                            std::swap(v, w);
                        // head_node.adds.push_back({endpoint, this->simplified_tree[left_of_tail_node.outgoing_edgeindices[1]].owner, my_tail_node.outgoing_weights[1]});
                        auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                        char old_val = relevant_edge.state.fetch_add(1);
                        // assert(old_val < 2); // todo
                        head_node.outgoing_edgeindices[0] = left_of_tail_node.outgoing_edgeindices[1];
                        head_node.outgoing_weights[0] = left_of_tail_node.outgoing_weights[1];
                        left_of_tail_node.outgoing_edgeindices[1] = -1;      
                        left_of_tail_available = false;                  
                        head_node.new_conn[0] = true;
                    } else if (left_of_tail_node.dummy_index != endpoint) {
                        auto& left_of_left_of_tail_node = this->simplified_tree[left_of_tail_node.outgoing_edgeindices[0]];
                        if(left_of_left_of_tail_node.dummy_index != endpoint) {
                            T v = endpoint;
                            T w = this->simplified_tree[left_of_left_of_tail_node.outgoing_edgeindices[1]].owner;
                            if(w < v)
                                std::swap(v, w);
                            // head_node.adds.push_back({endpoint, this->simplified_tree[left_of_left_of_tail_node.outgoing_edgeindices[1]].owner, left_of_left_of_tail_node.outgoing_weights[1]});
                            auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                            char old_val = relevant_edge.state.fetch_add(1);
                            // assert(old_val < 2); // todo
                            head_node.outgoing_edgeindices[0] = left_of_left_of_tail_node.outgoing_edgeindices[1];
                            head_node.outgoing_weights[0] = left_of_left_of_tail_node.outgoing_weights[1];
                            left_of_left_of_tail_node.outgoing_edgeindices[1] = -1;  
                            left_of_left_of_tail_available = false;
                            head_node.new_conn[0] = true;
                        }
                    }
                }

                if(head_node.outgoing_edgeindices[1] == -1) {
                    // There was a deletion
                    if(my_tail_node.outgoing_edgeindices[1] != -1 && tail_top_available) {
                        // update edge at this side
                        // try to steal from tail
                        T v = endpoint;
                        T w = this->simplified_tree[my_tail_node.outgoing_edgeindices[1]].owner;
                        if(w < v)
                            std::swap(v, w);
                        // head_node.adds.push_back({endpoint, this->simplified_tree[my_tail_node.outgoing_edgeindices[1]].owner, my_tail_node.outgoing_weights[1]});
                        auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                        char old_val = relevant_edge.state.fetch_add(1);
                        // assert(old_val < 2); // todo
                        head_node.outgoing_edgeindices[1] = my_tail_node.outgoing_edgeindices[1];
                        head_node.outgoing_weights[1] = my_tail_node.outgoing_weights[1];
                        head_node.new_conn[1] = true;
                        my_tail_node.outgoing_edgeindices[1] = -1;

                    } else if(my_tail_node.outgoing_edgeindices[2] != -1 && tail_right_available) { // try to steal from tail
                        T v = endpoint;
                        T w = this->simplified_tree[my_tail_node.outgoing_edgeindices[2]].owner;
                        if(w < v)
                            std::swap(v, w);
                        // head_node.adds.push_back({endpoint, this->simplified_tree[my_tail_node.outgoing_edgeindices[2]].owner, my_tail_node.outgoing_weights[2]});
                        auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                        char old_val = relevant_edge.state.fetch_add(1);
                        // assert(old_val < 2); // todo
                        head_node.outgoing_edgeindices[1] = my_tail_node.outgoing_edgeindices[2];
                        head_node.outgoing_weights[1] = my_tail_node.outgoing_weights[2];
                        head_node.new_conn[1] = true;
                        my_tail_node.outgoing_edgeindices[2] = -1;
                    } else if(left_of_tail_node.dummy_index != endpoint && left_of_tail_node.outgoing_edgeindices[1] != -1 && left_of_tail_available) {
                        T v = endpoint;
                        T w = this->simplified_tree[left_of_tail_node.outgoing_edgeindices[1]].owner;
                        if(w < v)
                            std::swap(v, w);
                        // head_node.adds.push_back({endpoint, this->simplified_tree[left_of_tail_node.outgoing_edgeindices[1]].owner, my_tail_node.outgoing_weights[1]});
                        auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                        char old_val = relevant_edge.state.fetch_add(1);
                        // assert(old_val < 2); // todo
                        head_node.outgoing_edgeindices[1] = left_of_tail_node.outgoing_edgeindices[1];
                        head_node.outgoing_weights[1] = left_of_tail_node.outgoing_weights[1];
                        head_node.new_conn[1] = true;
                        left_of_tail_node.outgoing_edgeindices[1] = -1;                        
                    } else if (left_of_tail_node.dummy_index != endpoint && left_of_left_of_tail_available) {
                        auto& left_of_left_of_tail_node = this->simplified_tree[left_of_tail_node.outgoing_edgeindices[0]];
                        if(left_of_left_of_tail_node.dummy_index != endpoint) {
                            T v = endpoint;
                            T w = this->simplified_tree[left_of_left_of_tail_node.outgoing_edgeindices[1]].owner;
                            if(w < v)
                                std::swap(v, w);
                            // head_node.adds.push_back({endpoint, this->simplified_tree[left_of_left_of_tail_node.outgoing_edgeindices[1]].owner, left_of_left_of_tail_node.outgoing_weights[1]});
                            auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                            char old_val = relevant_edge.state.fetch_add(1);
                            // assert(old_val < 2); // todo
                            head_node.outgoing_edgeindices[1] = left_of_left_of_tail_node.outgoing_edgeindices[1];
                            head_node.outgoing_weights[1] = left_of_left_of_tail_node.outgoing_weights[1];
                            head_node.new_conn[1] = true;
                            left_of_left_of_tail_node.outgoing_edgeindices[1] = -1;  
                        }
                    }
                }
                // std::cout << head_node << std::endl;
                // now that I have moved things to my head, let's consider the tail situation.
                return;
            });




            // collect del edges
            auto tail_dels_p1 = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {
                parlay::sequence<std::pair<T,T>> ret_del_seq;

                const T& endpoint = endpoints[i];

                // first check tail
                auto& head_node = this->simplified_tree[endpoint];
                
                if(head_node.outgoing_edgeindices[0] != -1) {
                    T v = endpoint;
                    T w = this->simplified_tree[head_node.outgoing_edgeindices[0]].owner;
                    if(w < v)
                        std::swap(v, w);
                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                    if(relevant_edge.state == 2) {
                        if(endpoint < this->simplified_tree[head_node.outgoing_edgeindices[0]].owner) {
                            // std::cout << "delpush.. " << relevant_edge.dummy_v << " " << relevant_edge.dummy_w << std::endl;
                            ret_del_seq.push_back({relevant_edge.dummy_v, relevant_edge.dummy_w}); 
                        }
                    } else if(relevant_edge.state == 1) {
                        char expected = 1;
                        char desired = 3;
                        bool success = relevant_edge.state.compare_exchange_strong(expected, desired);
                        if(success) {
                            // std::cout << "delpush " << relevant_edge.dummy_v << " " << relevant_edge.dummy_w << std::endl;
                            // std::cout << head_node << std::endl;
                            ret_del_seq.push_back({relevant_edge.dummy_v, relevant_edge.dummy_w}); 
                        }
                    }
                }
                if(head_node.outgoing_edgeindices[1] != -1) {
                    T v = endpoint;
                    T w = this->simplified_tree[head_node.outgoing_edgeindices[1]].owner;
                    if(w < v)
                        std::swap(v, w);
                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                    if(relevant_edge.state == 2) {
                        if(endpoint < this->simplified_tree[head_node.outgoing_edgeindices[1]].owner) {
                            // std::cout << "delpush.. " << relevant_edge.dummy_v << " " << relevant_edge.dummy_w << std::endl;
                            ret_del_seq.push_back({relevant_edge.dummy_v, relevant_edge.dummy_w}); 
                        }
                    } else if(relevant_edge.state == 1) {
                        char expected = 1;
                        char desired = 3;
                        bool success = relevant_edge.state.compare_exchange_strong(expected, desired);
                        if(success) {
                            // std::cout << "delpush " << relevant_edge.dummy_v << " " << relevant_edge.dummy_w << std::endl;
                            // std::cout << head_node << std::endl;
                            ret_del_seq.push_back({relevant_edge.dummy_v, relevant_edge.dummy_w}); 
                        }
                    }
                }
                return ret_del_seq;
            }));



            auto tail_adds_p1 = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {
                parlay::sequence<wedge> ret_seq;

                const T& endpoint = endpoints[i];
                auto& head_node = this->simplified_tree[endpoint];

                if(head_node.new_conn[0]) {
                    const short check_index = 0;
                    head_node.new_conn[check_index] = false;
                    T v = endpoint;
                    T w = this->simplified_tree[head_node.outgoing_edgeindices[check_index]].owner;
                    T their_endpoint = w;
                    if(w < v)
                        std::swap(v, w);
                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                    // assert(relevant_edge.state); // todo
                    // std::cout << relevant_edge << std::endl;
                    if(relevant_edge.state == 1 || relevant_edge.state == 3) { // that side should not be changed
                        // nothing fancy to do
                        if(this->simplified_tree[relevant_edge.dummy_v].owner == endpoint) {
                            relevant_edge.dummy_v = endpoint;
                            auto& their_node = this->simplified_tree[relevant_edge.dummy_w];
                            short ind = -1;
                            for(short a = 0; a < 3; ++a) {
                                // assert(a != 3); //
                                if(their_node.outgoing_edgeindices[a] != -1 && this->simplified_tree[their_node.outgoing_edgeindices[a]].owner == endpoint) {
                                    their_node.outgoing_edgeindices[a] = relevant_edge.dummy_v;
                                    break;
                                }
                            }   
                        }
                        else { 
                            relevant_edge.dummy_w = endpoint;
                            auto& their_node = this->simplified_tree[relevant_edge.dummy_v];
                            
                            for(short a = 0; a < 3; ++a) {
                                // assert(a != 3); // todo
                                if(their_node.outgoing_edgeindices[a] != -1 && this->simplified_tree[their_node.outgoing_edgeindices[a]].owner == endpoint) {
                                    their_node.outgoing_edgeindices[a] = relevant_edge.dummy_w;
                                    break;
                                }
                            }
                        }   
                        // std::cout << "Pushing " << relevant_edge.dummy_v << " " << relevant_edge.dummy_w << " " << head_node.outgoing_weights[check_index] << std::endl;
                        ret_seq.push_back({relevant_edge.dummy_v, relevant_edge.dummy_w, head_node.outgoing_weights[check_index]});
                    } else {
                        // two sided, at least I know my side and can make sure it is going to their head
                        head_node.outgoing_edgeindices[check_index] = their_endpoint;
                        if(endpoint < their_endpoint) {
                            // std::cout << "Pushing.. " << endpoint << " " << their_endpoint << " " << head_node.outgoing_weights[check_index] << std::endl;
                            ret_seq.push_back({endpoint, their_endpoint, head_node.outgoing_weights[check_index]});
                        }
                        if(this->simplified_tree[relevant_edge.dummy_v].owner == endpoint) {
                            relevant_edge.dummy_v = endpoint;
                        } else {
                            relevant_edge.dummy_w = endpoint;
                        }
                    }
                }
                if(head_node.new_conn[1]) {
                    const short check_index = 1;
                    head_node.new_conn[check_index] = false;
                    T v = endpoint;
                    T w = this->simplified_tree[head_node.outgoing_edgeindices[check_index]].owner;
                    T their_endpoint = w;
                    if(w < v)
                        std::swap(v, w);
                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                    // assert(relevant_edge.state); // todo
                    // std::cout << relevant_edge << std::endl;
                    if(relevant_edge.state == 1 || relevant_edge.state == 3) { // that side should not be changed
                        // nothing fancy to do
                        if(this->simplified_tree[relevant_edge.dummy_v].owner == endpoint) {
                            relevant_edge.dummy_v = endpoint;
                            auto& their_node = this->simplified_tree[relevant_edge.dummy_w];
                            short ind = -1;
                            for(short a = 0; a < 4; ++a) {
                                assert(a != 3); //
                                if(their_node.outgoing_edgeindices[a] != -1 && this->simplified_tree[their_node.outgoing_edgeindices[a]].owner == endpoint) {
                                    their_node.outgoing_edgeindices[a] = relevant_edge.dummy_v;
                                    break;
                                }
                            }   
                        }
                        else { 
                            relevant_edge.dummy_w = endpoint;
                            auto& their_node = this->simplified_tree[relevant_edge.dummy_v];
                            
                            for(short a = 0; a < 4; ++a) {
                                assert(a != 3); //
                                if(their_node.outgoing_edgeindices[a] != -1 && this->simplified_tree[their_node.outgoing_edgeindices[a]].owner == endpoint) {
                                    their_node.outgoing_edgeindices[a] = relevant_edge.dummy_w;
                                    break;
                                }
                            }
                        }   
                        // std::cout << "Pushing " << relevant_edge.dummy_v << " " << relevant_edge.dummy_w << " " << head_node.outgoing_weights[check_index] << std::endl;
                        ret_seq.push_back({relevant_edge.dummy_v, relevant_edge.dummy_w, head_node.outgoing_weights[check_index]});
                    } else {
                        // two sided, at least I know my side and can make sure it is going to their head
                        head_node.outgoing_edgeindices[check_index] = their_endpoint;
                        if(endpoint < their_endpoint) {
                            // std::cout << "Pushing " << endpoint << " " << their_endpoint << " " << head_node.outgoing_weights[check_index] << std::endl;
                            ret_seq.push_back({endpoint, their_endpoint, head_node.outgoing_weights[check_index]});
                        }
                        if(this->simplified_tree[relevant_edge.dummy_v].owner == endpoint) {
                            relevant_edge.dummy_v = endpoint;
                        } else {
                            relevant_edge.dummy_w = endpoint;
                        }
                    }
                }

                return ret_seq;
            }));

            parlay::parallel_for(0, endpoints.size(), [&] (T i) {
                // revert
                const T& endpoint = endpoints[i];
                auto& my_node = this->simplified_tree[endpoint];
                my_node.spill = -1; // reuse to keep track of extra index
                for(short a = 0; a < 2; a++) {
                    if(my_node.outgoing_edgeindices[a] != -1) {
                        T v = endpoint;
                        T w = this->simplified_tree[my_node.outgoing_edgeindices[a]].owner;
                        if(w < v)
                            std::swap(v, w);
                        auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                        relevant_edge.state = 0;
                    }
                }
            });

            // parlay::parallel_for(0, endpoints.size(), [&] (T i){ //todo remove this entire loop, it is unnecessary
            //     const T& endpoint = endpoints[i];
            //     T new_count = this->simplified_tree[endpoint].adjacency_count - this->new_connections[endpoint].size();
            //     if(new_count >= 3) {
            //         auto& node = this->simplified_tree[endpoint];
            //         assert(node.outgoing_edgeindices[0] != -1 && node.outgoing_edgeindices[1] != -1 && node.outgoing_edgeindices[2] != -1);   
            //     }
            // });

            // this->print_state();

            // this->print_state();

            auto tail_dels_p2 = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {
                parlay::sequence<std::pair<T,T>> ret_pairs;
                const T& endpoint = endpoints[i];
                auto& head_node = this->simplified_tree[endpoint];
                auto& tail_node = this->simplified_tree[head_node.tail_node_index];

                if(tail_node.dummy_index == endpoint) // I never spilled
                    return ret_pairs; 
                
                if(tail_node.outgoing_edgeindices[2] == -1 && tail_node.outgoing_edgeindices[1] == -1) {
                    // we need to remove the edge from of tail to tail, no tiebreak necessary
                    // ret_pairs.push_back({tail_node.dummy_index, tail_node.outgoing_edgeindices[0]});
                    T curr_tail_index = tail_node.dummy_index;
                    while(!this->simplified_tree[curr_tail_index].real && this->simplified_tree[curr_tail_index].outgoing_edgeindices[1] == -1) {
                        auto& tail_node = this->simplified_tree[curr_tail_index];
                        if(!tail_node.real)
                            ret_pairs.push_back({tail_node.dummy_index, tail_node.outgoing_edgeindices[0]});
                        // std::cout << tail_node << std::endl; 
                        curr_tail_index = tail_node.outgoing_edgeindices[0];
                        tail_node.outgoing_edgeindices[0] = tail_node.outgoing_edgeindices[1] = tail_node.outgoing_edgeindices[2] = -1; 
                        // assert(curr_tail_index != -1); // todo 
                    }
                    auto& tail_node = this->simplified_tree[curr_tail_index];
                    if(tail_node.real) {
                        tail_node.outgoing_edgeindices[2] = -1;
                        head_node.new_tail_node_index = curr_tail_index;
                        return ret_pairs; // just an empty tail
                    }

                    // assert(tail_node.marked_for_deletion == false); // todo
                    // whatever to my left is the real tail to consider
                    T payload = tail_node.outgoing_edgeindices[1];
                    // assert(payload != -1); // todo
                    T other_endpoint = this->simplified_tree[payload].owner;
                    // assert(other_endpoint != endpoint); // todo
                    // if(other_endpoint == -1) { // todo
                    //     std::cout << tail_node << std::endl;
                    // }
                    // assert(other_endpoint != -1); // todo
                    T v = endpoint;
                    T w = other_endpoint;
                    if(w < v)
                        std::swap(v, w);

                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                    relevant_edge.state.fetch_add(1);
                    // head_node.adds.push_back({v, w, tail_node.outgoing_weights[1]});
                    ret_pairs.push_back({curr_tail_index, tail_node.outgoing_edgeindices[0]});
                    head_node.new_tail_node_index = tail_node.outgoing_edgeindices[0];
                    this->simplified_tree[tail_node.outgoing_edgeindices[0]].outgoing_edgeindices[2] = other_endpoint;
                    this->simplified_tree[tail_node.outgoing_edgeindices[0]].outgoing_weights[2] = tail_node.outgoing_weights[1];
                    head_node.shortlisted = true;
                    head_node.spill = curr_tail_index;
                    tail_node.outgoing_edgeindices[0] = tail_node.outgoing_edgeindices[1] = tail_node.outgoing_edgeindices[2] = -1;
                    return ret_pairs;
                } else if(tail_node.outgoing_edgeindices[1] == -1) {
                    const short check_index = 2;
                    T curr_tail_index = tail_node.dummy_index;
                    T payload = tail_node.outgoing_edgeindices[check_index];
                    // assert(payload != -1); // todo
                    T other_endpoint = this->simplified_tree[payload].owner;
                    // assert(other_endpoint != endpoint); // todo
                    // assert(other_endpoint != -1);
                    D payload_data = tail_node.outgoing_weights[check_index];
                    // ret_pairs.push_back({tail_node.dummy_index, tail_node.outgoing_edgeindices[0]});

                    while(!this->simplified_tree[curr_tail_index].real && this->simplified_tree[curr_tail_index].outgoing_edgeindices[1] == -1) {
                        auto& tail_node = this->simplified_tree[curr_tail_index];
                        if(!tail_node.real)
                            ret_pairs.push_back({tail_node.dummy_index, tail_node.outgoing_edgeindices[0]});
                        curr_tail_index = tail_node.outgoing_edgeindices[0];
                        tail_node.outgoing_edgeindices[0] = tail_node.outgoing_edgeindices[1] = tail_node.outgoing_edgeindices[2] = -1;
                        // assert(curr_tail_index != -1); // todo 
                    }

                    // I found my new tail
                    head_node.new_tail_node_index = curr_tail_index;
                    T v = endpoint;
                    T w = other_endpoint;
                    if(w < v)
                        std::swap(v, w);
                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                    relevant_edge.state.fetch_add(1);

                    // head_node.adds.push_back({v, w, payload_data});

                    auto& tail_node = this->simplified_tree[curr_tail_index];
                    tail_node.outgoing_edgeindices[2] = other_endpoint;
                    tail_node.outgoing_weights[2] = payload_data;
                    head_node.shortlisted = true;

                } else if(tail_node.outgoing_edgeindices[2] == -1) {
                    const short check_index = 1;
                    T curr_tail_index = tail_node.dummy_index;
                    T payload = tail_node.outgoing_edgeindices[check_index];
                    // assert(payload != -1); // todo
                    T other_endpoint = this->simplified_tree[payload].owner;
                    // assert(other_endpoint != endpoint); // todo
                    // assert(other_endpoint != -1);
                    D payload_data = tail_node.outgoing_weights[check_index];
                    // ret_pairs.push_back({tail_node.dummy_index, tail_node.outgoing_edgeindices[0]});
                    bool first_time = true;

                    while(!this->simplified_tree[curr_tail_index].real && (this->simplified_tree[curr_tail_index].outgoing_edgeindices[1] == -1 || first_time)) {
                        first_time = false;
                        auto& tail_node = this->simplified_tree[curr_tail_index];
                        if(!tail_node.real)
                            ret_pairs.push_back({tail_node.dummy_index, tail_node.outgoing_edgeindices[0]});
                        curr_tail_index = tail_node.outgoing_edgeindices[0];
                        tail_node.outgoing_edgeindices[0] = tail_node.outgoing_edgeindices[1] = tail_node.outgoing_edgeindices[2] = -1;
                        // assert(curr_tail_index != -1); // todo 
                    }
                    // I found my new tail
                    head_node.new_tail_node_index = curr_tail_index;
                    T v = endpoint;
                    T w = other_endpoint;
                    if(w < v)
                        std::swap(v, w);
                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                    relevant_edge.state.fetch_add(1);
                    // head_node.adds.push_back({v, w, payload_data});
                    auto& tail_node = this->simplified_tree[curr_tail_index];
                    tail_node.outgoing_edgeindices[2] = other_endpoint;
                    tail_node.outgoing_weights[2] = payload_data;
                    head_node.shortlisted = true;
                }
                return ret_pairs;
            }));


            auto tail_dels_p3 = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {

                parlay::sequence<std::pair<T,T>> ret_pairs;

                const T& endpoint = endpoints[i];
                // does my new tails outgoing need to be updated
                auto& head_node = this->simplified_tree[endpoint];
                auto& tail_node = this->simplified_tree[head_node.new_tail_node_index];


                if(head_node.shortlisted == false)
                    return ret_pairs;

                T payload = tail_node.outgoing_edgeindices[2];
                D payload_val = tail_node.outgoing_weights[2];

                T v = endpoint;
                T w = payload;
                T other_endpoint = w;
                if(w < v)
                    std::swap(v, w);
                // assert(payload != -1); // todo

                // std::cout << "looking for " << v << " -- " << w << std::endl;
                auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                // assert(relevant_edge.state);
                if(relevant_edge.state == 0)
                    return ret_pairs;
                else if(relevant_edge.state == 1) {
                    ret_pairs.push_back({relevant_edge.dummy_v, relevant_edge.dummy_w});
                } else if(relevant_edge.state == 2) {
                    if(endpoint < other_endpoint) {
                        ret_pairs.push_back({relevant_edge.dummy_v, relevant_edge.dummy_w});
                    }
                }
                return ret_pairs;
            }));

            // checkDuplicates(tail_dels_p3); // todo

            // this->print_state();

            auto tail_adds_p2 = parlay::flatten(parlay::tabulate(endpoints.size(), [&] (T i) {
                parlay::sequence<wedge> ret_seq;

                const T& endpoint = endpoints[i];
                auto& head_node = this->simplified_tree[endpoint];
                auto& tail_node = this->simplified_tree[head_node.new_tail_node_index]; 
                if(head_node.shortlisted == false)
                    return ret_seq;
 
                T payload = tail_node.outgoing_edgeindices[2];
                D payload_val = tail_node.outgoing_weights[2];

                T v = endpoint;
                T w = payload;
                T other_endpoint = w;
                if(w < v)
                    std::swap(v, w);
                // assert(payload != -1); // todo
                
                auto& relevant_edge = this->find_basic_edge_in_ht(v, w);

                if(relevant_edge.state == 0)
                    return ret_seq;
                else if(relevant_edge.state == 1) {
                    char expected = 1;
                    char desired = 3;
                    bool success = relevant_edge.state.compare_exchange_strong(expected, desired);
                    if(!success) {
                        return ret_seq;
                    }
                    // We need to return an edge from my tail to wherever it is on their person
                    // since this edge is of type 1, it should be expected that their dummy node hasn't changed, only ours has
                    if(this->simplified_tree[relevant_edge.dummy_v].owner == endpoint) {
                        // oriented correctly
                        ret_seq.push_back({head_node.new_tail_node_index, relevant_edge.dummy_w, tail_node.outgoing_weights[2]});
                        for(short a = 0; a < 3; a++) {
                            // assert(a != 3); // todo
                            T& outgoing_index = this->simplified_tree[relevant_edge.dummy_w].outgoing_edgeindices[a];
                            if(outgoing_index == -1)
                                continue;
                            if(outgoing_index == relevant_edge.dummy_v){
                                outgoing_index = head_node.new_tail_node_index;
                                break;
                            }
                        }
                        tail_node.outgoing_edgeindices[2] = relevant_edge.dummy_w;
                        assert(tail_node.outgoing_edgeindices[2] == relevant_edge.dummy_w); // todo
                        relevant_edge.dummy_v = head_node.new_tail_node_index;
                    } else if(this->simplified_tree[relevant_edge.dummy_w].owner == endpoint) {
                        // flipped
                        ret_seq.push_back({head_node.new_tail_node_index, relevant_edge.dummy_v, tail_node.outgoing_weights[2]});
                        for(short a = 0; a < 3; a++) {
                            // assert(a != 3); // todo
                            T& outgoing_index = this->simplified_tree[relevant_edge.dummy_v].outgoing_edgeindices[a];
                            if(outgoing_index == -1)
                                continue;
                            if(outgoing_index == relevant_edge.dummy_w){
                                outgoing_index = head_node.new_tail_node_index;
                                break;
                            }
                        }
                        tail_node.outgoing_edgeindices[2] = relevant_edge.dummy_v;
                        // assert(tail_node.outgoing_edgeindices[2] == relevant_edge.dummy_v); // todo
                        relevant_edge.dummy_w = head_node.new_tail_node_index;
                    } else assert(false && "Shouldn't ever reach here");
                } else if(relevant_edge.state == 2) {
                    // I am guaranteed to be connected to their tail?
                    T their_tail_index = this->simplified_tree[other_endpoint].new_tail_node_index;                    
                    if(this->simplified_tree[relevant_edge.dummy_v].owner == endpoint)
                        relevant_edge.dummy_v = head_node.new_tail_node_index;
                    else if(this->simplified_tree[relevant_edge.dummy_w].owner == endpoint)
                        relevant_edge.dummy_w = head_node.new_tail_node_index;
                    else
                        assert(false && "Shouldn't ever reach here");
                    tail_node.outgoing_edgeindices[2] = their_tail_index;
                    if(endpoint < other_endpoint) {
                        ret_seq.push_back({head_node.new_tail_node_index ,their_tail_index, tail_node.outgoing_weights[2]});
                    }
                }
                return ret_seq;
            }));

            auto free_node_entries = parlay::filter(
                parlay::delayed_tabulate(all_dummy_endpoints.size(), [&] (T i) {
                    auto& node = this->simplified_tree[all_dummy_endpoints[i]];
                    if(node.real == false){
                        node.clear_data(this->identity);
                        return node.dummy_index;
                    }
                    return static_cast<T>(-1);
                }), [] (T i) {return i != -1;});

            auto extra_tail_frees = parlay::filter(parlay::delayed_tabulate(endpoints.size(), [&] (T i) {
                auto& node = this->simplified_tree[endpoints[i]];
                return node.spill;
            }), [] (T i) {return i != -1;});

            free_node_entries.append(extra_tail_frees);

            // for(auto& fr : free_node_entries)
            //     std::cout << fr << std::endl;

            // checkDuplicates(free_node_entries); // todo

            parlay::parallel_for(0, endpoints.size(), [&] (T i) {
                const T& endpoint = endpoints[i];
                auto& head_node = this->simplified_tree[endpoint];
                auto& tail_node = this->simplified_tree[head_node.new_tail_node_index];
                if(head_node.shortlisted) {
                    T v = endpoint;
                    T w = this->simplified_tree[tail_node.outgoing_edgeindices[2]].owner;
                    if(w < v)
                        std::swap(v, w);
                    auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
                    relevant_edge.state = 0;
                }
                head_node.tail_node_index = head_node.new_tail_node_index;
                head_node.spill = -1;
                head_node.marked_for_deletion = false;
                head_node.shortlisted = false;
                head_node.adjacency_count -= this->new_connections[endpoint].size(); 

                for(short a = 0; a < 2; a++) {
                    if(head_node.outgoing_edgeindices[a] == -1 && head_node.outgoing_edgeindices[a + 1] != -1) {
                        head_node.outgoing_edgeindices[a] = head_node.outgoing_edgeindices[a + 1];
                        head_node.outgoing_weights[a] = head_node.outgoing_weights[a+1]; 
                        head_node.outgoing_edgeindices[a + 1] = -1;
                    }
                }
                for(short a = 0; a < 2; a++) {
                    if(head_node.outgoing_edgeindices[a] == -1 && head_node.outgoing_edgeindices[a + 1] != -1) {
                        head_node.outgoing_edgeindices[a] = head_node.outgoing_edgeindices[a + 1];
                        head_node.outgoing_weights[a] = head_node.outgoing_weights[a+1]; 
                        head_node.outgoing_edgeindices[a + 1] = -1;
                    }
                }
                this->new_connections[endpoint].clear();
            });

            parlay::parallel_for(0, free_node_entries.size(), [&] (T i) {
                T free_entry = free_node_entries[i];
                T location = this->cursor - free_node_entries.size() + i;
                this->free_list[location] = free_entry;
            });

            this->cursor -= free_node_entries.size();
            // std::cout << "cursor decreased to " << cursor << std::endl;

            // std::cout << "Before size " << tail_dels_p2.size() << " " << linked_list_edges.size() << std::endl;
            // for(auto& dl : tail_dels_p2) {
            //     std::cout << red << dl.first << " " << dl.second << reset << std::endl;
            // }
            
            // for(auto& nl : tail_adds_p1) {
            //     std::cout << blue << std::get<0>(nl) << " " << std::get<1>(nl) << " " << std::get<2>(nl) << reset << std::endl;
            // }
            // remove_duplicates(tail_dels_p2, linked_list_edges);
            // for(auto& dl : tail_dels_p2) {
            //     std::cout << red << dl.first << " " << dl.second << reset << std::endl;
            // }
            
            // for(auto& nl : tail_adds_p1) {
            //     std::cout << blue << std::get<0>(nl) << " " << std::get<1>(nl) << " " << std::get<2>(nl) << reset << std::endl;
            // }
            
            // std::cout << "After size " << tail_dels_p2.size() << " " << linked_list_edges.size() << std::endl;

            ret_pairs.append(ternarized_edges);
            ret_pairs.append(del_iden_edges);
            ret_pairs.append(tail_dels_p1);
            ret_pairs.append(tail_dels_p2);
            ret_pairs.append(tail_dels_p3);
            
            
            ret_wedges.append(linked_list_edges);
            ret_wedges.append(tail_adds_p1); 
            ret_wedges.append(tail_adds_p2); 

            remove_duplicates(ret_pairs, ret_wedges);

            return {ret_pairs, ret_wedges};
        }

        std::pair<T,T> translate_edge(T v, T w) 
        {
            if(w < v)
                std::swap(v, w);
            auto& relevant_edge = this->find_basic_edge_in_ht(v, w);
            if(this->get_owner(relevant_edge.dummy_v) == this->get_owner(v))
                return {relevant_edge.dummy_v, relevant_edge.dummy_w};
            else
                return {relevant_edge.dummy_w, relevant_edge.dummy_v};

        }

        bool hasCycle(const parlay::sequence<wedge>& edges) {
            std::unordered_map<long, long> parent;
        
            auto find = [&](long x) {
                while (parent.count(x) && parent[x] != x)
                    x = parent[x] = parent[parent[x]];
                return parent.count(x) ? x : (parent[x] = x);
            };
        
            for (auto& e : edges) {
                long u, v; double w;
                std::tie(u, v, w) = e;
                long ru = find(u), rv = find(v);
                if (ru == rv) { //std::cout<<magenta<<"Warning, has cycle"<<reset<<std::endl;
                    return true;}
                parent[ru] = rv;
            }
            return false;
        }

        void verify_simple_tree(void)
        {
            
            
            parlay::parallel_for(0, this->num_vertices, [&] (T i) {
                auto& my_tern_node = this->simplified_tree[i];
                for(short I = 0; I < 3; I ++)
                {
                    T& neighbour_index = my_tern_node.outgoing_edgeindices[I];
                    if(neighbour_index == -1)
                        continue;

                    if(neighbour_index == my_tern_node.dummy_index) {
                        // std::cout << red << my_tern_node << reset << std::endl;
                        assert("self edge??" && false);
                    }
                    auto& neighbour_tern_node = this->simplified_tree[neighbour_index];
                    for(short j = 0; j < 3; j++)
                    {
                        if(neighbour_tern_node.outgoing_edgeindices[j] == my_tern_node.dummy_index)
                            return;
                    }
                    std::cout << my_tern_node  << std::endl;
                    std::cout <<  neighbour_tern_node  << std::endl;

                    // this->print_state();

                    assert("Neighbour doesn't have me?" && false);
                }
            });

            
            auto edge_list = parlay::flatten(parlay::delayed_tabulate(this->simplified_tree.size(),[&] (T I) {
                auto& node = this->simplified_tree[I];
                parlay::sequence<wedge> retwedge;
                for(int i = 0; i < 3; ++i) {
                    T me = I;
                    T them = node.outgoing_edgeindices[i];
                    if(me < them)
                        retwedge.push_back({me, them, this->identity});
                }
                return retwedge;
            }));

            parlay::parallel_for(0, this->hash_table_size, [&] (T i) {
                assert(this->key_edges[i].state == 0);
            });

            parlay::parallel_for(0, this->max_index, [&] (T i) {
                auto& head_node = this->simplified_tree[i];
                assert(head_node.real);

                T curr_index = head_node.dummy_index;

                while(curr_index != -1 && this->simplified_tree[curr_index].owner == head_node.dummy_index) {
                    auto& curr_node = this->simplified_tree[curr_index];
                    if(curr_node.real == false){
                        if(curr_node.outgoing_edgeindices[1] == -1)
                            std::cout << curr_node << std::endl;
                        assert(curr_node.outgoing_edgeindices[1] != -1);
                    }
                    if(curr_node.dummy_index == head_node.tail_node_index) {
                        if(head_node.adjacency_count > 3)
                            assert(curr_node.outgoing_edgeindices[2] != -1);
                    }

                    curr_index = curr_node.outgoing_edgeindices[2];
                }
            });

            parlay::parallel_for(0, this->hash_table_size, [&] (T i) {
                auto& relevant_edge = this->key_edges[i];
                if(this->key_valid[i] == false)
                    return;
                assert(relevant_edge.dummy_w != -1 && relevant_edge.dummy_v != -1);
                assert(this->get_owner(relevant_edge.dummy_v) == relevant_edge.v);
                assert(this->get_owner(relevant_edge.dummy_w) == relevant_edge.w);
            });

            // this->hasCycle(edge_list);

        }



        auto& get_simplified_tree(void) {
            return this->simplified_tree;
        }

};


#endif