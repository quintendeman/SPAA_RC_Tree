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
template <typename graph> 
auto remove_higher_degree(graph& G, const int max_degree)
{
    
    int n = G.size();
    auto vertices = parlay::tabulate<int>(n, [&] (int i) {return n;});

    parlay::parallel_for(0, vertices.size(), [&] (int i) {
        if(G[i].size() > max_degree)
            G[i] = G[i].subseq(0, max_degree); // Does this result in gaps?
    });

    return;
}