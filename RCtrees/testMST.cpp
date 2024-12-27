#include "../include/parlay/primitives.h"
#include "../include/parlay/sequence.h"
#include "../examples/kruskal.h"

using vertex = long;
using datatype = double;

using edge = std::tuple<vertex,vertex,double>;
using edgelist = parlay::sequence<edge>;

int main()
{
    edgelist EL;
    EL.push_back(edge(0,1,0.1));
    EL.push_back(edge(0,1,0.2));
    EL.push_back(edge(0,1,0.01));
    EL.push_back(edge(1,3,0.21));
    EL.push_back(edge(3,4,0.01));
    EL.push_back(edge(0,4,0.02));

    long totalvertices = 10;
    auto relevants = min_spanning_forest(EL, 10);

    std::cout << "Total relevant edges " << relevants.size() << std::endl;
    for(auto& relid : relevants)
    {
        std::cout << std::get<0>(EL[relid]) << " -- " << std::get<1>(EL[relid]) << ": " << std::get<2>(EL[relid]) << std::endl;
    }

    return 0;
}