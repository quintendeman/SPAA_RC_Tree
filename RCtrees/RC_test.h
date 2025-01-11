//testing functions for RC trees
//Pulling from RC.cpp because LCA code will use too 

#ifndef RCTEST_H
#define RCTEST_H



const double default_epsilon = 0.001;


static bool isNearlyEqual(double a, double b, double epsilon = default_epsilon) {
    return std::abs(a - b) < epsilon;
}

//given an RC tree, makes sure that various invariants hold
//the parent of a child of v is v
//the child of a parent of v is v
//for addition) makes sure that the value of a parent is the sum of the values of the children
//specific type of subtree query works
//if print is true, print out debugging info, otherwise don't
template<typename vertex, typename datatype> 
void test_rc_valid(const parlay::sequence<vertex>& parents, parlay::sequence<cluster<vertex, datatype>>& clusters, bool print=true)
{

    check_parents_children(clusters);
    check_children_values(clusters);

    parlay::sequence<std::tuple<vertex, vertex, datatype>> weighted_edges = parlay::tabulate(parents.size(), [&] (vertex i) {
        return std::tuple<vertex, vertex, datatype>(i, parents[i], (datatype) ( 1.0));
    });
    //only keep nondegenerate edges (where the two endpoints are different)
    weighted_edges = parlay::filter(weighted_edges, [&] (auto edge) {
        return std::get<0>(edge) != std::get<1>(edge);
    });
    
    //change all of the weights in the tree to 1
    batchModifyEdgeWeights(weighted_edges, [] (datatype a, datatype b) {
        return a + b;
    }, clusters);

    //TOD2* is this a good subtree query test? Because are we given that the direction giver will be an RC tree parent of the random cluster? 

    auto random_index = rand() % parents.size();
    if (print)
    std::cout << "Root: " << random_index << " Parent: " << parents[random_index] << std::endl;

    auto manual_sum_val = manual_subtree_sum(&clusters[random_index], &clusters[parents[random_index]], clusters);
    if (print)
    std::cout << "Manual subtree sum: " << red << manual_sum_val << reset << std::endl;

    
    auto sub_ret_val = subtree_query(&clusters[random_index], &clusters[parents[random_index]], (datatype) 0.0, [] (datatype a, datatype b) {
        return a+b;
    });

    if (!isNearlyEqual(manual_sum_val,sub_ret_val)) {
        std::cout << "manual and RC tree subtree sums different, aborting" << std::endl;
        exit(17);
    }
    if (print) {
        std::cout << "Subtree query returns: " << bold << red << sub_ret_val << reset << std::endl;
        std::cout << std::endl;

    }
   

    return;
}

#endif