
template <typename T, typename D, typename assocfunc>
D PathQueryHelper( cluster<T, D>* v,  cluster<T, D>* w, const D& defretval, assocfunc func)
{
    if(v == w)
        return defretval;

    bool use_v = false;

    cluster<T,D>* prev_boundary_v_l = nullptr;
    cluster<T,D>* prev_boundary_v_r = nullptr;
    D prev_val_till_v_l = defretval;
    D prev_val_till_v_r = defretval;
    
    cluster<T,D>* prev_boundary_w_l = nullptr;
    cluster<T,D>* prev_boundary_w_r = nullptr;
    D prev_val_till_w_l = defretval;
    D prev_val_till_w_r = defretval;


    while(true)
    {
        // Ascend using the lower-height cluster. If either has to ascend above root, they are not connected.
        if(v == nullptr && w == nullptr && prev_boundary_v_l == nullptr && prev_boundary_v_r == nullptr)
            return defretval;   
        else if(v == w)
            break;
        else if (v == nullptr && w != nullptr)
            use_v = false;
        else if (w == nullptr && v != nullptr)
            use_v = true;
        else if(v->get_height() < w->get_height())
            use_v = true;
        else
            use_v = false;
        
        if(use_v)
        {
            cluster<T,D>* &v_cluster = v;

            D lval_print,rval_print = defretval; // TODO remove

            // We are at the start probably can add a bunch of nullptr checks too TODO
            if (prev_val_till_v_l == defretval && prev_val_till_v_r == defretval) 
            {
                v_cluster->find_boundary_vertices(prev_boundary_v_l, prev_val_till_v_l, prev_boundary_v_r, prev_val_till_v_r, defretval);
                lval_print = prev_val_till_v_l;
                rval_print = prev_val_till_v_r;
            } 
            else 
            {
                D lval;
                D rval;
                cluster<T,D>* l;
                cluster<T,D>* r;

                v_cluster->find_boundary_vertices(l, lval, r, rval, defretval);

                // Covers both unary to unary and unary to binary case
                if (prev_boundary_v_l == v && prev_boundary_v_r == v) 
                {
                    lval_print = prev_val_till_w_l;
                    rval_print = prev_val_till_w_r;
                    
                    if(l != r) // ascending into unary
                    {
                        prev_val_till_v_l = func(prev_val_till_v_l, lval);
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                        
                    }
                    else
                    {
                        auto v_node_ptr = v_cluster->first_contracted_node;
                        for(auto& ptr : v_node_ptr->adjacents)
                        {
                            if(ptr != nullptr && ptr->state & (binary_cluster | base_edge))
                            {
                                prev_val_till_v_l = func(prev_val_till_v_l, ptr->cluster_ptr->data);
                                prev_val_till_v_r = func(prev_val_till_v_r, ptr->cluster_ptr->data);

                            }
                        }
                    }
                    prev_boundary_v_l = l;
                    prev_boundary_v_r = r;
                }
                // Covers binary ascending to unary/nullary
                else if (l == r) 
                {
                    if(prev_boundary_v_r == v)
                    {
                        prev_val_till_v_r = prev_val_till_v_l;
                    }
                    else
                        prev_val_till_v_l = prev_val_till_v_r;

                    prev_boundary_v_l = prev_boundary_v_r = l;
                }
                // Binary to binary
                else 
                {
                    if (prev_boundary_v_l == l) 
                    {
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    } 
                    else if (prev_boundary_v_r == l) 
                    {
                        std::swap(prev_boundary_v_r, prev_boundary_v_l);
                        std::swap(prev_val_till_v_l, prev_val_till_v_r);
                        
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    
                    } 
                    else if (prev_boundary_v_l == r) 
                    {
                        std::swap(r, l);
                        std::swap(rval, lval);    
                        
                        prev_boundary_v_r = r;
                        prev_val_till_v_r = func(prev_val_till_v_r, rval);
                    } 
                    else // prev_boundary_v_r == r 
                    {     
                        prev_boundary_v_l = l;
                        prev_val_till_v_l = func(prev_val_till_v_l, lval);
                    }
                }
                lval_print = lval;
                rval_print = rval;
            }
            // std::cout << bright_magenta << v->index << "[" << (short) v->get_height() << "]" << " ";
            // std::cout << prev_boundary_v_l->index << "(" << prev_val_till_v_l << ") " << prev_boundary_v_r->index << "(" << prev_val_till_v_r << ") ";
            // std::cout << "lval,rval=" << lval_print << "," << rval_print << " ";
            // auto& for_printing = v;
            // if(for_printing->adjacency.get_tail()->state & unary_cluster)
            //     std::cout << " unary ";
            // else if (for_printing->adjacency.get_tail()->state & binary_cluster)
            //     std::cout << " binary ";
            // else
            //     std::cout << " nullary ";
            // std::cout << reset << std::endl;
            v = v->get_parent();
        }
        else
        {
            cluster<T,D>* &w_cluster = w;
            D lval_print,rval_print = defretval;

            // We are at the start
            if (prev_val_till_w_l == defretval && prev_val_till_w_r == defretval) 
            {
                w_cluster->find_boundary_vertices(prev_boundary_w_l, prev_val_till_w_l, prev_boundary_w_r, prev_val_till_w_r, defretval);
                lval_print = prev_val_till_w_l;
                rval_print = prev_val_till_w_r;
            } 
            else 
            {
                D lval, rval;
                cluster<T,D>* l;
                cluster<T,D>* r;
                w_cluster->find_boundary_vertices(l, lval, r, rval, defretval);

                // Covers both unary to unary and unary to binary case
                if (prev_boundary_w_l == w && prev_boundary_w_r == w) 
                {

                    if(l != r)
                    {
                        prev_val_till_w_l = func(prev_val_till_w_l, lval);
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    }
                    else
                    {
                        auto w_node_ptr = w_cluster->first_contracted_node;
                        for(auto& ptr : w_node_ptr->adjacents)
                        {
                            if(ptr != nullptr && ptr->state & (binary_cluster | base_edge))
                            {
                                prev_val_till_w_l = func(prev_val_till_w_l, ptr->cluster_ptr->data);
                                prev_val_till_w_r = func(prev_val_till_w_r, ptr->cluster_ptr->data);

                            }
                        }
                    }
                    prev_boundary_w_l = l;
                    prev_boundary_w_r = r;
                }
                // Covers binary ascending to unary/nullary
                else if (l == r) 
                {
                    if(prev_boundary_w_r == w)
                    {
                        prev_val_till_w_r = prev_val_till_w_l;
                    }
                    else
                        prev_val_till_w_l = prev_val_till_w_r;

                    prev_boundary_w_l = prev_boundary_w_r = l;
                }
                // Binary
                else 
                {
                    if (prev_boundary_w_l == l) 
                    {
                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    } 
                    else if (prev_boundary_w_r == l) 
                    {
                        
                        std::swap(prev_boundary_w_l, prev_boundary_w_r);
                        std::swap(prev_val_till_w_l, prev_val_till_w_r);

                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);

                    } 
                    else if (prev_boundary_w_l == r) 
                    {
                        std::swap(l, r);
                        std::swap(lval, rval);

                        prev_boundary_w_r = r;
                        prev_val_till_w_r = func(prev_val_till_w_r, rval);
                    } 
                    else // prev_boundary_w_r == r 
                    { 
                        prev_boundary_w_l = l;
                        prev_val_till_w_l = func(prev_val_till_w_l, lval);
                    }
                }
                lval_print = lval;
                rval_print = rval;
            }
            // std::cout << bright_green << w->index << "[" << (short) w->get_height() << "]" << " ";
            // std::cout << prev_boundary_w_l->index << "(" << prev_val_till_w_l << ") " << prev_boundary_w_r->index << "(" << prev_val_till_w_r << ") ";
            // std::cout << "lval,rval=" << lval_print << "," << rval_print << " ";
            // auto& for_printing = w;
            // if(for_printing->adjacency.get_tail()->state & unary_cluster)
            //     std::cout << " unary ";
            // else if (for_printing->adjacency.get_tail()->state & binary_cluster)
            //     std::cout << " binary ";
            // else
            //     std::cout << " nullary ";
            // std::cout << reset << std::endl;
            w = w->get_parent();
        }
        
    }
    // std::cout << std::endl;

    if(prev_boundary_v_l == nullptr && prev_boundary_w_l == nullptr)
    {
        return defretval;
    }
    if(prev_boundary_v_l == nullptr)
    {
        if(w == prev_boundary_w_r)
            return prev_val_till_w_r;
        else
            return prev_val_till_w_l;   
    }
    if(prev_boundary_w_l == nullptr)
    {
        if(v == prev_boundary_v_l)
            return prev_val_till_v_l;
        else
            return prev_val_till_v_r;
    }
    
    D v_contrib, w_contrib;
    
    if(v == prev_boundary_v_l)
        v_contrib = prev_val_till_v_l;
    else
        v_contrib = prev_val_till_v_r;

    if(w == prev_boundary_w_l)
        w_contrib = prev_val_till_w_l;
    else
        w_contrib = prev_val_till_w_r;

    return func(v_contrib, w_contrib);

}

template<typename T, typename D, typename assocfunc>
D pathQuery(T v, T w, parlay::sequence<cluster<T,D>>& clusters, const D& identity, assocfunc func)
{
    return PathQueryHelper(&clusters[v], &clusters[w], identity, func);
}

// assumes you are using "+" as the assocfunc
template<typename T, typename D>
void testPathQueryValid(parlay::sequence<cluster<T,D>>& clusters, parlay::sequence<T>& parents, parlay::sequence<D>& weights, parlay::sequence<T>& random_perm_map, long base_size)
{
    parlay::parallel_for(0, 100000, [&] (T i) {
        T start_index = rand() % base_size;
        T index = start_index;
        T end_index;
        // doing sum!
        D manual_result = 0.0;
        while(parents[index] != index)
        {
            manual_result = manual_result + weights[index];
            index = parents[index];
        }
        end_index = index;

        D pathQueryResult = pathQuery(random_perm_map[start_index], random_perm_map[end_index], clusters, (D) 0.0, [] (D a, D b) {return a + b;});

        auto almost_equal = [] (D A, D B) -> bool {
        const double tolerance = 0.001;
            if(A == B)
                return true;
            if(A > B)
                if((A-B) < tolerance)
                    return true;
            if(B > A)
                if((B-A) < tolerance)
                    return true;
            return false;
        };

        if(!almost_equal(pathQueryResult, manual_result))
        {
            std::cout << "Path query from " << start_index << " to " << end_index << " invalid! " << std::endl;
            std::cout << pathQueryResult << " != " << manual_result << std::endl;
            assert("RC Path query and manual path query should be same" && false);
        }
       
    });

}
