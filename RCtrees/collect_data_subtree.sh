jemellocpath="/usr/lib/x86_64-linux-gnu/libjemalloc.so"


echo "num_threads,graph_size,ln,mean,dist,subtree_time,batched_subtree_time,time" > output_subtree.csv


CONST_GRAPH_SIZE=10000000 # graph size to use when running experiment with varying number of threads

NUM_THREADS=(48 32 16 8 1)

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./testBatchedSubtree.out --ln 0.1 --mean 20 --dist u --graph_size $CONST_GRAPH_SIZE  >> output_subtree.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./testBatchedSubtree.out --ln 0.9 --mean 40 --dist e --graph_size $CONST_GRAPH_SIZE  >> output_subtree.csv
done
