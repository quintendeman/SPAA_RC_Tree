jemellocpath="/usr/lib/x86_64-linux-gnu/libjemalloc.so"

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=48 ./testBatchedSubtree.out --ln 0.1 --mean 20 --dist u --graph_size 100000000  >> output_big_subtree.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=48 ./testBatchedSubtree.out --ln 0.9 --mean 40 --dist e --graph_size 100000000  >> output_big_subtree.csv