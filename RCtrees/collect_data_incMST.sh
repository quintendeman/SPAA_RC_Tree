jemellocpath="/usr/lib/x86_64-linux-gnu/libjemalloc.so"


# echo "num_threads,graph_size,ln,mean,dist,compressed_tree_gen_time,mst_gen_time,dynamic_insert_time,total_time" > output_incmst.csv

CONST_THREADS=48

CONST_GRAPH_SIZE=10000000 # graph size to use when running experiment with varying number of threads

NUM_THREADS=(48 32 16 8 1)

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS  ./testMST.out --ln 0.1 --mean 20 --dist e --graph_size $CONST_GRAPH_SIZE  >> output_incmst.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS  ./testMST.out --ln 0.1 --mean 40 --dist e --graph_size $CONST_GRAPH_SIZE  >> output_incmst.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS  ./testMST.out --ln 0.9 --mean 40 --dist e --graph_size $CONST_GRAPH_SIZE  >> output_incmst.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS  ./testMST.out --ln 0.9 --mean 40 --dist u --graph_size $CONST_GRAPH_SIZE  >> output_incmst.csv






for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT   ./testMST.out --ln 0.1 --mean 20 --dist e --graph_size $CONST_GRAPH_SIZE  --num-additions 1000000 >> output_incmst.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT   ./testMST.out --ln 0.1 --mean 40 --dist e --graph_size $CONST_GRAPH_SIZE  --num-additions 1000000 >> output_incmst.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT   ./testMST.out --ln 0.9 --mean 40 --dist e --graph_size $CONST_GRAPH_SIZE  --num-additions 1000000 >> output_incmst.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT   ./testMST.out --ln 0.9 --mean 40 --dist u --graph_size $CONST_GRAPH_SIZE  --num-additions 1000000 >> output_incmst.csv
done

