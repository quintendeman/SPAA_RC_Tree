jemellocpath="/usr/lib/x86_64-linux-gnu/libjemalloc.so"


echo "num_threads,graph_size,ln,mean,dist,batch_insert_size,time" > output_BatchInsert.csv

CONST_THREADS=48

CONST_GRAPH_SIZE=10000000 # graph size to use when running experiment with varying number of threads

NUM_THREADS=(48 32 16 8 1)

# run for ln 0.1 and mean 20 with exponential distribution
LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS   ./testBatchInsertion.out --ln 0.1 --mean 20 --dist e --graph_size $CONST_GRAPH_SIZE >> output_BatchInsert.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS   ./testBatchInsertion.out --ln 0.1 --mean 20 --dist u --graph_size $CONST_GRAPH_SIZE >> output_BatchInsert.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS   ./testBatchInsertion.out --ln 0.1 --mean 2 --dist e --graph_size $CONST_GRAPH_SIZE >> output_BatchInsert.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS   ./testBatchInsertion.out --ln 0.9 --mean 2 --dist e --graph_size $CONST_GRAPH_SIZE >> output_BatchInsert.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$CONST_THREADS   ./testBatchInsertion.out --ln 0.9 --mean 1.1 --dist u --graph_size $CONST_GRAPH_SIZE >> output_BatchInsert.csv


for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT   ./testBatchInsertion.out --ln 0.1 --mean 10 --dist e --graph_size $CONST_GRAPH_SIZE --num-insertions 611000  >> output_BatchInsert.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT   ./testBatchInsertion.out --ln 0.1 --mean 20 --dist e --graph_size $CONST_GRAPH_SIZE --num-insertions 1021000 >> output_BatchInsert.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT   ./testBatchInsertion.out --ln 0.9 --mean 20 --dist e --graph_size $CONST_GRAPH_SIZE --num-insertions 1021000  >> output_BatchInsert.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT   ./testBatchInsertion.out --ln 0.9 --mean 10 --dist u --graph_size $CONST_GRAPH_SIZE --num-insertions 611000  >> output_BatchInsert.csv
done

