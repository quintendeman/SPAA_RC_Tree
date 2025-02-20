jemellocpath="/usr/lib/x86_64-linux-gnu/libjemalloc.so"

echo "This script collects data for static and dynamic tree generation. Please set the LD_PRELOAD for jemelloc in this script, by default it is $jemellocpath"

echo "num_threads,graph_size,ln,mean,dist,static_gen_time,static_tern_time,dynamic_gen_time,dynamic_tern_time,randomize" > output.csv


CONST_GRAPH_SIZE=750 # graph size to use when running experiment with varying number of threads

graph_sizes=(500 600 700 800)
# graph_sizes=(1000 10000 100000 1000000 10000000 50000000 100000000)
NUM_THREADS=(1 2 4 8 16 32)

# run for ln 0.1 and mean 20 with exponential distribution
for graph_size in "${graph_sizes[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./testHashing.out --ln 0.1 --mean 20 --dist e --graph_size $graph_size >> output.csv
done

for graph_size in "${graph_sizes[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./testHashing.out --ln 0.1 --mean 20 --dist e --graph_size $graph_size --randomized >> output.csv
done

for graph_size in "${graph_sizes[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./testHashing.out --ln 0.1 --mean 40 --dist e --graph_size $graph_size >> output.csv
done

for graph_size in "${graph_sizes[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./testHashing.out --ln 0.9 --mean 40 --dist e --graph_size $graph_size >> output.csv
done

for graph_size in "${graph_sizes[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./testHashing.out --ln 0.9 --mean 40 --dist u --graph_size $graph_size >> output.csv
done



for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./testHashing.out --ln 0.1 --mean 20 --dist e --graph_size $CONST_GRAPH_SIZE >> output.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./testHashing.out --ln 0.1 --mean 20 --dist e --graph_size $CONST_GRAPH_SIZE  >> output.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./testHashing.out --ln 0.1 --mean 40 --dist e --graph_size $CONST_GRAPH_SIZE >> output.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./testHashing.out --ln 0.9 --mean 40 --dist e --graph_size $CONST_GRAPH_SIZE >> output.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./testHashing.out --ln 0.9 --mean 40 --dist u --graph_size $CONST_GRAPH_SIZE --randomized >> output.csv
done