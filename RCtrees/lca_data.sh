date

NUM_THREADS=(1 2 4 6 8 12 16 20 24 28 32)
make lca.out
touch output_lca.csv

date

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./lca.out -forr -4 -trials 10 -n 10000000 -k 10000000 >> output_lca.csv

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./lca.out -forr -5 -n 10000000 -k 500000 -trials 10 >> output_lca.csv
done

date