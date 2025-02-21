NUM_THREADS=(1 2 4 8 16 32)

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./lca.out -forr -6 >> output_lca.csv
done