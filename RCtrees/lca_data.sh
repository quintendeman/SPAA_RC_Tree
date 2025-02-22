NUM_THREADS=(1 2 4 8 16 32)
TRIES=(1 1 1 1 1)
make lca.out
for PNT in "${NUM_THREADS[@]}"; do
    for X in "${TRIES[@]}"; do
        LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./lca.out -forr -5 >> output_lca.csv
    done
done