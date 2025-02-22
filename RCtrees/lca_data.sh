NUM_THREADS=(1 2 4 6 8 12 16 20 24 28 32)
TRIES=(1) #note, perhaps increase to more tries
make lca.out
for PNT in "${NUM_THREADS[@]}"; do
    for X in "${TRIES[@]}"; do
        LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./lca.out -forr -5 >> output_lca.csv
    done
done