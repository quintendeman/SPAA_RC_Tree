NUM_THREADS=(1 2 4 6 8 12 16 20 24 28 32)
TRIES=(1 1 1) #note, perhaps increase to more tries
make lca.out
touch output_lca.csv

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./lca.out -forr -4 >> output_lca.csv

for PNT in "${NUM_THREADS[@]}"; do
    for X in "${TRIES[@]}"; do
        LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./lca.out -forr -5 -n 50000000 -k 500000 -trials 3 >> output_lca.csv
    done
done
#Tree structure really affects runtime? q
REV_ORDER=(32 28 24 20 16 12 8 6 4 2 1)
for PNT in "${REV_ORDER[@]}"; do
    for X in "${TRIES[@]}"; do
        LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./lca.out -forr -5 -n 50000000 -k 500000 -trials 3 >> output_lca.csv
    done
done

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./lca.out -forr -4 >> output_lca.csv
