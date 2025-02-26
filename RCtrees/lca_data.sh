#run with:
#bash lca_data.sh
date

#NUM_THREADS=(1 2 4 6 8 12 16 20 24 28 32) #could change to (1 2 4 8 16 32) if desired; the extra thread #s add some smoothness
NUM_THREADS=(1 2 4 8 16 32)

make lca.out
touch output_lca.csv

date

LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./lca.out -forr -4 -trials 10 -n 10000000 -k 10000000 -mean 20 -ln .1 >> output_lca.csv


LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=32 ./lca.out -forr -4 -trials 10 -n 10000000 -k 10000000 -mean 40 -ln .1 >> output_lca.csv

date

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./lca.out -forr -5 -n 10000000 -k 500000 -trials 10 -mean 20 -ln .1 >> output_lca.csv
done

for PNT in "${NUM_THREADS[@]}"; do
    LD_PRELOAD=$jemellocpath PARLAY_NUM_THREADS=$PNT ./lca.out -forr -5 -n 10000000 -k 500000 -trials 10 -mean 40 -ln .1 >> output_lca.csv
done

date
