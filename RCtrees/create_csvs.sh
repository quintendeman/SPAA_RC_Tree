#!/bin/bash

graphsizes=(100 1000 10000 100000 1000000 5000000 10000000)
querysizes=(100 1000 10000 100000 1000000 5000000 10000000)

rm n.csv
rm q.csv

max_n=10000000

for n in "${graphsizes[@]}"; do
    echo $(./RC --graph-size=$n --num-queries=100  --print-creation) >> n.csv

done

for q in "${querysizes[@]}"; do
    echo $(./RC --graph-size=$max_n --num-queries=$q  --print-query) >> q.csv

done