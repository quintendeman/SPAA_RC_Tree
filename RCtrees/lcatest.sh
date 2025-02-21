export PARLAY_NUM_THREADS=34
echo $PARLAY_NUM_THREADS > outfile13.txt
make lca.out >> outfile13.txt
nohup ./lca.out -forr -1 -pseed 72 >> outfile14.txt 2> outfile14err.txt #redirect stderr to stdout? 