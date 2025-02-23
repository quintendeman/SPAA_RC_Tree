export myfile=outfile19.txt
export myerrorfile=outfile19err.txt
export PARLAY_NUM_THREADS=50
echo $PARLAY_NUM_THREADS > $myfile
make lca.out >> $myfile
nohup ./lca.out -forr -1 -pseed 80 >> $myfile 2> $myerrorfile  #redirect stderr to stdout? 
