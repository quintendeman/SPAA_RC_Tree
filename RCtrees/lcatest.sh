export myfile=outfile17.txt
export myerrorfile=outfile17err.txt
export PARLAY_NUM_THREADS=33
echo $PARLAY_NUM_THREADS > $myfile
make lca.out >> $myfile
nohup ./lca.out -forr -1 -pseed 74 >> $myfile 2> $myerrorfile  #redirect stderr to stdout? 
