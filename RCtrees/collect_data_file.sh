filename=output_mst_from_file.csv

PARLAY_NUM_THREADS=1 ./testMSTfromFile.out sx-stackoverflow.txt 1 >> $filename
PARLAY_NUM_THREADS=1 ./testMSTfromFile.out sx-stackoverflow.txt 100 >> $filename
PARLAY_NUM_THREADS=1 ./testMSTfromFile.out sx-stackoverflow.txt 10000 >> $filename
PARLAY_NUM_THREADS=1 ./testMSTfromFile.out sx-stackoverflow.txt 100000 >> $filename
PARLAY_NUM_THREADS=1 ./testMSTfromFile.out sx-stackoverflow.txt 1000000 >> $filename

PARLAY_NUM_THREADS=8 ./testMSTfromFile.out sx-stackoverflow.txt 1 >> $filename
PARLAY_NUM_THREADS=8 ./testMSTfromFile.out sx-stackoverflow.txt 100 >> $filename
PARLAY_NUM_THREADS=8 ./testMSTfromFile.out sx-stackoverflow.txt 10000 >> $filename
PARLAY_NUM_THREADS=8 ./testMSTfromFile.out sx-stackoverflow.txt 100000 >> $filename
PARLAY_NUM_THREADS=8 ./testMSTfromFile.out sx-stackoverflow.txt 1000000 >> $filename

PARLAY_NUM_THREADS=16 ./testMSTfromFile.out sx-stackoverflow.txt 1 >> $filename
PARLAY_NUM_THREADS=16 ./testMSTfromFile.out sx-stackoverflow.txt 100 >> $filename
PARLAY_NUM_THREADS=16 ./testMSTfromFile.out sx-stackoverflow.txt 10000 >> $filename
PARLAY_NUM_THREADS=16 ./testMSTfromFile.out sx-stackoverflow.txt 100000 >> $filename
PARLAY_NUM_THREADS=16 ./testMSTfromFile.out sx-stackoverflow.txt 1000000 >> $filename

PARLAY_NUM_THREADS=32 ./testMSTfromFile.out sx-stackoverflow.txt 1 >> $filename
PARLAY_NUM_THREADS=32 ./testMSTfromFile.out sx-stackoverflow.txt 100 >> $filename
PARLAY_NUM_THREADS=32 ./testMSTfromFile.out sx-stackoverflow.txt 10000 >> $filename
PARLAY_NUM_THREADS=32 ./testMSTfromFile.out sx-stackoverflow.txt 100000 >> $filename
PARLAY_NUM_THREADS=32 ./testMSTfromFile.out sx-stackoverflow.txt 1000000 >> $filename

PARLAY_NUM_THREADS=48 ./testMSTfromFile.out sx-stackoverflow.txt 1 >> $filename
PARLAY_NUM_THREADS=48 ./testMSTfromFile.out sx-stackoverflow.txt 100 >> $filename
PARLAY_NUM_THREADS=48 ./testMSTfromFile.out sx-stackoverflow.txt 10000 >> $filename
PARLAY_NUM_THREADS=48 ./testMSTfromFile.out sx-stackoverflow.txt 100000 >> $filename
PARLAY_NUM_THREADS=48 ./testMSTfromFile.out sx-stackoverflow.txt 1000000 >> $filename