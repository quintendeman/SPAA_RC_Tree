filename=output_mst_from_file_two.csv

PARLAY_NUM_THREADS=1 ./testMSTfromFile.out USA-road-d.USA.gr 1 >> $filename
PARLAY_NUM_THREADS=1 ./testMSTfromFile.out USA-road-d.USA.gr 100 >> $filename
PARLAY_NUM_THREADS=1 ./testMSTfromFile.out USA-road-d.USA.gr 10000 >> $filename
PARLAY_NUM_THREADS=1 ./testMSTfromFile.out USA-road-d.USA.gr 100000 >> $filename
PARLAY_NUM_THREADS=1 ./testMSTfromFile.out USA-road-d.USA.gr 1000000 >> $filename

PARLAY_NUM_THREADS=8 ./testMSTfromFile.out USA-road-d.USA.gr 1 >> $filename
PARLAY_NUM_THREADS=8 ./testMSTfromFile.out USA-road-d.USA.gr 100 >> $filename
PARLAY_NUM_THREADS=8 ./testMSTfromFile.out USA-road-d.USA.gr 10000 >> $filename
PARLAY_NUM_THREADS=8 ./testMSTfromFile.out USA-road-d.USA.gr 100000 >> $filename
PARLAY_NUM_THREADS=8 ./testMSTfromFile.out USA-road-d.USA.gr 1000000 >> $filename

PARLAY_NUM_THREADS=16 ./testMSTfromFile.out USA-road-d.USA.gr 1 >> $filename
PARLAY_NUM_THREADS=16 ./testMSTfromFile.out USA-road-d.USA.gr 100 >> $filename
PARLAY_NUM_THREADS=16 ./testMSTfromFile.out USA-road-d.USA.gr 10000 >> $filename
PARLAY_NUM_THREADS=16 ./testMSTfromFile.out USA-road-d.USA.gr 100000 >> $filename
PARLAY_NUM_THREADS=16 ./testMSTfromFile.out USA-road-d.USA.gr 1000000 >> $filename

PARLAY_NUM_THREADS=32 ./testMSTfromFile.out USA-road-d.USA.gr 1 >> $filename
PARLAY_NUM_THREADS=32 ./testMSTfromFile.out USA-road-d.USA.gr 100 >> $filename
PARLAY_NUM_THREADS=32 ./testMSTfromFile.out USA-road-d.USA.gr 10000 >> $filename
PARLAY_NUM_THREADS=32 ./testMSTfromFile.out USA-road-d.USA.gr 100000 >> $filename
PARLAY_NUM_THREADS=32 ./testMSTfromFile.out USA-road-d.USA.gr 1000000 >> $filename

PARLAY_NUM_THREADS=48 ./testMSTfromFile.out USA-road-d.USA.gr 1 >> $filename
PARLAY_NUM_THREADS=48 ./testMSTfromFile.out USA-road-d.USA.gr 100 >> $filename
PARLAY_NUM_THREADS=48 ./testMSTfromFile.out USA-road-d.USA.gr 10000 >> $filename
PARLAY_NUM_THREADS=48 ./testMSTfromFile.out USA-road-d.USA.gr 100000 >> $filename
PARLAY_NUM_THREADS=48 ./testMSTfromFile.out USA-road-d.USA.gr 1000000 >> $filename