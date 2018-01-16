#!/bin/bash
upper_limit=1.0
for m in 1 2 3 4
do
	for p in {1..100}
	do
        for i in {1..100}
        do
        	mpirun -np $p ./main $m $upper_limit >> results${m}.txt
        done
    done
done 
echo Done
exit 0