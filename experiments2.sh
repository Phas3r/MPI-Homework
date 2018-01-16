#!/bin/bash
p=8
for m in 1 2 3 4
do
	for upper_limit in {1..1000}
	do
        for i in {1..100}
        do
        	mpirun -np $p ./main $m $upper_limit >> results2${m}.txt
        done
    done
done 
echo Done
exit 0