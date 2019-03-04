#!/bin/bash

{ echo "p, m, n, block, t"
for j in 4 8; do
	for i in 10 50 100 500 1000; do
		mpiexec --oversubscribe -n $j ./prog -b $i -n 4000 -m 2000
	done 
done } > a.out
