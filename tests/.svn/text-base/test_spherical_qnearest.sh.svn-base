#!/bin/bash

for THREADS in 1 2 4 8 16
do
  for POINTS in 10 100 1000 10000
  do
    ./test_spherical_qnearest ${POINTS} ${THREADS}
    if [ $? -ne 0 ]; then
	echo "Spherical q-nearest neighbors failed for ${POINTS} points with ${THREADS} threads."
	exit 1
    fi
  done
done


