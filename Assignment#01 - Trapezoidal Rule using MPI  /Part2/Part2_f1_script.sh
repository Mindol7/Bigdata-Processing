#!/bin/bash

for n in 1 2 4 8
do
    echo "Running with $n processes:"
    if [ $n -eq 8 ]; then
        mpiexec --oversubscribe -n $n MPI_Trapezoid_Part2_f1
    else
        mpiexec -n $n MPI_Trapezoid_Part2_f1
    fi
    echo "----------------------------------------"
done
